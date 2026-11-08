# -*- coding: utf-8 -*-
# Rebuilds the two panels of S5, for which no code was kept, from the original in one pass
#
# 왜 다시 만드나 - 앞서 pdf 를 여러 번 덧칠(redact + insert)했더니
# ⒜ 돌아간 loci frac 6 개 중 4 개가 지워졌고(원본 6, 내 판 2)
# ⒝ β (hardcall) 이 , (har + ll) 로 토막 났고
# ⒞ 눈금 숫자가 제자리를 벗어났다(오른쪽 정렬을 왼쪽 정렬로 다시 적었다).
# 고치는 길 - svg 를 정본으로 삼는다. 실측으로 svg 좌표 = pdf 좌표(둘 다 1:1, 왼위 기준)임을
# 확인했다(loci frac 6 개의 가운데 y 가 109.6, 222.3, … 한 자도 안 다름).
# ⇒ pdf 의 글자층(BT…ET)을 통째로 걷어내고 svg 가 말하는 자리에 다시 적는다.
#
# 하는 일
# ⑴ axes_33(오른쪽 아래 긴 그래프) 을 짝 맞춰 잘라낸다 - 앞 판은 닫는 </g> 를 같이 잘라
# XML 이 깨졌다.
# ⑵ 그래프 안 글(high mode …, low mode …, 그 그래프의 축 이름)을 뺀다.
# ⑶ 눈금은 축마다 처음과 끝만 남긴다(값 목록이 아니라 차례로 - 축마다 범위가 다르다).
# ⑷ 남은 글자를 T px 로(= 종이에서 3.5 pt).
# ⑸ A 패널 : 줄 이름(spacer-1 + …)을 왼쪽에서 오른쪽으로 옮기고
# 빈 왼쪽 여백만큼 통째로 민다 - 안 그러면 틀맞춤이 쪽 밖으로 끌고 나간다.
import html
import io
import os
import re
import shutil

import pymupdf


def _need(k):
    """The deposited copy keeps no default: it stops rather than read the wrong place."""
    v = os.environ.get(k)
    if not v:
        raise SystemExit("environment variable %s is not set (the "
                         "deposited copy has no default; see README.md,"
                         " section 'Environment variables')" % k)
    return v


# 기본 14 글꼴(Helvetica)에는 β 가 없다 - 앞 판에서 β (hardcall) 이 , (hardcall) 로
# 조용히 바뀌었다(내 기록 en-figure-arial-missing-glyphs 와 같은 함정).
# ⇒ 진짜 Arial 을 심어서 적는다. 폭도 Arial 로 재야 자리가 맞는다(svg 가 Arial 이다).
AR = os.path.join("C:" + os.sep, "Windows", "Fonts", "arial.ttf")
ARB = os.path.join("C:" + os.sep, "Windows", "Fonts", "arialbd.ttf")
FONT = {"ari": pymupdf.Font(fontfile=AR), "arib": pymupdf.Font(fontfile=ARB)}

SP = os.path.dirname(os.path.abspath(__file__))
# 원본 패널 넷을 저장소 안(figS5/inputs_panels)으로 들였다.
# search order : the S42_SRC environment value, then inside the repository, then the original.
_R = os.path.dirname(os.path.dirname(SP))
SRC = (os.environ.get("S42_SRC")
       or os.path.join(_R, "figS5", "inputs_panels"))
if not os.path.isdir(SRC):
    SRC = os.path.join(_need("RDNA_ARRAY_METH_DIR"), "03_figures",
                       "00_final_set_2026-10-16", "panels")
# PANEL_OUT 을 받는다 - 앞서는 폴더가 박혀 있어
# Run from the repository copy alone, the assembler does not find this panel.
DST = os.path.join(SP, os.environ.get("PANEL_OUT", "panelsplit_spec"), "_panels")
NL = chr(10)
os.makedirs(DST, exist_ok=True)
D1 = "FIG_D1_CODING32__panel_a"
D5 = "FIG_D5_a4_persite_tone_hardcall32"
T = float(os.environ.get("S42_T", "10.46"))      # 패널 안 글자 크기(px)
TEXT = re.compile(r"<text ([^>]*)>(.*?)</text>", re.S)
GRE = re.compile(r"<g\b|</g>")


def _rgb(h):
    """`#rrggbb` -> (r, g, b) 0~1"""
    return tuple(int(h[1 + 2 * i:3 + 2 * i], 16) / 255.0 for i in range(3))


def cut_block(s, tag):
    """`<g id="tag">` 를 짝을 세어 그 묶음만 잘라낸다"""
    i = s.find('<g id="%s">' % tag)
    if i < 0:
        raise SystemExit("못 찾음 : " + tag)
    j, depth = i, 0
    while True:
        m = GRE.search(s, j)
        if not m:
            raise SystemExit("짝이 안 맞는다 : " + tag)
        depth += 1 if m.group(0) == "<g" else -1
        j = m.end()
        if depth == 0:
            break
    return s[:i] + s[j:], j - i


def attr(a, name, default=None):
    m = re.search(r'%s="([^"]*)"' % name, a)
    return m.group(1) if m else default


def style(a, key, default=None):
    m = re.search(r"%s:\s*([^;\"]+)" % key, a)
    return m.group(1).strip() if m else default


def texts(s):
    out = []
    for m in TEXT.finditer(s):
        a, t = m.group(1), html.unescape(m.group(2)).strip()
        if not t:
            continue
        out.append(dict(span=m.span(), raw=m.group(0), x=float(attr(a, "x", "0")),
                        y=float(attr(a, "y", "0")),
                        size=float(style(a, "font-size").replace("px", "")),
                        anchor=style(a, "text-anchor", "start"),
                        bold="700" in (style(a, "font-weight", "") or ""),
                        rot="rotate(-90" in (attr(a, "transform", "") or ""),
                        text=t))
    return out


# 글귀 바꾸기
# ⒜ B 의 x 축 이름 : β (hardcall) -> Unit methylation β (하드콜이라는 사정은 legend 로)
# ⒝ 가운데점(,)은 그림 규칙상 금지 - 쉼표로
WORD = ((chr(0x3b2) + " (hardcall)", "Unit methylation " + chr(0x3b2)),
        (" " + chr(0xb7) + " ", ", "))


def _word(t):
    for a_, b_ in WORD:
        t = t.replace(a_, b_)
    return t


def _resize(m):
    r = re.sub(r"font-size: [\d.]+px", "font-size: %gpx" % T, m.group(0), count=1)
    return r[:r.index(">") + 1] + _word(r[r.index(">") + 1:r.rindex("<")])         + r[r.rindex("<"):]


# ──────────────────────────────── B 패널 (D5)
s = io.open(os.path.join(SRC, D5 + ".svg"), encoding="utf8").read()
n_before = s.count("<text")
s, ncut = cut_block(s, "axes_33")
print("D5 : axes_33 %d 자 잘라냄" % ncut)

keep = set()
ax = [m.start() for m in re.finditer(r'<g id="axes_\d+">', s)] + [len(s)]
for a0, a1 in zip(ax, ax[1:]):
    for kind in ("xtick", "ytick"):
        pos = [m.start() for m in re.finditer(r'<g id="%s_\d+">' % kind, s[a0:a1])]
        if not pos:
            continue
        for p in (pos[0], pos[-1]):
            tm = TEXT.search(s, a0 + p)
            if tm and tm.start() < a1:
                keep.add(tm.span())

drop = []
for t in texts(s):
    if t["text"].startswith(("high mode", "low mode")):
        drop.append(t)
    elif abs(t["size"] - 7.5) < 0.1:          # 뺀 그래프의 축 이름
        drop.append(t)
    elif abs(t["size"] - 6.5) < 0.2 and t["span"] not in keep:
        drop.append(t)
out = s
for t in sorted(drop, key=lambda z: -z["span"][0]):
    out = out[:t["span"][0]] + out[t["span"][1]:]
out = TEXT.sub(_resize, out)
io.open(os.path.join(DST, D5 + ".svg"), "w", encoding="utf8", newline=NL).write(out)
print("D5 : 글 %d -> %d (솎음 %d) · 크기 %.2f px" %
      (n_before, out.count("<text"), len(drop), T))

# ──────────────────────────────── A 패널 (D1)
s1 = io.open(os.path.join(SRC, D1 + ".svg"), encoding="utf8").read()
W1 = float(re.search(r'<svg[^>]*width="([\d.]+)pt"', s1).group(1))


def poly(m):
    xs, ys = [float(m.group(1))], [float(m.group(2))]
    for a, b in re.findall(r"L\s*([-\d.]+)\s+([-\d.]+)", m.group(3)):
        xs.append(float(a))
        ys.append(float(b))
    return xs, ys


# x 창을 유닛 끝에서 끊는다
# 실측 : 막대는 2240 bp(x 853.9)에서 끝나는데 축 틀은 2600 bp(x 959.2)까지 있다.
# 유닛보다 넓은 창은 없는 자리를 있는 것처럼 보이게 한다 ⇒ 틀,눈금,축 이름을 함께 옮긴다.
i1 = s1.find('<g id="axes_1">')
fx, fy = poly(re.search(r'<path d="M\s*([-\d.]+)\s+([-\d.]+)\s*(L[^"]*?)"',
                        s1[i1:i1 + 4000], re.S))
XR0, FY0, FY1 = max(fx), min(fy), max(fy)
XL0 = min(fx)
bars = []
for m in re.finditer(r'<path d="M\s*([-\d.]+)\s+([-\d.]+)\s*(L[^"]*?)z', s1, re.S):
    bx, _by = poly(m)
    if max(bx) < XR0 - 1:
        bars.append(max(bx))
XCUT = max(bars)
print("D1 : x 창 %.1f -> %.1f (데이터 %.0f -> %.0f bp)"
      % (XR0, XCUT, (XR0 - XL0) / 0.29340, (XCUT - XL0) / 0.29340))
# ⒜ 틀,배경의 오른쪽 선을 당긴다
blk = s1[i1:]
blk = re.sub(r"(?<![\d.])%s(?![\d.])" % re.escape("%g" % XR0), "%.4f" % XCUT, blk)
s1 = s1[:i1] + blk
# ⒝ 잘린 바깥의 눈금 묶음을 뺀다 (눈금 선과 숫자가 한 묶음이다)
while True:
    gone = None
    for m in re.finditer(r'<g id="(xtick_\d+)">', s1):
        tm = TEXT.search(s1, m.start())
        if tm and float(attr(tm.group(1), "x", "0")) > XCUT + 0.5:
            gone = m.group(1)
            break
    if not gone:
        break
    s1, _c = cut_block(s1, gone)
    print("D1 : 창 밖 눈금 %s 뺌" % gone)
# ⒞ 2026-10-20 : x 축 이름(array unit coordinate (v11 canonical; …))과
# 패널 제목(a) which part of the array unit …)을 뺀다.
# 까닭 : v11 은 학술 용어가 아니다(우리 내부 판 이름) - legend 에도 적지 않는다.
# 대신 어느 자리가 어느 구획인지 Fig 4 F 처럼 축 아래 색 띠로 적는다(아래 ⒠).
# 뒤에서 앞으로 지운다 - 앞에서 지우면 뒤 자리값이 밀려 엉뚱한 데를 자른다
# (실제로 그렇게 해서 </svg> 까지 잘려 나갔다)
for _t in sorted([t for t in texts(s1)
                  if "coordinate" in t["text"] or t["text"].startswith("a) which part")],
                 key=lambda z: -z["span"][0]):
    s1 = s1[:_t["span"][0]] + s1[_t["span"][1]:]
    print("D1 : 뺀 글 %r" % _t["text"][:46])

# ⒟ 막대 색을 논문 구획 팔레트로 (색을 아래,위와 맞출 것)
# 이 패널만 제 팔레트(주황~파랑 층)를 쓰고 있었다 ⇒ Fig 4 E,F,Fig 5 와 같은 색으로.
BARCOL = [("#f6d5a8", "#9ecae1", "spacer-1"), ("#cfe3d4", "#a6dbb4", "GA-rich"),
          ("#bcd6e8", "#d9d9d9", "spacer-2"), ("#9ec5e0", "#f2b705", "5'-flank"),
          ("#2a6faf", "#9e9e9e", "5S gene")]
for _o, _n, _nm in BARCOL:
    s1 = s1.replace("fill: %s" % _o, "fill: %s" % _n)

# ⒠ 구획 색 띠 + 이름 - Fig 4 F 와 같은 규칙(넓은 것은 가로, 좁은 것은 아래 줄 + 이음선)
ELN = [("AluY", 1, 318), ("spacer-1", 319, 1178), ("GA-rich", 1179, 1551),
       ("spacer-2", 1552, 2031), ("5'-flank", 2032, 2121), ("5S gene", 2122, 2240)]
ELC = {"AluY": "#c7b2de", "spacer-1": "#9ecae1", "GA-rich": "#a6dbb4",
       "spacer-2": "#d9d9d9", "5'-flank": "#f2b705", "5S gene": "#9e9e9e"}
BPX = (XCUT - XL0) / 2240.0                    # 1 bp 가 몇 pt 인가
BY0, BH = FY1 + 19.0, 9.0                      # 띠의 위 y, 높이
LINE1, LINE2 = BY0 + BH + T * 1.05, BY0 + BH + T * 2.35
BAND, BTX = [], []
_last = {LINE1: -1e9, LINE2: -1e9}
for _nm, _a, _b in ELN:
    x0, x1 = XL0 + (_a - 1) * BPX, XL0 + _b * BPX
    BAND.append((x0, BY0, x1, BY0 + BH, ELC[_nm]))
    _w = FONT["ari"].text_length(_nm, T)
    _cx = (x0 + x1) / 2.0
    # 이름이 띠보다 넓으면 이웃과 부딪친다 ⇒ 두 줄로 번갈아 놓고 이음선을 긋는다
    _y = LINE1 if _cx - _w / 2.0 > _last[LINE1] + T * 0.4 else LINE2
    if _cx - _w / 2.0 <= _last[_y] + T * 0.4:
        _y = LINE2 if _y == LINE1 else LINE1
    _last[_y] = _cx + _w / 2.0
    BTX.append((_nm, _cx, _y, _y != LINE1))
NEWH = LINE2 + T * 0.45
tt = texts(s1)
ann = [t for t in tt if "sites" in t["text"] and "loci" in t["text"]]
lab = [t for t in tt if t["anchor"] == "end"]
assert len(ann) == 5 and len(lab) == 5, (len(ann), len(lab))
annw = max(FONT["ari"].text_length(t["text"], T) for t in ann)
NEWX = ann[0]["x"] + annw + T * 1.2
new1 = s1
for t in sorted(lab, key=lambda z: -z["span"][0]):
    r = t["raw"]
    r = re.sub(r"text-anchor: \w+", "text-anchor: start", r)
    r = re.sub(r'x="[-\d.]+"', 'x="%.3f"' % NEWX, r, count=1)
    r = re.sub(r"rotate\(-0 [-\d.]+ ", "rotate(-0 %.3f " % NEWX, r)
    new1 = new1[:t["span"][0]] + r + new1[t["span"][1]:]
new1 = TEXT.sub(_resize, new1)
RIGHT = NEWX + max(FONT["ari"].text_length(t["text"], T) for t in lab) + 2.0
XMIN = min(t["x"] for t in texts(new1) if t["anchor"] == "middle")
DX = -max(0.0, XMIN - T * 1.6)
NEWW = RIGHT + DX + 2.0
# 구획 띠와 그 이름을 svg 에 넣는다 - 글자는 pdf 를 다시 깔 때 이 svg 에서 읽힌다
_band = []
for x0, y0, x1, y1, c in BAND:
    _band.append('   <rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" '
                 'style="fill: %s; stroke: none"/>' % (x0, y0, x1 - x0, y1 - y0, c))
for _nm, _cx, _y, _low in BTX:
    if _low:      # 아랫줄 이름은 띠에서 이음선을 긋는다
        _band.append('   <path d="M %.3f %.3f L %.3f %.3f" style="stroke: #9aa3ad; '
                     'stroke-width: 0.5; fill: none"/>'
                     % (_cx, BY0 + BH, _cx, _y - T * 0.95))
    _band.append('   <text style="font-size: %gpx; font-family: \'Arial\', '
                 '\'DejaVu Sans\', sans-serif; text-anchor: middle" x="%.3f" y="%.3f" '
                 'transform="rotate(-0 %.3f %.3f)">%s</text>'
                 % (T, _cx, _y, _cx, _y, html.escape(_nm)))
# 띠는 figure_1 밖에 붙으므로 통째 옮김(DX)을 제 transform 으로 직접 받는다
new1 = new1.replace("</svg>", NL.join(
    ['  <g id="element_band" transform="translate(%.3f,0)">' % DX] + _band + ["  </g>"])
    + NL + "</svg>", 1)
new1 = re.sub(r'(<svg[^>]*width=")[\d.]+(pt" height=")[\d.]+(pt")',
              r"\g<1>%.2f\g<2>%.2f\g<3>" % (NEWW, NEWH), new1, count=1)
new1 = re.sub(r'viewBox="0 0 [\d.]+ [\d.]+"',
              'viewBox="0 0 %.2f %.2f"' % (NEWW, NEWH), new1, count=1)
new1 = new1.replace('<g id="figure_1">',
                    '<g id="figure_1" transform="translate(%.3f,0)">' % DX, 1)
io.open(os.path.join(DST, D1 + ".svg"), "w", encoding="utf8", newline=NL).write(new1)
print("D1 : 줄 이름 x %.1f -> %.1f · 통째로 %.1f 밀고 폭 %.1f -> %.1f"
      % (lab[0]["x"], NEWX, DX, W1, NEWW))


# ──────────────────────────────── pdf : 글자층을 새로 깐다
# axes_33(오른쪽 아래 긴 그래프)의 그림(선,점)은 글자층이 아니라서 안 걷힌다
# ⇒ 그 자리를 지운다. 틀 실측 x 596.3~1104.8, y 638.1~707.6 (+ 눈금 여유)
CUT33 = (585.0, 625.0, 1116.0, 748.8)


def rebuild_pdf(name, dx=0.0, neww=None, newh=None, cut=None, box=None,
                band=None, leads=(), recolor=None):
    shutil.copy2(os.path.join(SRC, name + ".pdf"), os.path.join(DST, name + ".pdf"))
    q = os.path.join(DST, name + ".pdf")
    d = pymupdf.open(q)
    pg = d[0]
    oldh = pg.rect.height
    # 쪽을 아래로 넓힐 때는 내용을 그만큼 위로 올려야 한다 -
    # pdf 좌표는 아래 왼쪽이 원점이라, 안 올리면 새 자리가 위에 생긴다.
    dy = (newh - oldh) if newh else 0.0
    pg.wrap_contents()
    cont = pg.read_contents()
    n0 = len(re.findall(rb"BT", cont))
    cont = re.sub(rb"BT.*?ET", b"", cont, flags=re.S)       # 글자층을 걷어낸다
    # 막대 색 바꾸기는 svg 뿐 아니라 pdf 내용흐름의 r g b rg 도 고쳐야 한다
    # (png 는 pdf 에서 굽는다 - svg 만 고치면 두 꼴이 어긋난다)
    for _o, _n, _nm in (recolor or ()):
        _a = (" ".join("%.10g" % (int(_o[1 + 2 * i:3 + 2 * i], 16) / 255.0)
                       for i in range(3))).encode()
        _b = (" ".join("%.10g" % (int(_n[1 + 2 * i:3 + 2 * i], 16) / 255.0)
                       for i in range(3))).encode()
        if _a + b" rg" in cont:
            cont = cont.replace(_a + b" rg", _b + b" rg")
        else:
            print(" colour %s(%s) not found in the pdf" % (_nm, _o))
    cont = ("q 1 0 0 1 %.4f %.4f cm" % (dx, dy)).encode() + NL.encode() + cont \
        + (NL + "Q").encode()
    xref = pg.get_contents()[0]
    d.update_stream(xref, cont)
    pg.set_contents(xref)
    if neww or newh:
        pg.set_mediabox(pymupdf.Rect(0, 0, neww or pg.rect.width, newh or oldh))
        pg = d[0]
    if band:
        for x0, y0, x1, y1, c in band:
            pg.draw_rect(pymupdf.Rect(x0 + dx, y0, x1 + dx, y1),
                         color=None, fill=_rgb(c), width=0)
        print("   구획 색 띠 %d 칸" % len(band))
    for _cx, _y0, _y1 in leads:
        pg.draw_line(pymupdf.Point(_cx + dx, _y0), pymupdf.Point(_cx + dx, _y1),
                     color=(0.604, 0.639, 0.694), width=0.5)
    if cut:
        # redaction 은 선을 자르지 않는다 - IF_COVERED(기본)는 통째로 든 것만 지운다.
        # 긴 그래프(D5)는 그래서 지워졌지만, A 의 위,아래 틀 선은 걸치기만 해서 남았다.
        # ⇒ A 는 IF_TOUCHED 로 걸친 선까지 지우고 틀을 다시 그린다(box).
        pg.add_redact_annot(pymupdf.Rect(*cut))
        pg.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE,
                            graphics=(pymupdf.PDF_REDACT_LINE_ART_REMOVE_IF_TOUCHED
                                      if box else
                                      pymupdf.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED))
        print("   지운 자리 %s" % (cut,))
    if box:
        pg.draw_rect(pymupdf.Rect(*box), color=(0, 0, 0), width=0.8)
        print("   틀을 다시 그림 %s" % (box,))
    sv = io.open(os.path.join(DST, name + ".svg"), encoding="utf8").read()
    pg.insert_font(fontname="ari", fontfile=AR)
    pg.insert_font(fontname="arib", fontfile=ARB)
    put = 0
    for t in texts(sv):
        fn = "arib" if t["bold"] else "ari"
        w = FONT[fn].text_length(t["text"], t["size"])
        if t["rot"]:
            pg.insert_text(pymupdf.Point(t["x"] + dx, t["y"] + w / 2.0), t["text"],
                           fontname=fn, fontsize=t["size"], rotate=90)
        else:
            x = t["x"] - (w if t["anchor"] == "end" else
                          w / 2.0 if t["anchor"] == "middle" else 0.0)
            pg.insert_text(pymupdf.Point(x + dx, t["y"]), t["text"], fontname=fn,
                           fontsize=t["size"])
        put += 1
    d.save(os.path.join(DST, name + ".tmp.pdf"), deflate=True, garbage=3)
    d.close()
    os.replace(os.path.join(DST, name + ".tmp.pdf"), q)
    print("%-36s pdf : 옛 글자층 %d 덩이 걷고 %d 개 다시 적음" % (name, n0, put))


rebuild_pdf(D5, cut=CUT33)
# A : 창 밖(2240 bp 오른쪽)의 틀,눈금 그림을 지우고 새 오른쪽 선을 긋는다.
# 좌표는 통째 옮김(dx) 뒤의 자리로 준다 - redact,draw 는 옮김 다음에 일어난다.
rebuild_pdf(D1, dx=DX, neww=NEWW, newh=NEWH,
            cut=(XCUT + DX + 0.6, FY0 - 6.0, XR0 + DX + 3.0, FY1 + 6.0),
            box=(XL0 + DX, FY0, XCUT + DX, FY1), band=BAND, recolor=BARCOL,
            leads=[(cx, BY0 + BH, _y - T * 0.95) for _n, cx, _y, low in BTX if low])
