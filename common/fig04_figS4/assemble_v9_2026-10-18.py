# -*- coding: utf-8 -*-
# v9 (2026-10-18 밤) - 규격판 조립기
# ⒜ 줄 배율을 쓰지 않는다 - 배율은 글자 크기를 줄마다 다르게 만든다(규격 ⒝ 위반).
# 줄 길이는 칸 사이 틈에 남는 폭을 나눠 넣어 맞춘다(칸이 하나면 가운데).
# ⒝ 패널 글자 = 8.2 pt 대문자, 자리는 칸 구석이 아니라 틀 선에서
# 왼쪽 -22.5 pt, 기준선 틀 위 -1.5 pt.
# ⒞ 글자는 패널을 다 그린 뒤에 얹는다(먼저 얹으면 뒤 패널이 덮는다).
# ⒟ 폭은 환경변수 SPEC_W(기본 423.5 pt), 낼 곳은 OUT_SET.
# 패널은 이미 그 크기로 그려져 있어야 한다 - 여기서 늘이거나 줄이지 않는다.
import collections
import io
import os
import re
import shutil


def _need(k):
    """The deposited copy keeps no default: it stops rather than read the wrong place."""
    v = os.environ.get(k)
    if not v:
        raise SystemExit("environment variable %s is not set (the "
                         "deposited copy has no default; see README.md,"
                         " section 'Environment variables')" % k)
    return v


SP = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(SP, os.environ.get("PANEL_SRC", "panelsplit2"), "_panels")
BASE = os.path.join(_need("RDNA_ARRAY_METH_DIR"), "03_figures")
OLDD = os.path.join(BASE, "Fig4_submission_set_2026-10-09")
OLDPD = os.path.join(OLDD, "panels")
PREV = os.path.join(BASE, "00_final_set_2026-10-12")
OUTD = os.path.join(BASE, os.environ.get("OUT_SET", "00_최종세트_2026-10-16"))
PD = os.path.join(OUTD, "panels")
os.makedirs(PD, exist_ok=True)
NL = chr(10)
PAD, GAP, LABW = 10.0, 15.0, 13.0
TARGET_W = float(os.environ.get("SPEC_W", "423.5"))   # 발행 폭 (pt)
FS_LAB = 8.2                                          # 패널 글자
LAB_DX, LAB_DY = -22.5, -1.5                          # 틀 선에서
NUMRX = r"[-+]?[\d.]+"
INPAD = 6.0

FIG4 = [
    [("Fig4_a_a2_b_c2_c3_v3__a", "A"), ("Fig4_a_a2_b_c2_c3_v3__b", "B")],
    [("Fig4_a_a2_b_c2_c3_v3__a2", "C")],
    [("Fig4c_six_units_continuous__c", "D")],
    # G 삭제, E 와 F 가 그 줄의 반반
    [("Fig4c2c3c4_block_anatomy__c2", "E"), ("Claim_G_outside_array__b", "F")],
    # G 를 뺐으므로 뒤 글자를 한 칸씩 당긴다 (옛 H~P → 새 G~O)
    [("F39_5S_array__a", "G"), ("F39_5S_array__c", "H"), ("F39_5S_array__b", "I")],
    [("F39_D4Z4_4q_array__a", "J"), ("F39_D4Z4_4q_array__c", "K"),
     ("F39_D4Z4_4q_array__b", "L")],
    [("Claim_F2_show_hardcall_ext_EN__d10", "M"),
     ("Claim_F2_show_hardcall_ext_EN__d9", "N"),
     ("Claim_F2_show_hardcall_ext_EN__d11", "O")],
]

# S4 : 1 줄 = [위 A / 아래 C] + B, 그 아래 D,E,F, G,H, I,J,K
S41 = [
    [("COL", [("SuppS6_from_c", "A"), ("ClaimA_two_states__b", "C")]),
     ([["Fig4d_single_molecules__m1", "Fig4d_single_molecules__m2"],
       ["Fig4d_single_molecules__m3", "Fig4d_single_molecules__m4"]], "B")],
    [("ClaimA_two_states__c", "D"),
     ("Fig4c2c3c4_block_anatomy__c3", "E"), ("Fig4c2c3c4_block_anatomy__c4", "F")],
    [("F39_5S_array__d", "G"), ("F39_D4Z4_4q_array__d", "H")],
    [("F39_5S_array__e", "I"), ("F39_D4Z4_4q_array__e", "J"),
     ("Claim_F2_two_modes_of_transition__c", "K")],
]

S42 = [[("FIG_D1_CODING32__panel_a", "A")],
       [("FIG_D5_a4_persite_tone_hardcall32", "B")]]

ARRS = [("F39_D4Z4_10q_array", ["a", "c", "b", "e", "d"]),
        ("F39_DXZ4_array", ["a", "c", "b", "b-2", "e", "d"]),
        ("F39_RNU2_array", ["a", "c", "b", "e", "d"]),
        ("F39_SST1_array", ["a", "c", "b", "e", "d"])]
S44 = []
_lab = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
for _f, _p in ARRS:
    S44.append([("%s__%s" % (_f, q), next(_lab)) for q in _p if q not in ("d", "e")])
    S44.append([("%s__%s" % (_f, q), next(_lab)) for q in ("e", "d") if q in _p])


def fetch(b):
    got = []
    for ext in ("svg", "pdf"):
        s = os.path.join(SRC, b + "." + ext)
        for alt in (os.path.join(PREV, "panels", b + "." + ext),
                    os.path.join(OLDD, "panels_2026-10-11", b + "." + ext),
                    os.path.join(OLDPD, b + "." + ext)):
            if not os.path.exists(s):
                s = alt
        if os.path.exists(s):
            shutil.copy2(s, os.path.join(PD, b + "." + ext))
            got.append(ext)
    if "svg" not in got:
        raise SystemExit("패널 없음 : " + b)
    return "pdf" in got


def size_pt(p):
    s = io.open(p, encoding="utf8").read(400)
    return (float(re.search(r'width="([\d.]+)pt"', s).group(1)),
            float(re.search(r'height="([\d.]+)pt"', s).group(1)))


def frame_of(panel):
    """Top left of the axes rectangle (the axes_1 background patch) of the panel svg, in pt. None if it is not found."""
    s = io.open(os.path.join(PD, panel + ".svg"), encoding="utf8").read()
    i = s.find('<g id="axes_1">')
    if i < 0:
        return None
    m = re.search(r'<path d="M\s*(%s)\s+(%s)\s*(.*?)"' % (NUMRX, NUMRX),
                  s[i:i + 4000], re.S)
    if not m:
        return None
    xs, ys = [float(m.group(1))], [float(m.group(2))]
    for a, b in re.findall(r"L\s*(%s)\s+(%s)" % (NUMRX, NUMRX), m.group(3)):
        xs.append(float(a))
        ys.append(float(b))
    return min(xs), min(ys)


def inner(p, pref):
    s = io.open(p, encoding="utf8").read()
    i = s.find("<svg")
    j = s.find(">", i)
    body = s[j + 1:s.rfind("</svg>")]
    for k in sorted(set(re.findall(r'id="([^"]+)"', body)), key=len, reverse=True):
        nk = "%s_%s" % (pref, k)
        body = body.replace('id="%s"' % k, 'id="%s"' % nk)
        body = body.replace('#%s"' % k, '#%s"' % nk)
        body = body.replace('url(#%s)' % k, 'url(#%s)' % nk)
    return body


def as_grid(fil):
    if isinstance(fil, str):
        return [[fil]]
    return [r if isinstance(r, list) else [r] for r in fil]


def grid_size(g):
    w = h = 0.0
    for r in g:
        rw = sum(size_pt(os.path.join(PD, f + ".svg"))[0] for f in r) + INPAD * (len(r) - 1)
        rh = max(size_pt(os.path.join(PD, f + ".svg"))[1] for f in r)
        w = max(w, rw)
        h += rh + INPAD
    return w, h - INPAD


def cells_of(row):
    """줄 -> 칸 목록. 칸 = ('P', letter, grid) 또는 ('C', [(letter, grid), ...])"""
    out = []
    for a, b in row:
        if a == "COL":
            out.append(("C", [(lab, as_grid(f)) for f, lab in b]))
        else:
            out.append(("P", b, as_grid(a)))
    return out


def cell_size(c):
    if c[0] == "P":
        return grid_size(c[2])
    ws, hs = [], 0.0
    for _lab, g in c[1]:
        w, h = grid_size(g)
        ws.append(w)
        hs += h + GAP
    return max(ws), hs - GAP


def files_of(row):
    out = []
    for c in cells_of(row):
        gs = [c[2]] if c[0] == "P" else [g for _l, g in c[1]]
        for g in gs:
            for rr in g:
                out += rr
    return out


ROWS_OUT = []


def build(name, ROWS):
    need = [f for r in ROWS for f in files_of(r)]
    haspdf = {b: fetch(b) for b in need}
    left = [b for b in need
            if re.search(r">%s<" % re.escape(b.split("__")[-1]),
                         io.open(os.path.join(PD, b + ".svg"), encoding="utf8").read())]
    # ── 1 차 : 줄마다 자연 폭
    rows = [cells_of(r) for r in ROWS]
    natw = []
    for cs in rows:
        natw.append(sum(cell_size(c)[0] for c in cs) + PAD * (len(cs) - 1))
    target = max(natw + [TARGET_W - LABW])
    scales = [1.0 for _w in natw]          # 배율 금지 (규격 ⒝)
    slacks = [target - w for w in natw]    # 남는 폭은 틈으로 나눈다

    # ── 2 차 : 자리 잡기 (줄마다 배율)
    place, y, n = [], 0.0, 0
    for ri, (cs, s) in enumerate(zip(rows, scales)):
        gapx = PAD + (slacks[ri] / (len(cs) - 1) if len(cs) > 1 else 0.0)
        x = LABW + (slacks[ri] / 2.0 if len(cs) == 1 else 0.0)
        rh = 0.0
        for c in cs:
            cw, ch = cell_size(c)
            groups = [(c[1], c[2])] if c[0] == "P" else list(c[1])
            yy = y
            for lab, g in groups:
                gw, gh = grid_size(g)
                subs = []
                ty = yy
                for rr in g:
                    tx = x
                    rrh = 0.0
                    for f in rr:
                        fw, fh = size_pt(os.path.join(PD, f + ".svg"))
                        subs.append((f, tx, ty, fw * s, fh * s))
                        tx += fw * s + INPAD
                        rrh = max(rrh, fh * s)
                    ty += rrh + INPAD
                n += 1
                place.append((lab, subs, x, yy, gw * s, gh * s, n))
                yy += gh * s + GAP
            x += cw * s + gapx
            rh = max(rh, ch * s + (GAP * (len(groups) - 1)))
        y += rh + GAP
    H = y - GAP
    maxw = target + LABW

    out, idx = [], 0
    for lab, subs, x, yy, w, h, _i in place:
        for f, sx, sy, sw, sh in subs:
            idx += 1
            ow, oh = size_pt(os.path.join(PD, f + ".svg"))
            out.append('<g transform="translate(%.4f,%.4f)">' % (sx, sy))
            out.append('<svg width="%.4f" height="%.4f" viewBox="0 0 %.4f %.4f" '
                       'overflow="visible">' % (sw, sh, ow, oh))
            out.append(inner(os.path.join(PD, f + ".svg"), "r%d" % idx))
            out.append('</svg></g>')
    # 글자는 패널을 다 얹은 뒤에 - 그리고 자리는 틀 선에서 잰다
    LABXY = []
    for lab, subs, x, yy, w, h, _i in place:
        f0, sx0, sy0 = subs[0][0], subs[0][1], subs[0][2]
        fr = frame_of(f0)
        lx = max(0.5, (sx0 + fr[0] + LAB_DX) if fr else (x - LABW + 1.0))
        ly = (sy0 + fr[1] + LAB_DY) if fr else (yy + 10.0)
        LABXY.append((lab, lx, ly, fr is not None))
        out.append('<text x="%.4f" y="%.4f" font-family="Arial" font-size="%.1f" '
                   'font-weight="bold" fill="#000000">%s</text>'
                   % (lx, ly, FS_LAB, lab))
    txt = ('<?xml version="1.0" encoding="utf-8"?>' + NL
           + '<svg xmlns="http://www.w3.org/2000/svg" '
             'xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" '
             'width="%.2fpt" height="%.2fpt" viewBox="0 0 %.2f %.2f">' % (maxw, H, maxw, H)
           + NL + '<rect width="%.2f" height="%.2f" fill="#ffffff"/>' % (maxw, H) + NL
           + NL.join(out) + NL + "</svg>" + NL)
    io.open(os.path.join(OUTD, name + ".svg"), "w", encoding="utf8", newline=NL).write(txt)
    exp = sum(io.open(os.path.join(PD, f + ".svg"), encoding="utf8").read().count("<text")
              for _l, subs, _x, _y, _w, _h, _i in place for f, _a, _b, _c, _d in subs)
    got = txt.count("<text") - len(place)
    ids = re.findall(r'id="([^"]+)"', txt)
    dup = [k for k, v in collections.Counter(ids).items() if v > 1]
    refs = (set(re.findall(r'xlink:href="#([^"]+)"', txt))
            | set(re.findall(r'url\(#([^)]+)\)', txt)))
    nopdf = [b for b in need if not haspdf[b]]
    print("%-46s %7.1f x %7.1f pt · 칸 %2d (패널 %2d)" % (name, maxw, H, n, len(need)))
    print("     <text %d↔%d %s · id중복 %d · 끊긴참조 %d · 속letter %d%s"
          % (exp, got, "일치" if exp == got else "다름", len(dup),
             len(refs - set(ids)), len(left),
             " · pdf없음 %s" % nopdf if nopdf else ""))
    print("     폭 %.1f (목표 %.1f · 차 %+.1f) · 배율 없음(전부 1.00)"
          % (maxw, TARGET_W, maxw - TARGET_W))
    print("     줄 남는 폭(틈으로 나눔) : %s"
          % " · ".join("%.1f" % v for v in slacks))
    nofr = [l for l, _x, _y, ok in LABXY if not ok]
    print("     패널 글자 %.1f pt · 틀에서 잰 것 %d / %d%s"
          % (FS_LAB, len(LABXY) - len(nofr), len(LABXY),
             (" frame not found %s" % nofr) if nofr else ""))
    if not nopdf:
        import pymupdf
        doc = pymupdf.open()
        page = doc.new_page(width=maxw, height=H)
        page.draw_rect(pymupdf.Rect(0, 0, maxw, H), color=None, fill=(1, 1, 1))
        for lab, subs, x, yy, w, h, _i in place:
            for f, sx, sy, sw, sh in subs:
                src = pymupdf.open(os.path.join(PD, f + ".pdf"))
                page.show_pdf_page(pymupdf.Rect(sx, sy, sx + sw, sy + sh), src, 0)
                src.close()
        for lab, lx, ly, _ok in LABXY:          # 패널 뒤에 얹는다
            page.insert_text(pymupdf.Point(lx, ly), lab,
                             fontname="hebo", fontsize=FS_LAB, color=(0, 0, 0))
        doc.save(os.path.join(OUTD, name + ".pdf"))
        pg = pymupdf.open(os.path.join(OUTD, name + ".pdf"))[0]
        pg.get_pixmap(matrix=pymupdf.Matrix(150 / 72.0, 150 / 72.0)).save(
            os.path.join(OUTD, name + ".png"))
    for lab, subs, x, yy, w, h, _i in place:
        for f, sx, sy, sw, sh in subs:
            ROWS_OUT.append((name, lab, f, round(sx, 2), round(sy, 2),
                             round(sw, 2), round(sh, 2)))


ONLY = [x for x in os.environ.get("ONLY_FIG", "").split(",") if x]
for nm, R in (("Fig4_main_2026-10-16", FIG4),
              ("SuppFig_S4_blocks_and_molecules_2026-10-16", S41),
              ("SuppFig_S5_dispersed_5S_2026-10-16", S42),
              ("SuppFig_S6_other_macrosatellites_2026-10-16", S44)):
    if ONLY and not any(k in nm for k in ONLY):
        continue
    build(nm, R)

_tsv = os.path.join(OUTD, "00_panel_layout_2026-10-16.tsv")
_old = []
if ONLY and os.path.exists(_tsv):
    _old = [l for l in io.open(_tsv, encoding="utf8").read().split(NL)[1:]
            if l.strip() and not any(k in l.split(chr(9))[0] for k in ONLY)]
with io.open(_tsv, "w", encoding="utf8", newline=NL) as o:
    o.write(chr(9).join(["figure", "slot", "panel_file", "x_pt", "y_pt",
                         "w_pt", "h_pt"]) + NL)
    for l in _old:
        o.write(l + NL)
    for r in ROWS_OUT:
        o.write("\t".join(str(x) for x in r) + NL)
print()
print("배치표 %d 행" % len(ROWS_OUT))
