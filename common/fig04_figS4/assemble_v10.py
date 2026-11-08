# -*- coding: utf-8 -*-
# 조립기 v10 - 틀(축 사각형) 기준
#
# v9 까지는 칸(그림 전체 상자)으로 붙였다 ⇒ y 축 이름 길이가 다르면 그래프 선이 어긋났다.
# v10 은 틀로 붙인다 :
# , 같은 세로줄 : 틀 왼쪽 x, 오른쪽 x 가 같다 (맨 왼쪽 41.90, 맨 오른쪽 408.17)
# , 같은 가로줄 : 틀 위 y 가 같고 틀 높이가 같다
# , 패널 글자 : 그 틀에서 왼쪽 -22.48, 기준선 위 -1.53 (8.2 pt 대문자)
# 틀 크기는 paper_spec_hook 이 이미 못박아 두었다(SPEC_FRAME) - 여기서는 안 늘이고 안 줄인다.
import collections
import io
import json
import os
import re
import shutil
import sys

SP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SP)
from frame_layout import X_L, X_R, GAP_X, GAP_Y, INPAD, margins   # noqa


def _need(k):
    """The deposited copy keeps no default: it stops rather than read the wrong place."""
    v = os.environ.get(k)
    if not v:
        raise SystemExit("environment variable %s is not set (the "
                         "deposited copy has no default; see README.md,"
                         " section 'Environment variables')" % k)
    return v


# ㉧ : 출력 뿌리를 OUT_BASE 로 받는다 - 저장소만으로 돌릴 때
# So that nothing is written outside the output folder. The default is unchanged.
BASE = (os.environ.get("OUT_BASE")
        or os.path.join(_need("RDNA_ARRAY_METH_DIR"), "03_figures"))
PREV = os.path.join(BASE, "00_final_set_2026-10-12")
OLDD = os.path.join(BASE, "Fig4_submission_set_2026-10-09")
OUTD = os.path.join(BASE, os.environ.get("OUT_SET", "00_시험_423pt_2026-10-18"))
PD = os.path.join(OUTD, "panels")
SRC = os.path.join(SP, os.environ.get("PANEL_SRC", "panelsplit_spec"), "_panels")
os.makedirs(PD, exist_ok=True)
NL = chr(10)
W_TOTAL = float(os.environ.get("SPEC_W", "423.5"))
LAB_DX, LAB_DY, FS_LAB = -22.48, -1.53, 8.2
# 칸이 쪽(423.5)을 넘으면 오른쪽 틀 선을 그만큼 당긴다
XR_ADJ = [0.0]
# 남의 패널이라 다시 못 그리는 그림은 통째로 줄인다.
# 칸(그래프)만 줄이고 패널 글자는 8.2 pt 로 그대로 둔다.
# 글자 하한 4.5 는 깨진다 - 그 사실을 찍고 기록에 남긴다.
FIG_SCALE = {"S5": None}        # None 이면 쪽 폭에 맞춰 자동
# 칸마다 더 왼쪽으로 (세로축 이름이 길어 글자와 붙는 칸) - (그림, 글자) : 더할 값
LAB_DX_EXTRA = {("Fig4_main", "N"): -9.0,
                # 칸 글자가 왼쪽 이웃의 세로축 이름에 닿았다
                # (실측 E↔Low-methylation blocks, F↔High-…, J↔β)
                ("S4", "E"): -5.0, ("S4", "F"): -5.0, ("S4", "J"): -2.5,
                # 칸 글자가 제 틀 위의 범례,세로축 이름에 닿은 자리(크롭으로 눈 확인)
                ("Fig4_main", "F"): -3.0,
                ("S6", "C"): -3.0, ("S6", "G"): -6.5, ("S6", "H"): -3.0,
                ("S6", "R"): -3.0, ("S6", "S"): -3.0}
# 칸을 통째로 오른쪽으로 미는 값 - {(그림, 글자): pt}
# 줄마다 위로 당기는 값 - {(그림, 줄번호 0부터): pt}
ROW_DY = {("Fig4_main", 2): -10.0, ("Fig4_main", 3): -10.0,
          ("S4", 1): -10.0, ("S4", 2): -10.0,
          ("S6", 2): -10.0, ("S6", 5): -10.0}
CELL_DX = {("S6", "A"): 10.0, ("S6", "F"): 10.0, ("S6", "L"): 10.0,
           ("S6", "Q"): 10.0, ("Fig4_main", "G"): 10.0,
           ("Fig4_main", "H"): 6.0}

# ── 배치표는 v9 에서 그대로 읽는다 (한 벌만 둔다)
_src = io.open(os.path.join(SP, "assemble_v9_2026-10-18.py"), encoding="utf8").read()
_ns = {"__file__": os.path.join(SP, "assemble_v9_2026-10-18.py")}
exec(compile(_src[:_src.index("def fetch(")], "head", "exec"), _ns)
FIG4, S41, S42, S44 = _ns["FIG4"], _ns["S41"], _ns["S42"], _ns["S44"]


def fetch(b):
    got = []
    for ext in ("svg", "pdf"):
        s = os.path.join(SRC, b + "." + ext)
        for alt in (os.path.join(PREV, "panels", b + "." + ext),
                    os.path.join(OLDD, "panels_2026-10-11", b + "." + ext),
                    os.path.join(OLDD, "panels", b + "." + ext)):
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


def inner(p, pref):
    s = io.open(p, encoding="utf8").read()
    i = s.find("<svg")
    j = s.find(">", i)
    body = s[j + 1:s.rfind("</svg>")]
    for k in sorted(set(re.findall(r'id="([^"]+)"', body)), key=len, reverse=True):
        nk = "%s_%s" % (pref, k)
        body = body.replace('id="%s"' % k, 'id="%s"' % nk)
        body = body.replace('#%s"' % k, '#%s"' % nk)
        body = body.replace("url(#%s)" % k, "url(#%s)" % nk)
    return body


def columns_of(row):
    """줄 -> [(letter, [세로로 쌓인 (파일, 가로칸번호) ...])] 을 가로칸 차례로"""
    cols = []                      # [(letter, [[파일,...] (세로줄)], 가로칸 폭 몫)]
    for f, lab in row:
        if f == "COL":
            sub = [(s_f, s_l) for s_f, s_l in lab]
            cols.append(("COL", sub))
        elif isinstance(f, list):
            grid = [rr if isinstance(rr, list) else [rr] for rr in f]
            cols.append(("GRID", (grid, lab)))
        else:
            cols.append(("P", (f, lab)))
    return cols


def _cell_dx(fig, lab):
    for (f, l), v in CELL_DX.items():
        if f in fig and l == lab:
            return v
    return 0.0


def _extra(fig, lab):
    for (f, l), v in LAB_DX_EXTRA.items():
        if f in fig and l == lab:
            return v
    return 0.0


RETRY = [False]
SCALE = [1.0]


def build(name, rows):
    need = []
    for r in rows:
        for kind, v in columns_of(r):
            if kind == "COL":
                need += [f for f, _l in v]
            elif kind == "GRID":
                need += [x for rr in v[0] for x in rr]
            else:
                need.append(v[0])
    haspdf = {b: fetch(b) for b in need}
    MG = {b: margins(os.path.join(PD, b + ".svg")) for b in need}

    place, labs, y = [], [], 0.0
    for _ri, row in enumerate(rows):
        n0 = len(place)
        for (_f, _r), _v in ROW_DY.items():
            if _f in name and _r == _ri:
                y += _v          # 그 줄부터 위로 당긴다
        cols = columns_of(row)
        # ── 이 줄의 가로 칸 목록 (격자는 그 안에서 또 나뉜다)
        flat = []            # [(파일, 칸idx)]
        for ci, (kind, v) in enumerate(cols):
            if kind == "COL":
                flat += [(f, ci) for f, _l in v]
            elif kind == "GRID":
                flat += [(x, ci) for rr in v[0] for x in rr]
            else:
                flat.append((v[0], ci))
        FH = max(MG[f][5] for f, _ci in flat)
        topmax = max(MG[f][2] for f, _ci in flat)
        botmax = max(MG[f][3] for f, _ci in flat)
        fy = y + topmax                      # 이 줄의 틀 위 y (모두 같다)
        # 틈은 고정이 아니라 그 줄에서 계산한다 - 맨 오른쪽 틀이 X_R 에 닿아야 한다.
        # (틀을 못박은 뒤 여백이 조금 달라지므로 계획 때의 틈으로는 안 맞는다)
        seq = []                              # [(왼여백, 틀폭, 오른여백)]
        for ci, (kind, v) in enumerate(cols):
            if kind == "P":
                f = v[0]
                seq.append((MG[f][0], MG[f][4], MG[f][1]))
            elif kind == "COL":
                fs = [f for f, _l in v]
                seq.append((max(MG[f][0] for f in fs), max(MG[f][4] for f in fs),
                            max(MG[f][1] for f in fs)))
            else:
                grid = v[0]
                ncol = max(len(rr) for rr in grid)
                for j in range(ncol):
                    fs = [rr[j] for rr in grid if j < len(rr)]
                    seq.append((max(MG[f][0] for f in fs), max(MG[f][4] for f in fs),
                                max(MG[f][1] for f in fs)))
        _XR = X_R - XR_ADJ[0]
        need_w = sum(w for _l, w, _r in seq) + sum(seq[j - 1][2] + seq[j][0]
                                                 for j in range(1, len(seq)))
        GAP = ((_XR - X_L) - need_w) / (len(seq) - 1) if len(seq) > 1 else 0.0
        # 세로칸마다 틀 왼쪽 x 를 미리 잡는다 - 다음 칸의 왼쪽 여백까지 더해야
        # 맨 오른쪽 틀이 X_R 에 닿는다(앞 판에서 이 항을 빠뜨려 36 pt 모자랐다).
        XS, cur = [], X_L
        for j, (L, w, R) in enumerate(seq):
            XS.append(cur)
            if j + 1 < len(seq):
                cur = cur + w + R + GAP + seq[j + 1][0]
        # 한 줄 안에서 세로로 쌓인 칸들(COL, 격자)의
        # 같은 번째 줄은 틀 위 y 가 같아야 한다 - 그래야 아래 선도 맞는다
        # (S4 의 C 와 B 아래쪽이 어긋나 있었다).
        nsub = 1
        for kind, v in cols:
            if kind == "COL":
                nsub = max(nsub, len(v))
            elif kind == "GRID":
                nsub = max(nsub, len(v[0]))
        SUBY, yy0 = [], fy
        for j in range(nsub):
            SUBY.append(yy0)
            hmax = 0.0
            for kind, v in cols:
                fs = []
                if kind == "COL" and j < len(v):
                    fs = [v[j][0]]
                elif kind == "GRID" and j < len(v[0]):
                    fs = list(v[0][j])
                elif kind == "P" and j == 0:
                    fs = [v[0]]
                for f in fs:
                    hmax = max(hmax, MG[f][5] + MG[f][3] + MG[f][2])
            yy0 += hmax + GAP_Y * 0.5
        jj = 0
        for ci, (kind, v) in enumerate(cols):
            if kind == "P":
                f, lab = v
                L, R, T, B, fw, fh = MG[f]
                x = XS[jj] + _cell_dx(name, lab); jj += 1
                place.append((f, x - L, fy - T))
                labs.append((lab, x + LAB_DX + _extra(name, lab), fy + LAB_DY))
            elif kind == "COL":
                x = XS[jj]; jj += 1
                for j, (f, lab) in enumerate(v):
                    L, R, T, B, fw, fh = MG[f]
                    place.append((f, x - L, SUBY[j] - T))
                    labs.append((lab, x + LAB_DX + _extra(name, lab), SUBY[j] + LAB_DY))
            else:
                grid, lab = v
                ncol = max(len(rr) for rr in grid)
                xs_g = XS[jj:jj + ncol]; jj += ncol
                first = True
                for j, rr in enumerate(grid):
                    for j2, f in enumerate(rr):
                        L, R, T, B, fw, fh = MG[f]
                        place.append((f, xs_g[j2] - L, SUBY[j] - T))
                        if first:
                            labs.append((lab, xs_g[0] + LAB_DX + _extra(name, lab),
                                         SUBY[j] + LAB_DY))
                            first = False
        # 다음 줄 y 를 어림식이 아니라 이 줄에서 실제로 놓인 상자들의
        # 맨 아래로 잡는다 - 앞 판은 어림이 모자라 다음 줄이 올라타 겹쳤다.
        bot = 0.0
        for f, px, py in place[n0:]:
            _w, _h = size_pt(os.path.join(PD, f + ".svg"))
            bot = max(bot, py + _h)
        y = bot + GAP_Y
    H = y - GAP_Y

    # 쪽 폭을 무조건 423.5 로 쓰면 그보다 넓은 칸이 잘린다.
    # S5(남의 패널, 1,116 pt)가 그렇게 잘려 나갔다.
    # ⇒ 내용이 더 넓으면 그 폭을 쓴다 - 잘린 그림을 내보내느니 규격 밖이라고 적는다.
    WNEED = 0.0
    for f, px, py in place:
        WNEED = max(WNEED, px + size_pt(os.path.join(PD, f + ".svg"))[0])
    # (둘째 흠) : 틀맞춤은 칸을 왼쪽으로도 민다(px < 0).
    # 통째 축소가 오른쪽 끝(WNEED) 만 보면 왼쪽으로 나간 몫은 그대로 잘린다 -
    # S5 의 D1 칸이 -58.9 pt 에 놓여 spacer-1 + … 라벨 앞머리가 사라졌다.
    # ⇒ 왼쪽,위로 나간 몫을 폭,배율에 함께 넣고 그만큼 밀어 넣는다.
    XMIN = min([px for f, px, py in place] + [lx for l, lx, ly in labs] + [0.0])
    YMIN = min([py for f, px, py in place] + [0.0])
    # 통째 축소가 걸린 그림은 여기서 줄인다
    _sc = None
    for k, v in FIG_SCALE.items():
        if k in name:
            _sc = v if v else (W_TOTAL - 1.0) / (WNEED - XMIN)
    if _sc and _sc < 0.999:
        place = [(f, (px - XMIN) * _sc, (py - YMIN) * _sc) for f, px, py in place]
        labs = [(l, (lx - XMIN) * _sc, (ly - YMIN) * _sc) for l, lx, ly in labs]
        H = (H - YMIN) * _sc
        WNEED = (WNEED - XMIN) * _sc
        SCALE[0] = _sc
        if XMIN < -0.05 or YMIN < -0.05:
            print("     왼쪽 %.1f · 위 %.1f pt 나가 있던 것을 안으로 밀어 넣었다"
                  % (-XMIN * _sc, -YMIN * _sc))
        # 칸 글자(8.2 pt)는 축소를 안 받으므로 줄인 그림에서는 쪽 위로 삐져나간다
        # (실측 S5 의 A 가 y0 = -4.9) ⇒ 그 먹물 높이만큼 아래로 민다.
        _top = min([ly - FS_LAB * 1.05 for l, lx, ly in labs]
                   + [py for f, px, py in place] + [0.0])
        if _top < -0.05:
            _d = -_top + 1.0        # +1 = 패널 자기 상자를 넘는 제목 몫
            place = [(f, px, py + _d) for f, px, py in place]
            labs = [(l, lx, ly + _d) for l, lx, ly in labs]
            H += _d
            print("     칸 글자가 쪽 위로 %.1f pt 나가 있어 %.1f pt 내렸다" % (-_top, _d))
        print("     통째로 %.4f 배로 줄였다 (다시 못 그리는 패널) — 글자 하한 4.5 가 깨진다"
              % _sc)
    PAGEW = W_TOTAL
    if WNEED > W_TOTAL + 0.5 and not RETRY[0]:
        return WNEED            # 부르는 쪽이 오른쪽 선을 당겨 다시 부른다
    if WNEED > W_TOTAL + 0.5:
        PAGEW = round(WNEED + 1.0, 2)
        print("     칸이 쪽보다 넓다 %.1f > %.1f ⇒ 쪽을 %.1f 로 (규격 밖 · 잘리지 않게)"
              % (WNEED, W_TOTAL, PAGEW))
    out, idx = [], 0
    for f, px, py in place:
        idx += 1
        ow, oh = size_pt(os.path.join(PD, f + ".svg"))
        out.append('<g transform="translate(%.4f,%.4f)">' % (px, py))
        out.append('<svg width="%.4f" height="%.4f" viewBox="0 0 %.4f %.4f" '
                   'overflow="visible">' % (ow * SCALE[0], oh * SCALE[0], ow, oh))
        out.append(inner(os.path.join(PD, f + ".svg"), "r%d" % idx))
        out.append("</svg></g>")
    for lab, lx, ly in labs:
        out.append('<text x="%.4f" y="%.4f" font-family="Arial" font-size="%.1f" '
                   'font-weight="bold" fill="#000000">%s</text>'
                   % (max(0.5, lx), ly, FS_LAB, lab))
    txt = ('<?xml version="1.0" encoding="utf-8"?>' + NL
           + '<svg xmlns="http://www.w3.org/2000/svg" '
             'xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" '
             'width="%.2fpt" height="%.2fpt" viewBox="0 0 %.2f %.2f">'
             % (PAGEW, H, PAGEW, H)
           + NL + '<rect width="%.2f" height="%.2f" fill="#ffffff"/>' % (PAGEW, H)
           + NL + NL.join(out) + NL + "</svg>" + NL)
    io.open(os.path.join(OUTD, name + ".svg"), "w", encoding="utf8",
            newline=NL).write(txt)
    # : 패널 svg 의 태그가 하나라도 안 맞으면
    # 조립본이 XML 로 깨지고 렌더러는 첫 </svg> 에서 멈춘다(칸 글자가 사라진다).
    # 겉보기 png 는 멀쩡해 보이므로 반드시 파싱해서 찍는다.
    import xml.etree.ElementTree as _ET
    try:
        _ET.fromstring(txt.encode("utf8"))
    except Exception as _e:
        print("     조립본 svg 가 XML 로 깨졌다 : %s ⇒ 패널 svg 의 <g> 짝을 보라" % _e)
    print("%-46s %7.2f x %7.2f pt · 칸 %2d" % (name, PAGEW, H, len(labs)))
    # ── 보고용 실측
    if not all(haspdf[b] for b in need):
        print("     pdf 없는 패널 — png·pdf 는 건너뜀")
        return
    import pymupdf
    doc = pymupdf.open()
    page = doc.new_page(width=PAGEW, height=H)
    page.draw_rect(pymupdf.Rect(0, 0, PAGEW, H), color=None, fill=(1, 1, 1))
    for f, px, py in place:
        ow, oh = size_pt(os.path.join(PD, f + ".svg"))
        src = pymupdf.open(os.path.join(PD, f + ".pdf"))
        page.show_pdf_page(pymupdf.Rect(px, py, px + ow * SCALE[0],
                                        py + oh * SCALE[0]), src, 0)
        src.close()
    for lab, lx, ly in labs:
        page.insert_text(pymupdf.Point(max(0.5, lx), ly), lab,
                         fontname="hebo", fontsize=FS_LAB, color=(0, 0, 0))
    doc.save(os.path.join(OUTD, name + ".pdf"))
    pg = pymupdf.open(os.path.join(OUTD, name + ".pdf"))[0]
    # png 는 300 dpi - NAR 래스터 하한
    pg.get_pixmap(matrix=pymupdf.Matrix(300 / 72.0, 300 / 72.0)).save(
        os.path.join(OUTD, name + ".png"))


ONLY = [x for x in os.environ.get("ONLY_FIG", "").split(",") if x]
for nm, R in (("Fig4_main_2026-10-16", FIG4),
              ("SuppFig_S4_blocks_and_molecules_2026-10-16", S41),
              ("SuppFig_S5_dispersed_5S_2026-10-16", S42),
              ("SuppFig_S6_other_macrosatellites_2026-10-16", S44)):
    if ONLY and not any(k in nm for k in ONLY):
        continue
    XR_ADJ[0] = 0.0
    RETRY[0] = False
    SCALE[0] = 1.0
    for _try in range(4):
        _over = build(nm, R)
        if not _over:
            break
        _d = _over - W_TOTAL + 0.5
        XR_ADJ[0] += _d
        print("     오른쪽 틀 선을 %.1f pt 당긴다 (%d 차)" % (XR_ADJ[0], _try + 1))
        if _try == 2:
            RETRY[0] = True      # 마지막엔 어쩔 수 없으면 쪽을 넓혀 잘리지 않게 낸다
