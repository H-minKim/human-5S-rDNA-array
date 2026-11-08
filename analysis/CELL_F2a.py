# -*- coding: utf-8 -*-
r"""Main Fig 2 cell a: only panel A of FIG_UNITPOS.py, at cell size (215 x 76.5 pt)
   Data and computation are the same lines as FIG_UNITPOS.py (v11_io.load, H1 557, 51,801 units)
   Only the appearance changed: figure size, font size, line width, axis title line breaks, no panel title (the composite removes it anyway)
   Content check: CMPART.py (old fig_unitpos_by_region axes 0 vs this figure)
Output: $FIGOUT\cellF2_a_unitpos.png/svg, CELL_F2a_source_values.tsv
"""
import io, os, sys, collections
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import matplotlib.text
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir, "common", "fig02_figS2")))
import v11_io as G

_DISP = {"AluY": "AluY element", "NTS1": "spacer-1", "NTS2": "spacer-2",
         "GArich": "GA-rich", "GA-rich": "GA-rich",
         "promoter": "5′-flanking", "prom": "5′-flanking", "프로모터": "5′-flanking",
         "gene": "5S gene", "coding": "5S gene"}
def disp(s):
    return _DISP.get(s, s)


OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
BLOCKS = G.BLOCKS
ELEM = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA-rich", 1179, 1551),
        ("NTS2", 1552, 2031), ("promoter", 2032, 2121), ("gene", 2122, 2240)]
UL = 2240

def pos_of(t):
    if t.startswith("DB") or t.startswith("IB"):
        b = int(t[2:].split(":")[0]); return BLOCKS[b - 1][0] if 1 <= b <= len(BLOCKS) else None
    if t[0] == "D": return int(t[1:].split(":")[0])
    if t[0] == "I": return int(t[1:].split(":")[0])
    if t[0] == "R": return int(t[1:-2])
    try: return int(t[:-2])
    except Exception: return None

U, _, _ = G.load("primary", verbose=False)
def zone5(ap, ad):
    if ap <= ad:
        return "근위1-10" if ap <= 10 else ("근위11-30" if ap <= 30 else "가운데")
    return "원위1-10" if ad <= 10 else ("원위11-30" if ad <= 30 else "가운데")
Z5 = ["근위1-10", "근위11-30", "가운데", "원위11-30", "원위1-10"]

SITE = {z: collections.Counter() for z in Z5}
SNV = {z: collections.Counter() for z in Z5}
BLK = {z: collections.Counter() for z in Z5}
NZ = collections.Counter()
TYPES = collections.defaultdict(set)
occ = 0
for k in sorted(U):
    for r in U[k]:
        z = zone5(r[0], r[1]); NZ[z] += 1
        seen = set(); seenb = set(); seen2 = set()
        for t in r[4]:
            p = pos_of(t)
            if p is None or not (1 <= p <= UL): continue
            TYPES[p].add(t); occ += 1
            if p not in seen:
                SITE[z][p] += 1; seen.add(p)
            if t.startswith("DB") or t.startswith("IB"):
                _b = int(t[2:].split(":")[0])
                if _b not in seenb:
                    BLK[z][_b] += 1; seenb.add(_b)
            elif p not in seen2:
                SNV[z][p] += 1; seen2.add(p)
ALLP = sorted(set().union(*[set(SITE[z]) for z in Z5]))
print("유닛 %d · 변이 출현 %d · 자리 %d · 종류 %d" % (sum(NZ.values()), occ, len(ALLP), sum(len(v) for v in TYPES.values())))

NALL = sum(NZ[z] for z in Z5)
_SNV = {}
for z in Z5:
    for q_, c_ in SNV[z].items(): _SNV[q_] = _SNV.get(q_, 0) + c_
_BLK = {}
for z in Z5:
    for q_, c_ in BLK[z].items(): _BLK[q_] = _BLK.get(q_, 0) + c_
psS = sorted(_SNV); vS = np.array([_SNV[q_] / NALL for q_ in psS], float)
BARS = sorted(((G.BLOCKS[i - 1][0], G.BLOCKS[i - 1][1], n / NALL, i)
               for i, n in _BLK.items() if 1 <= i <= len(G.BLOCKS)), key=lambda x: -x[2])
SRC = io.open(os.path.join(OUT, "CELL_F2a_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("종류\tbp_또는_블록\t시작\t끝\t비율\t출처\n")
for q_, v_ in zip(psS, vS):
    SRC.write("블록 밖\t%d\t%d\t%d\t%.6f\tv11_io.load('primary')\n" % (q_, q_, q_, v_))
for _a, _b2, _v, _i in BARS:
    SRC.write("블록 안\t%d\t%d\t%d\t%.6f\tv11_io.load('primary')\n" % (_i, _a, _b2, _v))
SRC.close()

# ---------- cell figure (215 x 76.5 pt) ----------
fig, ax = CS.cell_at("a")   # margins and size from the single cellstyle.ROWS table (plot frame alignment)
# only indicators needed to know what it is in the figure, the legend only as far as 'variants outside blocks, per base' and 'variants inside repeat blocks',
#   the bracketed explanations (number of units, drawn at block width, positions cannot be resolved) and the value marks for blocks 46, 56 and 58 go to the Legend text (values stay in the value table CELL_F2a)
# legend text as actually counted, blue line = all substitutions (including substitutions inside blocks) + deletions and insertions outside blocks (all tokens that are not DB or IB),
#   orange = only deletions and insertions inside blocks (DB, IB). The old 'variants outside blocks, per base' and 'variants inside repeat blocks' read as if substitutions inside blocks were orange too
ax.plot(psS, vS, "-", lw=0.55, color="#2471a3", alpha=0.95, label="반복 블록 밖 결실과 삽입, 그리고 모든 치환")
for _a, _b2, _v, _i in BARS:
    ax.add_patch(mp.Rectangle((_a, 0.0), max(_b2 - _a, 1), _v, fc="#E8A33D", ec="none",
                              alpha=0.80, zorder=3))
ax.plot([], [], "s", color="#E8A33D", alpha=0.80, ms=3.5, label="반복 블록 안 결실과 삽입")
# move the unit element strip outside the axes, below the x tick labels
#   (bottom margin of the first row in cellstyle.ROWS 17 -> 26 pt, the strip, transcription arrow, names and leader lines keep the second version's shape, positions set in pt and converted by the axes height ratio)
YTOP = max(0.95, float(max(vS.max(), BARS[0][2])) * 1.30)
AXH = ax.get_position().height * fig.get_figheight() * 72.0   # axes height in pt (cell a = 60.5)
YBOT = -1.2 * YTOP / (AXH - 1.2)                              # only 1.2 pt so the 0 line does not touch the bottom axis line
for _, a, b in ELEM:
    ax.axvline(b, color="#909497", lw=0.3, ls=":", zorder=0)
ax.set_xlim(0, UL); ax.set_ylim(YBOT, YTOP)
ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8]); ax.set_ylim(YBOT, YTOP)   # ticks keep the old values (no negative ticks)
ax.set_xticks(range(0, UL + 1, 250))
ax.set_ylabel(CS.L("변이율", "Variant frequency"))   # 'something like variant frequency' (old: fraction differing from the reference / Fraction != reference, meaning in the Legend)   # 2026-10-18 night: axis title on one line, short (old: fraction of units differing from the reference sequence at that position)
_lg = ax.legend(loc="upper left", frameon=False)   # remove the legend frame
# the legend text got longer, count whether blue line points or orange bars fall inside the legend box (must be 0)
fig.canvas.draw(); _lb = _lg.get_window_extent(fig.canvas.get_renderer())
_lxy = ax.transData.transform(np.column_stack([psS, vS]))
print("범례 상자 겹침 : 선 점 %d · 막대 %d · 상자 오른끝 %.1f / 축 오른끝 %.1f px" % (
      sum(1 for x_, y_ in _lxy if _lb.contains(x_, y_)),
      sum(1 for p_ in ax.patches if p_.get_window_extent().overlaps(_lb)), _lb.x1, ax.get_window_extent().x1))
# position below the axes: measure the bottom of the x tick labels and place the strip below it
fig.canvas.draw(); _rd = fig.canvas.get_renderer(); PXPT = fig.dpi / 72.0
_axb = ax.get_window_extent(_rd)
TICK_PT = (_axb.y0 - min(t.get_window_extent(_rd).y0 for t in ax.get_xticklabels() if t.get_text())) / PXPT
GAP_PT, BAR_PT, ARR_PT = 1.5, 6.2, 3.4                         # gap below tick labels, strip, transcription arrow row (name row below it)
B_TOP = TICK_PT + GAP_PT; B_BOT = B_TOP + BAR_PT; Y_ARR = B_BOT + 0.70 * ARR_PT; Y_LAB = B_BOT + ARR_PT   # pt from the bottom axis line
Y_BASE = Y_LAB + 0.95 * CS.FS_ANN + 0.8                        # common baseline of the two names below the strip (baseline of the old top alignment ~ Y_LAB + 0.95 FS, 0.8 pt below that)
fy = lambda pt: -pt / AXH                                      # pt below the bottom axis line (axes fraction)
XT = ax.get_xaxis_transform()                                  # x = bp, y = axes fraction
INK = "#17202a"                                                # leader line and arrow colour (grey 566573 was faint, width tied to 0.45)
BELOW = {"promoter": 1985, "gene": UL}                        # right end (bp) of the names below the strip, so they do not overlap on one row
LEADX = {"promoter": 2050, "gene": 2215}                      # where the 5S gene leader line meets the strip (bp, inside the box), avoiding the arrow (2122-2185)
TSS, TSS_END = 2122, 2185                                     # 5S gene start = transcription start, arrow end
# strip colours from the unit element colours of Fig 1 A (fig01_A_B_array_and_pipeline_schematic.py ELEM, box line = its EDGE)
ECOL = {"AluY": (0.725, 0.804, 0.894), "NTS1": (0.890, 0.925, 0.965), "GA-rich": (0.953, 0.875, 0.722),
        "NTS2": (0.890, 0.925, 0.965), "promoter": (0.835, 0.878, 0.812), "gene": (0.918, 0.722, 0.682)}
EEDGE = (0.29, 0.478, 0.69)
_RECT = {}
for n_, a, b in ELEM:
    _RECT[n_] = ax.add_patch(mp.Rectangle((a, fy(B_BOT)), b - a, BAR_PT / AXH, transform=XT, fc=ECOL[n_], ec=EEDGE,
                                          lw=0.4, clip_on=False, zorder=5))
    if n_ not in BELOW:   # the four wide ones: name inside the strip (enlarged to FS_ANN on save, so that size from the start)
        ax.text((a + b) / 2, fy(B_TOP + BAR_PT / 2), disp(n_), transform=XT, ha="center", va="center", fontsize=CS.FS_ANN, zorder=6, clip_on=False)
    else:                 # the two narrow ones: one row below the strip, straight dark leader line
        # "5'-flanking and 5S gene very slightly lower, same y for both": va top -> baseline (with top the ascender differs per text, so the two rows misalign), 0.8 pt lower
        ax.text(BELOW[n_], fy(Y_BASE), disp(n_), transform=XT, ha="right", va="baseline", fontsize=CS.FS_ANN, zorder=7, clip_on=False)
        _yg = Y_BASE - 0.73 * CS.FS_ANN                        # top of the text (above the baseline by the cap and digit height)
        if n_ == "promoter":   # straight from the top right corner of the name to the centre of the box
            ax.plot([BELOW[n_] + 4, (a + b) / 2], [fy(_yg), fy(B_BOT)], "-", transform=XT, lw=0.45, color=INK, clip_on=False, zorder=7, label="_nolegend_")
        else:                  # straight from above the name into the box (avoiding the arrow)
            ax.plot([LEADX[n_], LEADX[n_]], [fy(_yg), fy(B_BOT)], "-", transform=XT, lw=0.45, color=INK, clip_on=False, zorder=7, label="_nolegend_")
# the same bent arrow as Fig 1 A, down from the gene start and then right (transcription direction)
ax.plot([TSS, TSS], [fy(B_BOT), fy(Y_ARR)], "-", transform=XT, lw=0.45, color=INK, clip_on=False, zorder=8, solid_capstyle="butt", label="_nolegend_")
ax.annotate("", xy=(TSS_END, fy(Y_ARR)), xytext=(TSS, fy(Y_ARR)), xycoords=XT, textcoords=XT, annotation_clip=False, zorder=8,
            arrowprops=dict(arrowstyle="-|>", mutation_scale=3.2, lw=0.45, color=INK, shrinkA=0, shrinkB=0))
# Checks: names inside the strip stay in their box, names below the strip do not overlap, lines and arrows do not touch other names, the schematic is below the tick labels, nothing goes outside the figure
fig.canvas.draw(); _rd = fig.canvas.get_renderer(); _bad = []
_lab = [t for t in ax.texts if t.get_text() in [disp(n) for n, _, _ in ELEM]]
_low = {}
for t in _lab:
    bb = matplotlib.text.Text.get_window_extent(t, _rd)   # an annotation includes its leader line in the box, so text only
    n_ = next(e for e in ELEM if disp(e[0]) == t.get_text())[0]
    if n_ not in BELOW:
        rb = _RECT[n_].get_window_extent(_rd)
        if not (bb.y0 >= rb.y0 - 0.01 and bb.y1 <= rb.y1 + 0.01 and bb.x0 >= rb.x0 and bb.x1 <= rb.x1): _bad.append(("띠 밖", t.get_text()))
    else:
        _low[t.get_text()] = bb
        if bb.y0 < 0.5 or bb.x1 > fig.bbox.x1: _bad.append(("그림 밖", t.get_text(), round(bb.y0, 2)))
_bl = list(_low.values())
if len(_bl) == 2 and _bl[0].overlaps(_bl[1]): _bad.append(("아래 이름 겹침",))
_ink = [l for l in ax.lines if l.get_color() == INK] + [t.arrow_patch for t in ax.texts if getattr(t, 'arrow_patch', None) is not None and t.get_text() == '']
for k, o in enumerate(_ink):
    for nm, tb in _low.items():
        if o in ax.lines:
            if (nm == disp("promoter") and abs(o.get_xdata()[0] - BELOW["promoter"] - 4) < 1e-9) or (nm == disp("gene") and abs(o.get_xdata()[0] - LEADX["gene"]) < 1e-9): continue
            _xy = o.get_transform().transform(list(zip(o.get_xdata(), o.get_ydata())))   # measuring a diagonal by its box gives false alarms, 200 points on the segment
            _hit = any(tb.contains(*(p0 + (p1 - p0) * u)) for p0, p1 in zip(_xy[:-1], _xy[1:]) for u in [i / 199.0 for i in range(200)])
        else:
            _hit = o.get_window_extent(_rd).overlaps(tb)
        if _hit: _bad.append(("선이 이름에 닿음", k, nm))
_tk = min(t.get_window_extent(_rd).y0 for t in ax.get_xticklabels() if t.get_text())
if max(r.get_window_extent(_rd).y1 for r in _RECT.values()) > _tk + 0.01: _bad.append(("띠가 눈금 글자에 닿음",))
_low_all = min([r.get_window_extent(_rd).y0 for r in _RECT.values()] + [b.y0 for b in _low.values()])
print("띠 글자 검사 : 축 아래 · 눈금 글자 %.1f pt · 띠 %.1f pt · 이름 %d · 그림 아래 여유 %.1f pt · 문제 %d %s" % (TICK_PT, BAR_PT, len(_lab), _low_all / PXPT, len(_bad), _bad))
CS.save(fig, OUT, "cellF2_a_unitpos")
print("저장", OUT)
