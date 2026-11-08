# -*- coding: utf-8 -*-
r"""Main Fig 2 cell p, only panel B of FIG_T2.py (including the overlaid right axis), at cell size (93.1 x 82 pt)
   Data and computation are the same lines as FIG_T2.py panel B (W2_marker_by_depth.tsv, S1_per_position_profile.tsv)
   Only the appearance changes: figure size, font size (shared by main cells), line and point sizes, axis label line breaks, text box line breaks, no panel title
   Content check: CMPART.py (old fig_T2_arrangement axes 1 and 2 <-> this figure axes 0 and 1)
Output : $FIGOUT\cellF2_p_depth.png/svg, CELL_F2p_source_values.tsv
"""
import io, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS
TAB = os.environ.get("FIGTAB", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/tables"))
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])


def read(nm):
    H, R = None, []
    for L in io.open(os.path.join(TAB, nm), encoding="utf8"):
        f = L.rstrip("\n").split("\t")
        if not f or f == [""]:
            continue
        if H is None:
            H = {c: i for i, c in enumerate(f)}; continue
        R.append(f)
    return H, R


SRC = io.open(os.path.join(OUT, "CELL_F2p_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("panel\t항목\t값\t출처파일\n")
def rec(p, k, v, f):
    SRC.write("%s\t%s\t%s\t%s\n" % (p, k, v, f)); return v


CT, CO, CG = "#c0392b", "#2471a3", "#909497"
# ================= B: by depth (same lines as FIG_T2.py) =================
H, R = read("W2_marker_by_depth.tsv")
d = [int(r[H["depth"]]) for r in R]
fD = [float(r[H["frac_term"]]) for r in R]   # fraction of terminal-type units (10 or more of the 12 diagnostic variants, termdef), the old frac_D1096_1 is the variant D1096:1 itself
nU = [int(r[H["n_units"]]) for r in R]
for r in R:
    rec("B", "깊이 %s" % r[H["depth"]], "말단형 %s · D1096:1 %s · 635CT %s · 1096TG %s · 1164GA %s (하플로 %s)"
        % (r[H["frac_term"]], r[H["frac_D1096_1"]], r[H["frac_635CT"]], r[H["frac_1096TG"]], r[H["frac_1164GA"]], r[H["n_units"]]),
        "W2_marker_by_depth.tsv")
fig, ax = CS.cell_at("p")   # margins and sizes from the single cellstyle.ROWS table (graph frame alignment)
ax.bar(d, fD, width=0.66, color=[CT if x > 0.10 else "#e59866" for x in fD], zorder=3)
# value text removed (values go to the value table) (old: '%.3f' for depths above 0.10 and for depth 1)
# 'depth 1 = %s (%d of %s)' belongs in the Results text, removed from the figure (values go to the value table)
rec("B", "깊이 1 (그림에서 뺀 글 · Result 로)", "깊이 1 = %s (%s 중 %d)" % (R[0][H["frac_term"]], R[0][H["n_units"]], round(fD[0] * nU[0])),
    "W2_marker_by_depth.tsv")
ax.set_xticks(d)
ax.tick_params(axis="x", labelrotation=90, pad=1.0)   # 15 ticks, written vertically because the cell is narrow (all numbers unchanged)
ax.set_xlabel(CS.L("원위 깊이 (유닛)", "Distal depth (units)"))   # axis label on one line, shorter (old: depth from the distal junction (units))
ax.set_ylabel(CS.L("말단형 보유 몫", "Frac. with terminal unit"))   # axis label on one line, shorter (old: fraction of haplotypes with a terminal-type unit at that depth)
ax.set_yticks([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])     # original ticks as is
ax.set_ylim(0, 0.35)
_SP = os.path.join(TAB, "S1_per_position_profile.tsv")
_H2, _R2 = None, []
for _L in io.open(_SP, encoding="utf8"):
    _f = _L.rstrip("\n").split("\t")
    if not _f or _f == [""]: continue
    if _H2 is None: _H2 = {c: i for i, c in enumerate(_f)}; continue
    _R2.append(_f)
_D2 = {}
for _r in _R2:
    if _r[_H2["stratum"]] != "all" or _r[_H2["metric"]] != "n_var" or _r[_H2["side"]] != "dist": continue
    _D2[int(_r[_H2["unit_from_junction"]])] = (float(_r[_H2["mean"]]), float(_r[_H2["ci_lo"]]),
                                               float(_r[_H2["ci_hi"]]))
_dk = [x for x in d if x in _D2]
_y = [_D2[x][0] for x in _dk]
ax2 = ax.twinx()
ax2.fill_between(_dk, [_D2[x][1] for x in _dk], [_D2[x][2] for x in _dk], color="#2471a3", alpha=0.16, lw=0, zorder=4)
ax2.plot(_dk, _y, "-o", color="#2471a3", ms=1.4, lw=0.6, zorder=5)
_pk = [x for x in _dk[1:-1] if _D2[x][0] > _D2[x - 1][0] and _D2[x][0] > _D2[x + 1][0]]
_pk = sorted(sorted(_pk, key=lambda x: -_D2[x][0])[:3])
# remove the blue inverted triangles (peak markers) (peak depths go to the value table)
for _x in _pk:
    rec("B", "변이 봉우리 깊이 %d" % _x, "유닛당 %.4f [%.4f, %.4f]" % _D2[_x], "S1_per_position_profile.tsv")
ax2.set_ylabel("유닛당 변이 수 (선)", color="#2471a3", labelpad=0.8)   # 1.5 -> 0.8 (Helvetica, the right axis label of the Korean edition is 1 pt inside the cell edge)
ax2.tick_params(axis="y", labelcolor="#2471a3")
ax2.set_yticks([6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0])   # original ticks as is (range from the original formula below)
ax2.set_ylim(min(_y) - 0.35, max(_y) + 1.15)
ax2.set_zorder(ax.get_zorder() + 1); ax.patch.set_visible(True); ax2.patch.set_visible(False)
rec("B", "변위 봉우리 깊이 (전부)", ",".join(str(x) for x in _pk), "S1_per_position_profile.tsv")
# The text box 'haplotypes %d (%d with depth 12 or more)' is not drawn on the figure; it is in the legend
rec("B", "분모 (그림에서 뺀 글 상자 · Legend 로)", "하플로타입 %d (깊이 12 이상은 %d)" % (nU[0], min(nU)), "W2_marker_by_depth.tsv")
SRC.close()
CS.save(fig, OUT, "cellF2_p_depth")
print("깊이 1 %s · 봉우리 %s" % (fD[0], _pk))
