# -*- coding: utf-8 -*-
r"""Main Fig 2 panel k - only version B of FIG_G3.py, at cell size (80.2 × 82 pt)
   Data and computation are the same lines as version B of FIG_G3.py (S14_length.tsv)
   Only the appearance changed: figure size, font size (common to main text cells), line width, line breaks of axis labels and notes, no panel title (the assembly removes it anyway)
   Content check: CMPART.py (old fig_G3_length_invariance axis 1 <-> this figure)
Output: $FIGOUT\cellF2_k_length.png/svg, CELL_F2k_source_values.tsv
"""
import io, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS

TAB = os.environ.get("FIGTAB", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/tables"))
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])


def read(nm):
    H, rows = None, []
    for L in __import__("taben").open_ko(os.path.join(TAB, nm), encoding="utf8"):
        f = L.rstrip("\n").split("\t")
        if not f or f == [""]:
            continue
        if H is None:
            H = {c: i for i, c in enumerate(f)}; continue
        rows.append(f)
    return H, rows


SRC = io.open(os.path.join(OUT, "CELL_F2k_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("panel\t항목\t값\t출처파일\n")


def rec(p, k, v, f):
    SRC.write("%s\t%s\t%s\t%s\n" % (p, k, v, f)); return v


# ================= B : by compartment (same lines as FIG_G3.py) =================
H, R = read("S14_length.tsv")
fig, ax = CS.cell_at("k")   # margins and sizes come from the single cellstyle.ROWS table (graph frame alignment)
BAS = sorted({r[H["basis"]] for r in R})
GRP = ["근위구획", "가운데구획", "원위구획"]
W = 0.36
for bi, b in enumerate(BAS):
    xs, ys, el, eh = [], [], [], []
    for gi, gname in enumerate(GRP):
        row = next((r for r in R if r[H["basis"]] == b and r[H["group"]] == gname), None)
        if row is None:
            continue
        v = float(row[H["slope"]])
        rec("B", "기준%s · %s" % (b, gname), "기울기 %s [%s,%s] · 평균크기 %s · 유닛 %s · rho %s"
            % (row[H["slope"]], row[H["slope_lo"]], row[H["slope_hi"]], row[H["mean_size"]],
               row[H["n_units"]], row[H["spearman"]]), "S14_length.tsv")
        xs.append(gi + (bi - 0.5) * W); ys.append(v)
        el.append(v - float(row[H["slope_lo"]])); eh.append(float(row[H["slope_hi"]]) - v)
    ax.bar(xs, ys, width=W * 0.9, color=["#aab7b8", "#17202a"][bi], label="기준 %s" % b, zorder=3)
    ax.errorbar(xs, ys, yerr=[el, eh], fmt="none", ecolor="#566573", elinewidth=0.5, capsize=1.2, zorder=4)
    for xx, yy in zip(xs, ys):
        ax.text(xx, yy + 0.035, "%.3f" % yy, ha="center", va="bottom", rotation=90, fontsize=CS.FS_ANN - 0.6,   # vertical (the two bars are adjacent so horizontal text overlaps, positions are the same as the original)
                weight="bold", color=["#566573", "#17202a"][bi])
ax.axhline(0, color="#2c3e50", lw=0.4)
# 'remove the right and top borders of j, m and n (mostly keep only left and bottom)'
for _s in ("top", "right"):
    ax.spines[_s].set_visible(False)
_EN = bool(os.environ.get("EN_DICT"))
_ENG = {"근위구획": "Proximal", "가운데구획": "Middle", "원위구획": "Distal"}   # the English tick labels overlapped, so they are shortened (the meaning 'compartment' goes in the Legend)
ax.set_xticks(range(len(GRP))); ax.set_xticklabels([CS.L(g.replace("구획", "\n구획"), _ENG[g]) for g in GRP], linespacing=1.1, ha="center")   # two lines because the cell is narrow (same text)
ax.set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])      # original ticks as is (range from the original formula below)
ax.set_ylim(min(0, min(float(r[H["slope_lo"]]) for r in R)) - 0.06, 1.30)   # the vertical 0.979 text touches the top of the frame -> top margin (old 1.14)
ax.set_ylabel(CS.L("길이 +1 당 늘어난 유닛", "Units per +1 length"))   # axis label on one line and shortened (old: number of units that compartment gains when array length grows by 1 unit)
ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), borderaxespad=0.0, ncol=2, frameon=False, handlelength=1.0,   # legend left = left line of the graph
          borderpad=0.1, columnspacing=0.8)   # one line in the cell margin above the axis (inside the axis it hides the middle bar)
# The note 'only the middle grows, both ends stay near 0' is not drawn on the figure; it is in the legend (it hides bars in an 80 pt cell)
SRC.close()
CS.save(fig, OUT, "cellF2_k_length")
print("기준 %s" % BAS)
