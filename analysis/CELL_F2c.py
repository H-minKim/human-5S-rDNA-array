# -*- coding: utf-8 -*-
r"""Main Fig 2 cell c - only panel a of the P08_FIG.py fig2, at cell size (101 x 82 pt)
   Data are the same rows as panel a of P08_FIG.py fig2 (reads only XB2_band_all and XB2_band_noodd from OUT_P08)
   Only the look changed: figure size, font size (shared by main-text cells), font (Malgun Gothic also in the English edition), no panel letter, ticks pinned to the original values
   English text is the EN dictionary of the P08 code as is (English run = when env var EN_DICT is set)
   Content check: CMPART.py (old FIG_P08_2_distance axis 0 ↔ this figure)
Output: $FIGOUT\cellF2_c_distance.png/svg, CELL_F2c_source_values.tsv
"""
import os, sys, io
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS

OUTP = os.environ.get("P08_OUT", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/handover_P08_homogenisation_2026-10-10/tables_OUT_P08"))   # read only
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
EN_RUN = bool(os.environ.get("EN_DICT"))
BLUE, GREY, ORANGE = "#2a78d6", "#8a8f98", "#eb6834"
EN = {"귀무 띠": "Global envelope", "관측": "Observed", "발산 사분위 제외": "Divergent quartile removed",
      "유닛 거리 (copy)": "Separation (repeat units)", "쌍별 평균 차이 (bp)": "Mean pairwise diff. (bp)"}


def T(s):
    return EN.get(s, s) if EN_RUN else s


BA = np.load(os.path.join(OUTP, "XB2_band_all.npz"))
BN = np.load(os.path.join(OUTP, "XB2_band_noodd.npz"))
with io.open(os.path.join(OUT, "CELL_F2c_source_values.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("층\t거리\t관측\t귀무_하\t귀무_상\t출처\n")
    for nm, B, src in (("전체", BA, "XB2_band_all.npz"), ("발산 사분위 제외", BN, "XB2_band_noodd.npz")):
        for k, o, lo, hi in zip(B["KS"], B["obs"], B["lo"], B["hi"]):
            g.write("%s\t%g\t%.6f\t%.6f\t%.6f\t%s\n" % (nm, k, o, lo, hi, src))

fig, ax = CS.cell_at("c")   # margins and size come from the single cellstyle.ROWS table (aligned plot frames)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.fill_between(BA["KS"], BA["lo"], BA["hi"], color="#e2e5e9", edgecolor="none", zorder=2, label=T("귀무 띠"))
ax.fill_between(BN["KS"], BN["lo"], BN["hi"], color="#fbdccf", edgecolor="none", zorder=2)
ax.plot(BA["KS"], BA["obs"], "-o", color=BLUE, ms=1.4, lw=0.7, zorder=4, label=T("관측"))
ax.plot(BN["KS"], BN["obs"], "-o", color=ORANGE, ms=1.4, lw=0.7, zorder=4, label=T("발산 사분위 제외"))
ax.set_xlabel(T("유닛 거리 (copy)"))
ax.set_ylabel(T("쌍별 평균 차이 (bp)"))
ax.set_xticks([0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5])          # original ticks as is (visible ones are 2.5 to 15.0)
ax.set_yticks([4.0, 6.0, 8.0, 10.0, 12.0, 14.0])                    # original ticks as is (visible ones are 6 to 12)
ax.set_xlim(0.4, 15.6)
ax.set_ylim(5.4, 12.8)
ax.legend(frameon=False, loc="center right", handlelength=1.0, borderpad=0.1)
# horizontal grid lines removed (old: ax.grid(axis="y", color="#dddddd", linewidth=0.3))
ax.set_axisbelow(True)
CS.save(fig, OUT, "cellF2_c_distance")
print("관측 거리 1 %.3f · 15 %.3f" % (BA["obs"][0], BA["obs"][-1]))
