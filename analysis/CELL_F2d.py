# -*- coding: utf-8 -*-
r"""Main Fig 2 cell d, only panel b of P08_FIG.py fig4, at cell size (96 x 82 pt)
   Data are the same lines as P08_FIG.py fig4 panel b (reads only X4_lag2_and_patch[B] and X9_patch_floor from OUT_P08)
   Only the appearance changes: figure size, font size (shared by main cells), font (Malgun Gothic also for the English edition), no panel letter, ticks pinned to the original values
   English text is the EN dictionary of the P08 code as is (English run = when the environment variable EN_DICT is set)
   Content check: CMPART.py (old FIG_P08_4_tract axis 1 <-> this figure)
Output : $FIGOUT\cellF2_d_tract.png/svg, CELL_F2d_source_values.tsv
"""
import os, sys, io
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS

OUTP = os.environ.get("P08_OUT", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/handover_P08_homogenisation_2026-10-10/tables_OUT_P08"))   # read only
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
EN_RUN = bool(os.environ.get("EN_DICT"))
BLUE, GREY, ORANGE = "#2a78d6", "#8a8f98", "#eb6834"
EN = {"유닛 거리 (copy)": "Separation (repeat units)", "패치 길이 추정 (bp)": "Patch length (bp)",
      "교환 없는 바닥 (무관 개체)": "No-exchange floor (unrelated)"}


def T(s):
    return EN.get(s, s) if EN_RUN else s


def read_sec(path, sec):
    rows, on, hdr = [], False, None
    for Ln in __import__("taben").open_ko(os.path.join(OUTP, path), encoding="utf8"):
        if Ln.startswith("["):
            on = Ln.startswith("[" + sec + "]")
            if on:
                hdr = Ln.rstrip("\n").split("\t")[1:]
            continue
        if on and Ln.startswith("\t"):
            rows.append(Ln.rstrip("\n").split("\t")[1:])
    return hdr, rows


hdrB2, rB2 = read_sec("X4_lag2_and_patch.tsv", "B")
pl = [(int(x[0]), float(x[3]), float(x[4]), float(x[5])) for x in rB2 if x[0].isdigit()]
pg = np.array([p[0] for p in pl])
pm = np.array([p[1] for p in pl])
plo = np.array([p[2] for p in pl])
phi_ = np.array([p[3] for p in pl])
flo = None
for Ln in __import__("taben").open_ko(os.path.join(OUTP, "X9_patch_floor.tsv"), encoding="utf8"):
    f_ = Ln.rstrip().split(chr(9))
    if f_[0] == "unrelated_all":
        flo = float(f_[4])
with io.open(os.path.join(OUT, "CELL_F2d_source_values.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("거리\t패치길이\t하\t상\t출처\n")
    for row in zip(pg, pm, plo, phi_):
        g.write("%d\t%.6f\t%.6f\t%.6f\tX4_lag2_and_patch[B]\n" % row)
    g.write("바닥(무관 전체)\t%s\t\t\tX9_patch_floor[unrelated_all]\n" % flo)

fig, ax = CS.cell_at("d")   # margins and sizes from the single cellstyle.ROWS table (graph frame alignment)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.fill_between(pg, plo, phi_, color="#fbdccf", edgecolor="none", zorder=2)
ax.plot(pg, pm, "-o", color=ORANGE, ms=1.4, lw=0.7, zorder=4)
ax.set_xlabel(T("유닛 거리 (copy)"))
ax.set_ylabel(T("패치 길이 추정 (bp)"))
if flo:
    ax.axhline(flo, color=GREY, lw=0.5, ls=(0, (3, 2)), zorder=3)
    ax.text(25, flo - 22, T("교환 없는 바닥 (무관 개체)"), ha="right", va="top", fontsize=CS.FS_ANN, color="#333333")
ax.set_xticks([-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0])        # original ticks as is (visible range 0 to 25)
ax.set_yticks([0.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0])
ax.set_xlim(-0.2, 26.2)
ax.set_ylim(0, 660)
# horizontal grid removed (old ax.grid(axis="y", color="#dddddd", linewidth=0.3)), the dashed no-exchange floor is data and is kept
ax.set_axisbelow(True)
CS.save(fig, OUT, "cellF2_d_tract")
print("패치 길이 거리 1 %.1f · 바닥 %s" % (pm[0], flo))
