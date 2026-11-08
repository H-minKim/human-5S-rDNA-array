# -*- coding: utf-8 -*-
r"""Main Fig 2 cell b: only panel a of P08_FIG_4bars.py fig1, at cell size (177.4 x 76.5 pt)
   Data and computation are the same lines as P08_FIG_4bars.py fig1 (reads only Q1_results and X1_population_control of OUT_P08)
   Only the appearance changed: figure size, font size (common to main-text cells), font (Malgun Gothic in the English version too, so cells match), no panel letter (the composite removes it anyway)
   English text is the EN dictionary of the P08 code as is (English run = when the environment variable EN_DICT is set)
   Content check: CMPART.py (old FIG_P08_1_scope axes 0 vs this figure)
Output: $FIGOUT\cellF2_b_scope.png/svg, CELL_F2b_source_values.tsv
"""
import os, sys, io
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS

OUTP = os.environ.get("P08_OUT", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/handover_P08_homogenisation_2026-10-10/tables_OUT_P08"))   # read only
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
EN_RUN = bool(os.environ.get("EN_DICT"))
BLUE, GREY, ORANGE = "#2a78d6", "#8a8f98", "#eb6834"
EN = {"닮은 상대 보유": "Units with a match", "완전 동일": "Identical", "5 변이 이내": "Within 5 variants",
      "array 내": "Within\narray", "상동": "Homologous", "무관\n같은 세부집단": "Unrelated\nsame pop.",
      "무관\n전체": "Unrelated\nall"}          # only the lines this panel uses, copied as is from the P08_FIG_4bars.py EN dictionary


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


hdr, r13 = read_sec("Q1_results.tsv", "Q1-13")
tw = {x[0]: x for x in r13}
hdrD, rD = read_sec("X1_population_control.tsv", "D")
dd = {x[0]: x for x in rD}
labs = ["array 내", "상동", "무관\n같은 세부집단", "무관\n전체"]
ks = ["상동 (같은 개인)", "무관 · 같은 세부집단", "무관 · 전체"]
f0 = [100 * float(tw["T1 array내"][7])] + [100 * float(dd[k][3]) for k in ks]
f5 = [100 * float(tw["T1 array내"][9])] + [100 * float(dd[k][7]) for k in ks]
with io.open(os.path.join(OUT, "CELL_F2b_source_values.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("무리\t완전동일_pct\t5변이이내_pct\t출처\n")
    for l_, a_, b_ in zip(labs, f0, f5):
        g.write("%s\t%.6f\t%.6f\tQ1_results[Q1-13] · X1_population_control[D]\n" % (l_.replace("\n", " "), a_, b_))

fig, ax = CS.cell_at("b")   # margins and size from the single cellstyle.ROWS table (plot frame alignment)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
x = np.arange(len(labs))
w = 0.30
b1 = ax.bar(x - w / 2, f0, w, color=GREY, edgecolor="none", label=T("완전 동일"), zorder=3)
b2 = ax.bar(x + w / 2, f5, w, color=BLUE, edgecolor="none", label=T("5 변이 이내"), zorder=3)
for bars, v in ((b1, f0), (b2, f5)):
    for bb, vv in zip(bars, v):
        ax.text(bb.get_x() + bb.get_width() / 2, vv + 2.0, "%.1f" % vv, ha="center", va="bottom", fontsize=CS.FS_ANN)
ax.set_xticks(x)
ax.set_xticklabels([T(s) for s in labs], linespacing=1.2)
ax.set_ylim(0, 112)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_ylabel(T("닮은 상대 보유") + " (%)")
ax.legend(frameon=False, loc="upper right", handlelength=1.0, borderpad=0.1)
# removed the old ax.grid(axis="y", color="#dddddd", linewidth=0.3, zorder=0)
ax.set_axisbelow(True)
CS.save(fig, OUT, "cellF2_b_scope")
print("완전 동일 %s · 5 변이 이내 %s" % (["%.1f" % v for v in f0], ["%.1f" % v for v in f5]))
