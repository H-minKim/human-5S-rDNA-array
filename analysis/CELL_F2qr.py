# -*- coding: utf-8 -*-
r"""Main Fig 2 cells q and r - panels A and B of FIG_F8m.py, each at cell size (93.1 x 82 pt)
   Data and computation are the same rows as FIG_F8m.py panel A (ZP_gap_shape) and B (F8X2_pairs, F8X2_null, F8X2_sites)
   Only the look changed: figure size, font size (shared by main-text cells), line width, axis label line breaks, annotation position, no panel title
   Content check: CMPART.py (old fig_F8_mechanism axes 0 and 1 ↔ these two figures)
Output: $FIGOUT\cellF2_q_gap.png/svg, cellF2_r_exchange.png/svg, CELL_F2qr_source_values.tsv
"""
import io, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS
TAB = os.environ.get("FIGTAB", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/tables"))
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
_EN = bool(os.environ.get("EN_DICT"))


def read(nm):
    H, rows = None, []
    for L in __import__("taben").open_ko(os.path.join(TAB, nm), encoding="utf8"):  # English tables
        f = L.rstrip("\n").split("\t")
        if not f or f == [""] or f[0].startswith("#"):
            continue
        if H is None:
            H = {c: i for i, c in enumerate(f)}; continue
        rows.append(f)
    return H, rows


SRC = io.open(os.path.join(OUT, "CELL_F2qr_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("panel\t항목\t값\t출처파일\n")


def rec(p, k, v, f):
    SRC.write("%s\t%s\t%s\t%s\n" % (p, k, v, f)); return v


CT, CO, CM = "#c0392b", "#2471a3", "#7d3c98"
KC = {"말단형-말단형": CT, "비말단형-비말단형": CO, "말단형-비말단형": CM}

# ================= q = A : gap shape (same rows as FIG_F8m.py) =================
H, R = read("ZP_gap_shape.tsv")
gg = [int(r[H["gap"]]) for r in R]
obs = [float(r[H["obs_frac"]]) for r in R]
nf = [float(r[H["null_free"]]) for r in R]
n2 = [float(r[H["null_min2"]]) for r in R]
for r in R:
    rec("A", "간격 %s" % r[H["gap"]], "관측 %s (%s 건) · 귀무자유 %s · 귀무최소2 %s"
        % (r[H["obs_frac"]], r[H["obs"]], r[H["null_free"]], r[H["null_min2"]]), "ZP_gap_shape.tsv")
fig, ax = CS.cell_at("q")   # margins and size come from the single cellstyle.ROWS table (aligned plot frames)
W = 0.27; x = np.arange(len(gg))
ax.bar(x - W, obs, W, color=CT, label=CS.L("관측", "Observed"), zorder=3)
ax.bar(x, nf, W, color="#aab7b8", label=CS.L("귀무 ① 자유 배치", "Null (1) free"), zorder=3)   # English legend was longer than the frame, so shortened (meaning given in the Legend)
ax.bar(x + W, n2, W, color="none", edgecolor="#566573", lw=0.5, label=CS.L("귀무 ② 최소 2 칸 띄움", "Null (2) ≥2 apart"), zorder=3)
# 'remove the numbers above the bars' - value text removed (values are in the value table) (old: count '%d' above the observed bars)
i1 = gg.index(1)
# The arrow annotation 'gap 1 has only %s cases (%.4f)' is not drawn on the figure; it is in the legend
rec("A", "간격 1 (그림에서 뺀 주석 · Legend 로)", "간격 1 은 %s 건뿐 (%.4f)" % (R[i1][H["obs"]], obs[i1]), "ZP_gap_shape.tsv")
ax.set_xticks(x); ax.set_xticklabels([str(v) for v in gg])
ax.set_xlabel(CS.L("말단형 사이 간격 (유닛)", "Terminal-unit gap")); ax.set_ylabel(CS.L("간격 몫", "Fraction of gaps"))   # axis label on one line and shorter (old: share taken by that gap)
ax.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])          # original ticks as is (range from the original formula below)
ax.set_ylim(0, max(max(obs), max(nf), max(n2)) * 1.26)
# user 'go with a lower limit of 4.8' : legend text grew to 4.8 pt, so in the Korean edition the last line 'Null (2) at least 2 units apart' overflows the cell on the right by 3.9 pt → in the Korean edition only, the legend's left edge moves 5.6 pt left of the frame's left line (panel letter K sits at -22 to -16 pt from the frame's left, so they do not touch)
ax.legend(loc="lower left", bbox_to_anchor=(0.0 if _EN else -0.11, 1.01), borderaxespad=0.0, frameon=False,   # above the plot, aligned to the left line (so it does not cover the bar numbers)
          fontsize=CS.FS_LEG - (0.8 if _EN else 0.4), handlelength=1.0, borderpad=0.0, labelspacing=0.15)   # borderpad 0.25 → 0 (Helvetica, the Korean-edition legend overflowed the top of the cell, and with a frameless legend only the inner padding shrinks)
CS.legend_fit_above(fig, ax, ax.get_legend())
CS.save(fig, OUT, "cellF2_q_gap")

# ================= r = B : exchange rate (gap ≤2) (same rows as FIG_F8m.py) =================
H, R = read("F8X2_pairs.tsv")
H0, R0 = read("F8X2_null.tsv")
S0 = {r[0]: r[1] for r in R0}
NSITE = S0.get("n_sites", "?")
V = {(r[H["부류"]], int(r[H["간격"]])): r for r in R}
KS = ["말단형-말단형", "비말단형-비말단형", "말단형-비말단형"]
fig, ax = CS.cell_at("r")   # margins and size come from the single cellstyle.ROWS table (aligned plot frames)
xs = np.arange(len(KS))
for ki, k in enumerate(KS):
    r = V.get((k, 0))
    if r is None:
        continue
    v = float(r[H["자카드평균"]])
    rec("B", "%s · 간격 ≤2" % k, "자카드 %s [%s,%s] · 짝 %s · 하플로 %s · 공유자리 %s"
        % (r[H["자카드평균"]], r[H["ci_lo"]], r[H["ci_hi"]], r[H["짝수"]], r[H["하플로수"]],
           r[H["공유자리평균"]]), "F8X2_pairs.tsv")
    ax.bar([ki], [v], 0.56, color=KC[k], zorder=3)
    ax.errorbar([ki], [v], yerr=[[v - float(r[H["ci_lo"]])], [float(r[H["ci_hi"]]) - v]],
                fmt="none", ecolor="#2c3e50", elinewidth=0.5, capsize=1.5, zorder=4)
    ax.text(ki, v + 0.03, "%.3f\n짝 %s" % (v, format(int(r[H["짝수"]]), ",")), ha="center",
            fontsize=CS.FS_ANN - 0.4, color=KC[k], weight="bold", linespacing=1.1)
_ENTK = {"말단형-말단형": "term.–\nterm.", "비말단형-비말단형": "non-t.–\nnon-t.", "말단형-비말단형": "term.–\nnon-t."}   # the cell got wider, so two horizontal lines (old: rotated 90°)   # short ticks for the English edition (meaning given in the Legend)
ax.set_xticks(xs); ax.set_xticklabels([CS.L(k.replace("-", "\n"), _ENTK[k]) for k in KS], linespacing=1.1,
                                      rotation=0, ha="center")   # tick number size = shared by the figure   # in the English edition 'non-terminal-type' is long, so it is laid flat
ax.set_xlim(-0.50, len(KS) - 1 + 0.64)                       # the first bar's value text (pair count) touched the y axis → horizontal margin
ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])              # original ticks as is
ax.set_ylim(0, 1.18)
ax.set_ylabel(CS.L("이웃 교환율 (자카드)", "Exchange rate (Jaccard)"))   # axis label on one line and shorter (old: neighbour exchange rate at sites polymorphic on both sides (Jaccard))
Hs2, Rs2 = read("F8X2_sites.tsv")
_tot = len(Rs2)
_keep = sum(1 for r in Rs2 if r[Hs2["채택"]] == "1")
_DIAG = {"D1096:1", "635CT", "1164GA", "740GA", "DB43:2", "DB16:17", "488GC",
         "DB46:10", "DB58:1", "R1752CG", "R405CG", "R347CG"}
_dk = sum(1 for r in Rs2 if r[Hs2["토큰"]] in _DIAG and r[Hs2["채택"]] == "1")
rec("B", "자리 고르기", "전체 %d 자리 → 양쪽 다형성 %d 자리" % (_tot, _keep), "F8X2_sites.tsv")
rec("B", "옛 진단 12 개 중 남은 것", "%d / %d" % (_dk, len(_DIAG)), "F8X2_sites.tsv")
_r1 = V.get(("말단형-비말단형", 1)); _r2 = V.get(("말단형-비말단형", 2))
if _r1 and _r2:
    rec("B", "말단형-비말단형 간격별", "간격1 %s (짝 %s) · 간격2 %s (짝 %s)"
        % (_r1[H["자카드평균"]], _r1[H["짝수"]], _r2[H["자카드평균"]], _r2[H["짝수"]]), "F8X2_pairs.tsv")
# The top note 'pairs within gap 2 (neighbours) were counted together' is not drawn on the figure; it is in the legend
SRC.close()
CS.save(fig, OUT, "cellF2_r_exchange")
print("A 간격 1 %s · B 자리 %d → %d" % (obs[i1], _tot, _keep))
