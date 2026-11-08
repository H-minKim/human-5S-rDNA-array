# -*- coding: utf-8 -*-
r"""Main Fig 2 cell o, only panel B of FIG_X12E.py, at cell size (93.1 x 82 pt)
   Data and computation are the same lines as FIG_X12E.py panel B (X12_POP_VAR_DOWNSAMPLE.tsv, jitter random seed 20261011, same order)
   Only the appearance changes: figure size, font size (shared by main cells), point and line sizes, no panel title (the composite removes it anyway)
   Content check: CMPART.py (old fig_X12e_what axis 1 <-> this figure)
Output : $FIGOUT\cellF2_o_poptoken.png/svg, CELL_F2o_source_values.tsv
"""
import io, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS

TAB = os.environ.get("FIGTAB", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/tables"))
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])


def read(nm):
    H, rows = None, []
    for L in __import__("taben").open_ko(os.path.join(TAB, nm), encoding="utf8"):  # English tables
        if L.startswith("#"): continue
        f = L.rstrip("\n").split("\t")
        if not f or f == [""]: continue
        if H is None: H = {c: i for i, c in enumerate(f)}; continue
        rows.append(f)
    return H, rows


SRC = io.open(os.path.join(OUT, "CELL_F2o_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("panel\t항목\t값\t출처파일\n")
rec = lambda p, k, v, f: SRC.write("%s\t%s\t%s\t%s\n" % (p, k, v, f))
CB, CG, CO = "#2471a3", "#95a5a6", "#E8A33D"

# ===== B: tokens found only in that population, same lines as FIG_X12E.py =====
H2, R2 = read("X12_POP_VAR_DOWNSAMPLE.tsv")
POPS = [c for c in ("AFR", "SAS", "AMR", "EAS", "EUR") if c in H2]
V = {p: np.array([float(r[H2[p]]) for r in R2]) for p in POPS}
NT = np.array([float(r[H2["검정가능토큰"]]) for r in R2])
fig, ax = CS.cell_at("o")   # margins and sizes from the single cellstyle.ROWS table (graph frame alignment)
# 'remove the right and top borders of j, m, n (mostly keep only left and bottom)'
for _s in ("top", "right"):
    ax.spines[_s].set_visible(False)
ORD = sorted(POPS, key=lambda p: -np.median(V[p]))
x = np.arange(len(ORD))
rng = np.random.default_rng(20261011)
for i, p in enumerate(ORD):
    ax.scatter(x[i] + rng.normal(0, 0.07, len(V[p])), V[p], s=0.8, color=CG, alpha=0.45, lw=0, zorder=3)
    m = float(np.median(V[p]))
    ax.plot([x[i] - 0.28, x[i] + 0.28], [m, m], color=CO if p == "AFR" else CB, lw=1.0, zorder=4)
    # median text removed (old ax.annotate("%.1f" % m, ...), values go to the value table via rec below)
    rec("B", "%s 중앙값" % p, "%.1f (200 회 · 사분위 %.1f-%.1f)"
        % (m, np.percentile(V[p], 25), np.percentile(V[p], 75)), "X12_POP_VAR_DOWNSAMPLE.tsv")
rec("B", "검정 가능 토큰 수", "중앙 %.0f (범위 %.0f-%.0f)" % (np.median(NT), NT.min(), NT.max()),
    "X12_POP_VAR_DOWNSAMPLE.tsv")
ax.set_xticks(x); ax.set_xticklabels(ORD)
ax.set_ylabel(CS.L("집단 고유 토큰 수", "Private tokens"))   # axis label on one line, shorter (old: number of tokens found only in that population)
ax.set_ylim(0, max(float(V[p].max()) for p in ORD) * 1.18)
ax.grid(axis="y", linewidth=0.3); ax.set_axisbelow(True)
# read P of the individual-level label shuffling permutation (x12_poptoken_perm.py) from the table and print it
H3, R3 = read("X12_POP_VAR_PERM.tsv")
_pr = [r for r in R3 if r[H3["통계량"]].startswith("AFR 고유")][0]
_p = float(_pr[H3["P"]])
ax.text(0.97, 0.97, "P = %.3f" % _p, transform=ax.transAxes, ha="right", va="top", fontsize=CS.FS_ANN, color="#1b2631")
rec("B", "치환 검정 (개인 단위 표지 섞기)", "AFR 고유 토큰 중앙 관측 %s · 귀무 중앙 %s (범위 %s-%s) · %s 회 · P %.6g (가장 작은 값 %.6g) · 하플로타입 %s · 개인 %s"
    % (_pr[H3["관측"]], _pr[H3["귀무중앙"]], _pr[H3["귀무최소"]], _pr[H3["귀무최대"]], _pr[H3["치환수"]], _p,
       float(_pr[H3["가장작은P"]]), _pr[H3["하플로타입"]], _pr[H3["개인"]]), "X12_POP_VAR_PERM.tsv")
SRC.close()
_xl = ax.get_xlim(); ax.set_xlim(_xl[0], _xl[1] + 0.35)   # right margin because the EUR value text touched the right axis line
CS.save(fig, OUT, "cellF2_o_poptoken")
print("중앙값 %s" % " ".join("%s %.1f" % (p, np.median(V[p])) for p in ORD))
