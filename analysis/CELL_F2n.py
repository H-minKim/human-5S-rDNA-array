# -*- coding: utf-8 -*-
r"""Main Fig 2 cell n: only panel A of FIG_X12T.py, at cell size (93.1 x 82 pt)
   Data and computation are the same lines as panel A of FIG_X12T.py (X12_TYPE_GROUPS.tsv, reference 'majority', 5 or more haplotypes)
   Only the appearance changed: figure size, font size (common to main-text cells), legend position (inside the cell below the axes), no panel title (the composite removes it anyway)
   Content check: CMPART.py (old fig_X12t_type_groups axes 0 vs this figure)
Output: $FIGOUT\cellF2_n_typepop.png/svg, CELL_F2n_source_values.tsv
"""
import io, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS

TAB = os.environ.get("FIGTAB", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/tables"))
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
MINHAP = 5


def read(nm):
    H, rows = None, []
    for L in __import__("taben").open_ko(os.path.join(TAB, nm), encoding="utf8"):  # English tables
        if L.startswith("#"): continue
        f = L.rstrip("\n").split("\t")
        if not f or f == [""]: continue
        if H is None: H = {c: i for i, c in enumerate(f)}; continue
        rows.append(f)
    return H, rows


SRC = io.open(os.path.join(OUT, "CELL_F2n_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("panel\t항목\t값\t출처파일\n")
rec = lambda p, k, v, f: SRC.write("%s\t%s\t%s\t%s\n" % (p, k, v, f))
H, R = read("X12_TYPE_GROUPS.tsv")
G = [r for r in R if r[H["기준"]] == "다수형" and int(r[H["하플로"]]) >= MINHAP]
G.sort(key=lambda r: -int(r[H["하플로"]]))
POP = ["AFR", "AMR", "EAS", "EUR", "SAS"]
PC = {"AFR": "#D4A017", "AMR": "#1E8449", "EAS": "#1F618D", "EUR": "#922B21", "SAS": "#6C3483"}
PC["미표기"] = "#b3b6b7"

# ===== A: population composition (including unlabelled), the same lines as FIG_X12T.py =====
fig, ax = CS.cell_at("n")   # margins and size from the single cellstyle.ROWS table (plot frame alignment)
# 'remove the right and top borders of j, m and n (mostly only left and bottom)'
for _s in ("top", "right"):
    ax.spines[_s].set_visible(False)
x = np.arange(len(G))
bot = np.zeros(len(G))
NH = np.array([int(r[H["하플로"]]) for r in G], float)
UN = np.array([int(r[H["하플로"]]) - sum(int(r[H["n_" + p]]) for p in POP) for r in G], float)
for p in POP + ["미표기"]:
    v = (np.array([int(r[H["n_" + p]]) for r in G], float) if p != "미표기" else UN) / NH
    ax.bar(x, v, 0.76, bottom=bot, color=PC[p], edgecolor="none", zorder=3,
           label=p if p != "미표기" else "집단 못 찾음")
    bot += v
assert np.allclose(bot, 1.0), "구성 합이 1 이 아니다"
for i, r in enumerate(G):
    ax.text(x[i], 1.015, r[H["하플로"]], ha="center", va="bottom", fontsize=CS.FS_ANN - 0.8, color="#566573")
    rec("A", "무리 %s" % r[H["무리"]],
        "하플로타입 %s · " % r[H["하플로"]]
        + " · ".join("%s %s" % (p, r[H["n_" + p]]) for p in POP)
        + " · 미표기 %d" % UN[i], "X12_TYPE_GROUPS.tsv")
rec("A", "미표기 합", "%d / %d" % (UN.sum(), NH.sum()), "X12_TYPE_GROUPS.tsv")
ax.set_xticks(x); ax.set_xticklabels([r[H["무리"]] for r in G])   # tick number size = common to the figure
ax.set_xlabel("array 유형 무리")   # the bracket '(numbers above are haplotype counts)' goes to the Legend text
ax.set_ylabel("집단 구성")
# vertical ticks keep the original values, the original Korean version (0.2 step) and English version (0.25 step) differed (automatic ticks), each version follows its own
ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0, 1.25] if os.environ.get("EN_DICT") else [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
ax.set_ylim(0, 1.09)
ax.legend(frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0.0, -0.24), borderaxespad=0.0,   # legend left = left line of the plot frame
          handlelength=1.0, columnspacing=0.7, handletextpad=0.3, fontsize=CS.FS_LEG - 0.4)
SRC.close()
CS.save(fig, OUT, "cellF2_n_typepop")
print("무리 %d · 하플로타입 %d · 미표기 %d" % (len(G), NH.sum(), UN.sum()))
