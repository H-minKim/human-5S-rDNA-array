# -*- coding: utf-8 -*-
r"""FIG TRIOPOS - trio transmission compared down to the unit position (new on 2026-10-13: replaces supplementary a)
  Data: <trio_current>/T1_pair_summary.tsv (H1 557, position-wise comparison after unit alignment, best alignment direction and shift per pair)
  transmitted     = not unrelated + pairs with position-wise mean Jaccard >= 0.95 (user approved 2026-11-05. The next value is 0.838, and the pairs are the same 5 as at 0.99) (hg002|2 <- hg003|2, whose relation carries the 'excluded' tag, is dropped)
  non-transmitted = the other haplotype of the same parent of each transmitted pair (best alignment only)
  unrelated       = the 150 unrelated pairs (best alignment only)
  A  fraction of units whose variant set is exactly identical at each position (frac_exact)
  B  mean position-wise Jaccard (mean_jaccard)
Output: $FIGOUT\fig_trio_unitmatch.png/svg, FIG_TRIOPOS_source_values.tsv
"""
import io, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as _fm
for _p in ((__import__("os").environ.get("RDNA_FONT_DIR", __import__("os").path.join(__import__("os").environ.get("SystemRoot", ""), "Fonts")) + "/malgun.ttf"), (__import__("os").environ.get("RDNA_FONT_DIR", __import__("os").path.join(__import__("os").environ.get("SystemRoot", ""), "Fonts")) + "/gulim.ttc")):
    if os.path.exists(_p):
        _fm.fontManager.addfont(_p)
        plt.rcParams["font.family"] = _fm.FontProperties(fname=_p).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["svg.fonttype"] = "none"
TAB = os.environ.get("FIGTAB", (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout/tables"))
OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
SRCF = os.path.join(TAB, "trio_current_reference", "T1_pair_summary.tsv")

L = [l.rstrip("\n").split("\t") for l in __import__("taben").open_ko(SRCF, encoding="utf8")]
R = [dict(zip(L[0], r)) for r in L[1:]]
best = {}
for r in R:
    k = (r["child"], r["child_slot"], r["parent"], r["parent_slot"])
    if k not in best or float(r["mean_jaccard"]) > float(best[k]["mean_jaccard"]):
        best[k] = r
TR = [r for r in best.values() if r["relation"] != "무관" and float(r["mean_jaccard"]) >= 0.95
      and "제외" not in r["relation"]]
NT = []
for r in TR:
    o = "2" if r["parent_slot"] == "1" else "1"
    b = best.get((r["child"], r["child_slot"], r["parent"], o))
    if b is not None: NT.append(b)
UN = [r for r in best.values() if r["relation"] == "무관"]
assert (len(TR), len(NT), len(UN)) == (5, 5, 150), (len(TR), len(NT), len(UN))
assert all(r["child_N"] == r["parent_N"] for r in TR), "전달 짝의 유닛 수가 다르다"
n_ex = sum(round(float(r["frac_exact"]) * int(r["n_pairs"])) for r in TR)
n_al = sum(int(r["n_pairs"]) for r in TR)

SRC = io.open(os.path.join(OUT, "FIG_TRIOPOS_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("group\tchild\tchild_slot\tparent\tparent_slot\tchild_N\tparent_N\tn_aligned\tmean_jaccard\tfrac_exact\tsource\n")
for g, V in (("전달", TR), ("비전달", NT), ("무관", UN)):
    for r in V:
        SRC.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tT1_pair_summary.tsv\n"
                  % (g, r["child"], r["child_slot"], r["parent"], r["parent_slot"], r["child_N"], r["parent_N"],
                     r["n_pairs"], r["mean_jaccard"], r["frac_exact"]))
SRC.write("# 전달 5 짝 유닛 %d 개 중 자리별 완전 일치 %d (%.4f)\n" % (n_al, n_ex, n_ex / n_al))
SRC.close()

C = {"전달": "#2471a3", "비전달": "#e67e22", "무관": "#95a5a6"}
G = [("전달", TR), ("비전달", NT), ("무관", UN)]
W = 414.4 / 72.0
fig, axs = plt.subplots(1, 2, figsize=(W, W * 0.42), constrained_layout=True)
rng = np.random.default_rng(20261013)
for ax, col, ylab, title in ((axs[0], "frac_exact", "자리별로 변이가 완전히 같은 유닛의 비율", "A   자리까지 같은 유닛"),
                             (axs[1], "mean_jaccard", "자리별 변이 자카드 (유닛 평균)", "B   자리별 닮음")):
    for i, (g, V) in enumerate(G):
        v = np.array([float(r[col]) for r in V])
        x = i + (rng.uniform(-0.17, 0.17, len(v)) if len(v) > 10 else np.linspace(-0.12, 0.12, len(v)))
        ax.scatter(x, v, s=14 if len(v) > 10 else 22, color=C[g], alpha=0.55 if len(v) > 10 else 0.9,
                   edgecolors="none", zorder=3)
        md = float(np.median(v))
        ax.plot([i - 0.28, i + 0.28], [md, md], color="#17202a", lw=1.4, zorder=4)
        ax.text(i + 0.31, md, "%.3f" % md, va="center", ha="left", fontsize=6.4, color="#17202a")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["전달", "비전달", "무관"], fontsize=7)
    for i, (g, V) in enumerate(G):   # n is kept apart from the Korean text (so no number enters the English-edition dictionary keys)
        ax.text(i, -0.13, "n = %d" % len(V), ha="center", va="top", fontsize=6.2, color="#566573",
                transform=ax.get_xaxis_transform())
    ax.set_ylabel(ylab, fontsize=7.2)
    ax.set_ylim(-0.03, 1.08); ax.set_xlim(-0.55, 2.75)
    ax.tick_params(labelsize=7)
    ax.grid(axis="y", lw=0.4, alpha=0.5); ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=10)
axs[0].text(0, 1.035, "%d / %d" % (n_ex, n_al), ha="center", va="bottom", fontsize=6.4, color=C["전달"])
fig.suptitle("삼자 전달, 유닛 자리까지 맞춘 비교", fontsize=10.5)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "fig_trio_unitmatch." + ext), dpi=300, bbox_inches="tight")
print("전달 %d · 비전달 %d · 무관 %d · 전달 유닛 %d 중 완전 일치 %d" % (len(TR), len(NT), len(UN), n_al, n_ex))
