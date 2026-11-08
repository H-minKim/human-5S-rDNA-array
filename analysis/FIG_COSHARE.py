# -*- coding: utf-8 -*-
r"""FIG COSHARE - two units that share a rare variant also share other variants (new 2026-10-16: replaces supplementary e)
  Data: COSHARE_PERHAP.tsv (produced by COSHARE_BSTRAT.py, H1 557, middle window 0.2-0.8, 533 haplotypes with 23 or more window units)
    A = rare variant (cohort unit frequency < 1%), pair = two units in the same array
    cls 0 = pairs that share A, cls 1 = pairs at the same distance where only one unit carries A
    value = fraction of other variants B, carried by at least one of the two units, that are in both units, by |B - A| bin
  Panel A: B is rare (< 1%)   Panel B: B is common (≥ 10%)
  bars = haplotype cluster bootstrap 95% (500 times). Bins at or below 100 bp are not drawn (a single event spanning several positions would be mixed in)
Output: $FIGOUT\fig_coshare.png/svg, FIG_COSHARE_source_values.tsv
"""
import io, os, collections
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
# The English version uses this table without a hook (2026-10-16, the translation hook folder in scratchpad disappeared, so only this figure translates itself)
EN = os.environ.get("FIG_LANG", "ko") == "en"
TR = {"희귀 변이를 함께 가진 두 유닛은 다른 변이도 함께 가진다": "Two units that share a rare variant also share other variants",
      "같은 거리, A 를 한쪽만 가진 두 유닛": "Same distance, only one unit carries A",
      "희귀 변이 A 를 함께 가진 두 유닛": "Both units carry rare variant A",
      "다른 변이를 두 유닛이 모두 가진 비율": "Fraction of other variants carried by both units",
      "A 에서 그 변이까지의 거리 (bp)": "Distance from A to the other variant (bp)",
      "A   다른 변이가 희귀할 때": "A   Other variant rare",
      "B   다른 변이가 흔할 때": "B   Other variant common"}
T_ = lambda s: TR[s] if EN else s
BINS = ["101-250", "251-500", "501-1000", "1001-2240"]
XLAB = ["101–250", "251–500", "501–1,000", "1,001–2,240"]
STR = {"B 희귀 (<1%)": 0, "B 흔함 (>=10%)": 1}

L = [l.rstrip("\n").split("\t") for l in __import__("taben").open_ko(os.path.join(TAB, "COSHARE_PERHAP.tsv"), encoding="utf8")
     if not l.startswith("#")]
H = {c: i for i, c in enumerate(L[0])}
haps = sorted({(r[H["sample"]], r[H["slot"]]) for r in L[1:]})
HI = {h: i for i, h in enumerate(haps)}
A = np.zeros((len(haps), 2, len(BINS), 2, 2))
for r in L[1:]:
    st = STR.get(r[H["B_stratum"]]); b = r[H["bin_bp"]]
    if st is None or b not in BINS:
        continue
    A[HI[(r[H["sample"]], r[H["slot"]])], st, BINS.index(b), int(r[H["cls"]])] = (int(r[H["n_shared"]]), int(r[H["n_cand"]]))
assert len(haps) == 541, len(haps)   # the alignment set is the full 541 (was 531, 10 slots missing from the FASTA map)
rng = np.random.default_rng(20261016); NB = 500


def rates(idx, si):
    T = A[idx, si].sum(0)
    return T[:, :, 0] / T[:, :, 1]


SRC = io.open(os.path.join(OUT, "FIG_COSHARE_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("panel\tB_stratum\tbin_bp\tcls\tn_cand\tshare\tlo\thi\tfold_together_over_onesided\tsource\n")
C = {0: "#e67e22", 1: "#7f8c8d"}
W = 414.4 / 72.0
fig, axs = plt.subplots(1, 2, figsize=(W, W * 0.44), constrained_layout=True)
for ax, si, pan, ttl, ylim in ((axs[0], 0, "A", T_("A   다른 변이가 희귀할 때"), (0, 0.50)),
                               (axs[1], 1, "B", T_("B   다른 변이가 흔할 때"), (0.6, 1.0))):
    r = rates(np.arange(len(haps)), si)
    bs = np.array([rates(rng.integers(0, len(haps), len(haps)), si) for _ in range(NB)])
    lo, hi = np.percentile(bs, 2.5, axis=0), np.percentile(bs, 97.5, axis=0)
    x = np.arange(len(BINS))
    for cls, lab, dx in ((0, T_("희귀 변이 A 를 함께 가진 두 유닛"), -0.06), (1, T_("같은 거리, A 를 한쪽만 가진 두 유닛"), 0.06)):
        ax.errorbar(x + dx, r[:, cls], yerr=[r[:, cls] - lo[:, cls], hi[:, cls] - r[:, cls]],
                    fmt="-o", color=C[cls], ms=4, lw=1.4, capsize=2.5, label=lab)
        for i in range(len(x)):
            ax.text(x[i] + dx + 0.09, r[i, cls], "%.3f" % r[i, cls], ha="left", va="center", fontsize=5.6, color=C[cls])
            ncand = int(A[:, si, i, cls, 1].sum())
            SRC.write("%s\t%s\t%s\t%d\t%d\t%.5f\t%.5f\t%.5f\t%.3f\tCOSHARE_PERHAP.tsv\n"
                      % (pan, [k for k, v in STR.items() if v == si][0], BINS[i], cls, ncand, r[i, cls], lo[i, cls], hi[i, cls],
                         r[i, 0] / r[i, 1]))
    if si == 0:
        for i in range(len(x)):
            ax.text(x[i], ylim[1] * 0.985, "×%.0f" % (r[i, 0] / r[i, 1]), ha="center", va="top", fontsize=6.4, color="#17202a")
    ax.set_xticks(x); ax.set_xticklabels(XLAB, fontsize=6.4); ax.set_xlim(-0.35, len(x) - 0.35)
    ax.set_xlabel(T_("A 에서 그 변이까지의 거리 (bp)"), fontsize=7)
    ax.set_ylabel(T_("다른 변이를 두 유닛이 모두 가진 비율"), fontsize=7)
    ax.set_ylim(*ylim); ax.tick_params(labelsize=6.6)
    ax.grid(axis="y", lw=0.4, alpha=0.5); ax.set_axisbelow(True)
    ax.set_title(ttl, loc="left", fontsize=9)
SRC.close()
h_, l_ = axs[0].get_legend_handles_labels()
fig.legend(h_, l_, fontsize=6.6, frameon=False, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.07))
fig.suptitle(T_("희귀 변이를 함께 가진 두 유닛은 다른 변이도 함께 가진다"), fontsize=10)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "fig_coshare." + ext), dpi=300, bbox_inches="tight")
r0 = rates(np.arange(len(haps)), 0)
print("하플로타입 %d · A(희귀 B) 배수 %s" % (len(haps), " ".join("%.1f" % (r0[i, 0] / r0[i, 1]) for i in range(len(BINS)))))
