# -*- coding: utf-8 -*-
r"""FIG F9, population and history: carriage of the terminal-type lineage differs strongly between populations, AFR is lowest

All numbers are read from tables. None are typed into the code.
  A  ZM_population_carrier.tsv (super)  carriage by superpopulation, haplotypes and individuals
  B  ZM_population_carrier.tsv (sub)    subpopulations, ESN 0/41
  C  ZN_source_population.tsv           confound exclusion, the same difference after splitting by assembly source
  D  ZN_depthset_census.tsv             carriage pattern (depth set) x population, smallest sets first (replaced 2026-10-09)
Output : $FIGOUT\fig_F9_population.png/svg, FIG_F9_source_values.tsv
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


def read(nm, keep_hash=False):
    H, rows, notes = None, [], []
    for L in io.open(os.path.join(TAB, nm), encoding="utf8"):
        f = L.rstrip("\n").split("\t")
        if not f or f == [""]:
            continue
        if f[0].startswith("#"):
            notes.append(f); continue
        if H is None:
            H = {c: i for i, c in enumerate(f)}; continue
        rows.append(f)
    return (H, rows, notes) if keep_hash else (H, rows)


SRC = io.open(os.path.join(OUT, "FIG_F9_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("panel\t항목\t값\t출처파일\n")


def rec(p, k, v, f):
    SRC.write("%s\t%s\t%s\t%s\n" % (p, k, v, f)); return v


POPC = {"AFR": "#1a5276", "AMR": "#148f77", "EAS": "#ca6f1e", "EUR": "#7d3c98", "SAS": "#c0392b"}
fig = plt.figure(figsize=(15.0, 9.4))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.24, width_ratios=[1.0, 1.25])

H, R, NOTE = read("ZM_population_carrier.tsv", keep_hash=True)
SUP = [r for r in R if r[H["level"]] == "super"]
SUB = [r for r in R if r[H["level"]] == "sub"]
ST = {}
for f in NOTE:
    ST[f[0].lstrip("#")] = f[1:]

# ================= A: superpopulations =================
ax = fig.add_subplot(gs[0, 0])
SUP = sorted(SUP, key=lambda r: float(r[H["frac_hap"]]))
x = np.arange(len(SUP)); W = 0.38
fh = [float(r[H["frac_hap"]]) for r in SUP]
fi = [float(r[H["frac_indiv"]]) for r in SUP]
for r in SUP:
    rec("A", r[H["group"]], "하플로 %s/%s = %s · 개체 %s/%s = %s"
        % (r[H["n_carrier_hap"]], r[H["n_haplotypes"]], r[H["frac_hap"]],
           r[H["n_indiv_with_any"]], r[H["n_individuals"]], r[H["frac_indiv"]]),
        "ZM_population_carrier.tsv")
ax.bar(x - W / 2, fh, W, color=[POPC[r[H["group"]]] for r in SUP], zorder=3, label="하플로타입")
ax.bar(x + W / 2, fi, W, color=[POPC[r[H["group"]]] for r in SUP], alpha=0.45, zorder=3, label="개체(하나라도)")
for i, r in enumerate(SUP):
    ax.text(i - W / 2, fh[i] + 0.018, "%.3f" % fh[i], ha="center", fontsize=8.6, weight="bold")
    ax.text(i + W / 2, fi[i] + 0.018, "%.3f" % fi[i], ha="center", fontsize=8.6, color="#566573")
ax.set_xticks(x)
ax.set_xticklabels(["%s\n%s / %s" % (r[H["group"]], r[H["n_haplotypes"]], r[H["n_individuals"]])
                    for r in SUP], fontsize=9.0)
ax.set_xlabel("대집단  (아래 숫자 = 하플로타입 수 / 개체 수)", fontsize=9.2)
ax.set_ylim(0, max(fi) * 1.30); ax.set_ylabel("말단형 계통 보유 비율")
ax.legend(fontsize=9, loc="upper left")
chi = ST.get("chi2_hap", [])
perm = ST.get("perm_indiv", [])
if chi:
    txt = "하플로타입 χ² %s · p %s · 자유도 %s" % (chi[0], chi[2], chi[4])
    rec("A", "χ² 검정", txt, "ZM_population_carrier.tsv")
    if perm:
        txt += "\n개체 순열 %s 벌 · 폭 %s · p %s" % (perm[5], perm[1], perm[3])
        rec("A", "개체 순열", "%s 벌 · 폭 %s · p %s" % (perm[5], perm[1], perm[3]),
            "ZM_population_carrier.tsv")
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=8.6, color="#17202a")
ax.set_title("A   집단마다 크게 다르다 — AFR 이 가장 낮다", loc="left", fontsize=12.2)

# ================= B: subpopulations =================
ax = fig.add_subplot(gs[0, 1])
SUB2 = sorted(SUB, key=lambda r: float(r[H["frac_hap"]]))
xb = np.arange(len(SUB2))
POPOF = {"ESN": "AFR", "GWD": "AFR", "MSL": "AFR", "ACB": "AFR", "PUR": "AMR", "PEL": "AMR",
         "GBR": "EUR", "CHS": "EAS", "KHV": "EAS", "JPT": "EAS", "BEB": "SAS", "PJL": "SAS"}
for r in SUB2:
    rec("B", r[H["group"]], "%s/%s = %s (개체 %s)" % (r[H["n_carrier_hap"]], r[H["n_haplotypes"]],
                                                   r[H["frac_hap"]], r[H["n_individuals"]]),
        "ZM_population_carrier.tsv")
vb = [float(r[H["frac_hap"]]) for r in SUB2]
cb = [POPC.get(POPOF.get(r[H["group"]], ""), "#566573") for r in SUB2]
ax.bar(xb, vb, 0.62, color=cb, zorder=3)
for i, r in enumerate(SUB2):
    ax.text(i, vb[i] + 0.012, "%s/%s" % (r[H["n_carrier_hap"]], r[H["n_haplotypes"]]),
            ha="center", fontsize=7.8, color="#17202a")
ax.set_xticks(xb); ax.set_xticklabels([r[H["group"]] for r in SUB2], fontsize=8.6)
ax.set_ylim(0, max(vb) * 1.22); ax.set_ylabel("하플로타입 보유 비율")
hs = [plt.Rectangle((0, 0), 1, 1, fc=POPC[k]) for k in ("AFR", "AMR", "EAS", "EUR", "SAS")]
ax.legend(hs, ["AFR", "AMR", "EAS", "EUR", "SAS"], fontsize=8.6, ncol=5, loc="upper left", frameon=False)
z = SUB2[0]
ax.annotate("%s 는 %s/%s" % (z[H["group"]], z[H["n_carrier_hap"]], z[H["n_haplotypes"]]),
            xy=(0, vb[0]), xytext=(0.6, max(vb) * 0.52), fontsize=9.6, color="#1a5276", weight="bold",
            arrowprops=dict(arrowstyle="->", color="#1a5276", lw=1.2))
ax.set_title("B   하위집단 — 같은 AFR 안에서도 보유가 거의 없다", loc="left", fontsize=12.2)

# ================= C: confound exclusion (assembly source) =================
H2, R2 = read("ZN_source_population.tsv")
ax = fig.add_subplot(gs[1, 0])
SRCS = []
for r in R2:
    if r[H2["source"]] not in SRCS:
        SRCS.append(r[H2["source"]])
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
W = 0.8 / max(len(SRCS), 1)
for si, s in enumerate(SRCS):
    xs, ys, nn = [], [], []
    for pi, pname in enumerate(POPS):
        row = next((r for r in R2 if r[H2["source"]] == s and r[H2["superpop"]] == pname), None)
        if row is None:
            continue
        rec("C", "%s · %s" % (s, pname), "%s/%s = %s" % (row[H2["n_carrier"]], row[H2["n_hap"]],
                                                         row[H2["frac"]]), "ZN_source_population.tsv")
        xs.append(pi + (si - (len(SRCS) - 1) / 2.0) * W); ys.append(float(row[H2["frac"]]))
        nn.append(row[H2["n_hap"]])
    ax.bar(xs, ys, W * 0.9, label="출처 %s" % s, zorder=3,
           color=["#17202a", "#566573", "#aab7b8"][si % 3])
    for xx, yy, n in zip(xs, ys, nn):
        ax.text(xx, yy + 0.015, "n%s" % n, ha="center", fontsize=6.8, color="#566573")
ax.set_xticks(range(len(POPS))); ax.set_xticklabels(POPS, fontsize=9.6)
ax.set_ylim(0, 1.06); ax.set_ylabel("하플로타입 보유 비율")
ax.legend(fontsize=8.8, loc="upper left")
ax.text(0.98, 0.97, "조립 출처를 나눠도 집단 차이가 그대로다\n= 조립 방법 탓이 아니다",
        transform=ax.transAxes, ha="right", va="top", fontsize=9.0, color="#c0392b")
ax.set_title("C   교란 배제 — 조립 출처로 나눠도 같다", loc="left", fontsize=12.2)

# ================= D: depth set size x population (redrawn as bars 2026-10-11) =================
H3, R3 = read("ZN_depthset_census.tsv")
ax = fig.add_subplot(gs[1, 1])
D3 = []
for r in R3:
    cnt = {p_: int(r[H3[p_]]) for p_ in POPS}
    nh = int(r[H3["n_hap"]])
    D3.append((int(r[H3["size"]]), nh, r[H3["depth_set"]], cnt, nh - sum(cnt.values())))
SZ = sorted({d[0] for d in D3})
GRP = [q for q in ("AFR", "SAS", "AMR", "EAS", "EUR") if q in POPS]   # The 2 unlabelled have a denominator of 2 and dominate the axis, excluded here
GC = dict(POPC); GC["미표기"] = "#b3b6b7"


def _n(d, g):
    return d[4] if g == "미표기" else d[3][g]


TOT = {g: sum(_n(d, g) for d in D3) for g in GRP}
MSZ = {g: (sum(d[0] * _n(d, g) for d in D3) / TOT[g] if TOT[g] else float("nan")) for g in GRP}
NHAP = sum(d[1] for d in D3)
NPOP = sum(TOT[g] for g in GRP if g != "미표기")
# Check: per population the size fractions sum to 1, and the haplotype total matches
for g in GRP:
    if TOT[g]:
        _s = sum(sum(_n(d, g) for d in D3 if d[0] == s) / TOT[g] for s in SZ)
        assert abs(_s - 1.0) < 1e-9, "%s 몫 합이 1 이 아니다" % g
UNL = sum(d[4] for d in D3)            # Population unlabelled (not drawn in this panel)
assert sum(TOT.values()) + UNL == NHAP, "하플로타입 총합이 안 맞는다"
print("검산 5 절 : 패턴 %d · 하플로타입 %d (집단 표기 %d · 미표기 %d) · 크기 %s · 집단마다 비율 합 1"
      % (len(D3), NHAP, NPOP, UNL, SZ))
# Denominator = total haplotypes of that population (ZM_population_carrier), unlabelled are counted from the cohort
HZ, RZ = read("ZM_population_carrier.tsv")
TOTALL = {}
for r in RZ:
    if r[HZ["level"]] == "super":
        TOTALL[r[HZ["group"]]] = int(r[HZ["n_haplotypes"]])
TOTALL["미표기"] = TOT.get("미표기", 0)     # (not used in this panel)
for g in GRP:
    assert g in TOTALL, "전체 하플로타입 수가 없다 : %s" % g
    assert TOTALL[g] >= TOT[g], "분모가 분자보다 작다 : %s" % g
print("분모(전체 하플로타입) : " + " · ".join("%s %d (말단형 보유 %d)" % (g, TOTALL[g], TOT[g]) for g in GRP))
W_ = 0.8 / len(GRP)
for j, g in enumerate(GRP):
    if not TOT[g]:
        continue
    v = [sum(_n(d, g) for d in D3 if d[0] == s) / TOTALL[g] for s in SZ]
    xs = [k + (j - (len(GRP) - 1) / 2.0) * W_ for k in range(len(SZ))]
    ax.bar(xs, v, W_ * 0.92, color=GC[g], edgecolor="none", zorder=3)
    for s, x_, a_ in zip(SZ, xs, v):
        rec("D", "%s 크기 %d" % (g, s), "몫 %.4f (%d / %d)"
            % (a_, sum(_n(d, g) for d in D3 if d[0] == s), TOTALL[g]), "ZN_depthset_census.tsv")
import matplotlib.patches as _mp
ax.legend(handles=[_mp.Patch(color=GC[g], label="%s  \uc804\uccb4 %d \u00b7 \ub9d0\ub2e8\ud615 \ubcf4\uc720 %d"
                             % (g, TOTALL[g], TOT[g])) for g in GRP if TOT[g]],
          fontsize=7.0, loc="upper right", frameon=False, labelspacing=0.3, handletextpad=0.4)
ax.set_xticks(range(len(SZ)))
ax.set_xticklabels([str(s) for s in SZ], fontsize=9)
ax.set_xlabel("말단형이 들어앉은 깊이 집합의 크기 (자리 수)")
ax.set_ylabel("그 집단 전체 하플로타입 중 비율")
ax.set_ylim(0, max(sum(_n(d, g) for d in D3 if d[0] == s) / TOTALL[g]
                   for g in GRP if TOT[g] for s in SZ) * 1.55)
ax.grid(axis="y"); ax.set_axisbelow(True)
ax.set_title("D   말단형 보유 패턴의 크기 × 집단  (패턴 %d 개를 크기로 묶었다)" % len(D3),
             loc="left", fontsize=12.2)
# the bottom text 'haplotypes with a population label %d, the 2 unlabelled ...' is removed from the figure and goes to the Legend text (values go to the value table)
rec("D", "그림에서 뺀 글 (Legend 로)", "집단이 표기된 하플로타입 %d · 미표기 2 는 분모가 2 라 이 패널에서 뺐다 (그 둘은 array 유형 그림에 '집단 못 찾음' 으로 보인다)" % NPOP,
    "ZN_depthset_census.tsv")
rec("D", "집단별 평균 집합 크기",
    " · ".join("%s %.3f (n=%d)" % (g, MSZ[g], TOT[g]) for g in GRP if TOT[g]),
    "ZN_depthset_census.tsv")


fig.suptitle("말단형 계통 보유의 집단 분포",
             fontsize=13.0, y=0.975)
SRC.close()
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, "fig_F9_population." + ext), dpi=160, bbox_inches="tight")
print("저장", OUT, "· 대집단 %d · 하위 %d · 출처 %d · 패턴 %d" % (len(SUP), len(SUB2), len(SRCS), len(D3)))
