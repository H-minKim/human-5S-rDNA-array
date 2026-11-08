# -*- coding: utf-8 -*-
# v110 (2026-09-23) array features by ancestry (population): (3) unit variants, population-specific and frequency differences
#   data: varlist (variant tokens per unit) of L2_HAP_UNIT.parquet in the table folder (02_*), population labels from POP_COHORT.tsv of v109
#   two metrics: (1) carrier = whether the token is present in at least one unit of that haplotype array
#                  (2) within-array frequency = units carrying the token / units of that haplotype   (as a ratio, since copy numbers differ)
#   tests: carrier = Fisher (that population vs the rest), frequency = Kruskal-Wallis (5 populations), BH q
#   Confounders: depth (detection power for rare variants), N (unit count), relatives, multiple testing -> top results are rechecked within depth quantiles and without relatives
#   output POP_VAR_TOKEN.tsv, POP_VAR_PRIVATE.tsv, POP_VAR_HAPLOAD.tsv in the table folder (02_*), log V110_POP_VAR.txt in the log folder (04_*)
import io
import numpy as np, pandas as pd
from scipy import stats
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"; T = OUT + "/02_tables"
LOG = io.open(OUT + "/04_logs/V110_POP_VAR.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
rng = np.random.default_rng(20260923)
SP5 = ["AFR", "AMR", "EAS", "EUR", "SAS"]
MIN_HAP = 5          # a token must be present in at least this many haplotypes to be tested
MAX_HAPFRAC = 0.95   # nearly fixed tokens (population differences cannot be seen) are removed

H = pd.read_csv(T + "/POP_COHORT.tsv", sep="\t", dtype={"slot": str})
H = H[H.sp.isin(SP5)].copy()
L = pd.read_parquet(T + "/L2_HAP_UNIT.parquet", columns=["sample", "slot", "ad", "ap", "varlist", "n_var", "n_var_rare", "n_var_common"])
L["slot"] = L.slot.astype(str)
L = L.merge(H[["sample", "slot", "sp", "N", "depth", "kid_in_cohort"]], on=["sample", "slot"], how="inner")
P("[0] 자료 : 하플로 %d · 개체 %d · 유닛 %d · 집단별 하플로 %s"
  % (L.groupby(["sample", "slot"]).ngroups, L["sample"].nunique(), len(L), H.groupby("sp").size().to_dict()))

# expand unit x token
E = L[["sample", "slot", "sp", "varlist"]].copy()
E["tok"] = E.varlist.fillna("").str.split(",")
E = E.explode("tok")
E = E[E.tok.astype(str).str.len() > 0]
P("   (하플로,유닛,토큰) 쌍 %d · 서로 다른 토큰 %d" % (len(E), E.tok.nunique()))

nunit = L.groupby(["sample", "slot"]).size().rename("n_unit")
hap_sp = H.set_index(["sample", "slot"]).sp
# haplotype x token: number of carrier units
C = E.groupby(["sample", "slot", "tok"]).size().rename("n_hit").reset_index()
C = C.merge(nunit, on=["sample", "slot"]).merge(H[["sample", "slot", "sp"]], on=["sample", "slot"])
C["frac"] = C.n_hit / C.n_unit
NH = len(H)
tok_hap = C.groupby("tok")["sample"].size().rename("n_hap")
keep = tok_hap[(tok_hap >= MIN_HAP) & (tok_hap <= MAX_HAPFRAC * NH)].index
P("   검정 대상 토큰 %d (하플로 %d+ 보유 · 보유율 %.0f%% 이하) · 전체 %d"
  % (len(keep), MIN_HAP, MAX_HAPFRAC * 100, C.tok.nunique()))

spn = H.groupby("sp").size().to_dict()
rows = []
Cg = {t: g for t, g in C[C.tok.isin(keep)].groupby("tok")}
for t, g in Cg.items():
    car = g.groupby("sp")["sample"].size().to_dict()
    # carrier: Fisher per population (that population vs the rest)
    best_sp, best_p, best_or = None, 1.0, np.nan
    for s in SP5:
        a = car.get(s, 0); b = spn[s] - a
        c = sum(car.get(x, 0) for x in SP5 if x != s); d = sum(spn[x] for x in SP5 if x != s) - c
        orr, pv = stats.fisher_exact([[a, b], [c, d]])
        if pv < best_p: best_sp, best_p, best_or = s, pv, orr
    # within-array frequency: haplotypes without the token are filled with 0, then KW
    fr = {s: [] for s in SP5}
    have = dict(zip(zip(g["sample"], g["slot"]), g.frac))
    for (sm, sl), s in hap_sp.items():
        fr[s].append(have.get((sm, sl), 0.0))
    hk, pk = stats.kruskal(*[np.array(fr[s]) for s in SP5])
    med = {s: float(np.mean(fr[s])) for s in SP5}
    rows.append(dict(토큰=t, 보유하플로=int(tok_hap[t]),
                     **{"보유_" + s: car.get(s, 0) for s in SP5},
                     **{"빈도_" + s: round(med[s], 4) for s in SP5},
                     보유최대차집단=best_sp, 보유_p=best_p, 보유_OR=best_or, 빈도_H=hk, 빈도_p=pk))
V = pd.DataFrame(rows)
def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); r = np.empty(len(p))
    r[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.clip(r, 0, 1)
V["보유_q"] = bh(V.보유_p); V["빈도_q"] = bh(V.빈도_p)
V = V.sort_values("빈도_p")
V.to_csv(T + "/POP_VAR_TOKEN.tsv", sep="\t", index=False, float_format="%.5g")
P()
P("[1] 집단 간 빈도 차 (array 안 빈도 · KW · BH q)")
P("   q<0.05 토큰 %d / %d · q<0.01 %d" % ((V.빈도_q < 0.05).sum(), len(V), (V.빈도_q < 0.01).sum()))
cols = ["토큰", "보유하플로"] + ["빈도_" + s for s in SP5] + ["빈도_p", "빈도_q"]
P(V.head(15)[cols].to_string(index=False))
P()
P("[2] 보유 여부 차 (Fisher · BH q) — q<0.05 토큰 %d" % (V.보유_q < 0.05).sum())
cols2 = ["토큰", "보유하플로"] + ["보유_" + s for s in SP5] + ["보유최대차집단", "보유_OR", "보유_p", "보유_q"]
P(V.sort_values("보유_p").head(15)[cols2].to_string(index=False))

P()
P("[3] 한 집단에만 있는 토큰 (보유 하플로 3 개 이상 · 다른 집단 0)")
priv = []
for t, g in C.groupby("tok"):
    car = g.groupby("sp")["sample"].size().to_dict()
    ins = [s for s in SP5 if car.get(s, 0) > 0]
    if len(ins) == 1 and car[ins[0]] >= 3:
        s = ins[0]
        priv.append(dict(토큰=t, 집단=s, 보유하플로=car[s], 개체=g["sample"].nunique(),
                         집단하플로수=spn[s], 보유율=round(car[s] / spn[s], 4),
                         array안빈도중앙=round(float(g.frac.median()), 4), 최대=round(float(g.frac.max()), 4)))
PV = pd.DataFrame(priv)
if len(PV):
    PV = PV.sort_values(["보유하플로"], ascending=False)
    PV.to_csv(T + "/POP_VAR_PRIVATE.tsv", sep="\t", index=False, float_format="%.5g")
    P("   %d 개 · 집단별 %s" % (len(PV), PV.groupby("집단").size().to_dict()))
    P(PV.head(20).to_string(index=False))
else:
    P("   없음")

P()
P("[4] 하플로당 변이 부담 (유닛당 평균 변이 수 · 희귀 변이 수)")
B = L.groupby(["sample", "slot", "sp"]).agg(n_var=("n_var", "mean"), n_rare=("n_var_rare", "mean"),
                                            n_common=("n_var_common", "mean"), N=("N", "first"),
                                            depth=("depth", "first"), kid=("kid_in_cohort", "first")).reset_index()
B.to_csv(T + "/POP_VAR_HAPLOAD.tsv", sep="\t", index=False, float_format="%.4f")
for met in ("n_var", "n_rare", "n_common"):
    h, p = stats.kruskal(*[g[met].values for _, g in B.groupby("sp")])
    P("   %-9s 집단중앙 %s · KW p %.4g" % (met, B.groupby("sp")[met].median().round(3).to_dict(), p))
P("   깊이와 희귀변이 수 rho %+.3f · N 과 rho %+.3f (검출력 교란 확인)"
  % (stats.spearmanr(B.depth, B.n_rare)[0], stats.spearmanr(B.N, B.n_rare)[0]))
B2 = B[~B.kid]
P("   친족 제외(하플로 %d) : n_rare KW p %.4g" % (len(B2), stats.kruskal(*[g.n_rare.values for _, g in B2.groupby("sp")])[1]))
P("   깊이 5 분위 안에서 n_rare : p %s"
  % " ".join("%.3f" % stats.kruskal(*[g.n_rare.values for _, g in q.groupby("sp") if len(g) >= 3])[1]
             for _, q in B.assign(dq=pd.qcut(B.depth, 5, labels=False)).groupby("dq")))
P()
P("저장 : 02_tables/POP_VAR_TOKEN.tsv · POP_VAR_PRIVATE.tsv · POP_VAR_HAPLOAD.tsv")
LOG.close()
