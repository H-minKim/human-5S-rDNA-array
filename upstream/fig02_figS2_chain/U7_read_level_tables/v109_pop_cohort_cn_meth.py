# -*- coding: utf-8 -*-
# v109 (2026-09-23) array features by ancestry (population): (1) cohort definition (2) CN (unit count N) (3) methylation zone metrics
#   data: L2_HAP_UNIT.parquet in the table folder (02_*) (combined haplotype x unit table, includes pop/super_pop)
#   labels: existing pop.panel (1kGP 2504) + 1kGP_3202_ped_population_2026-09-23.txt in the reference folder (05_*) (extended 3202, public ftp)
#          filled only after first confirming that the roster agrees with the 170 existing labels (170/170).
#   Confounder control: depth, relatives (child and parent both included), multiple testing (BH), label-shuffling null
#   Interpretation caveat: basecaller version information is not available per sample (a remaining confounder for absolute beta comparison) -> focus on state fractions (f_hi).
#   output POP_COHORT.tsv, POP_CN.tsv, POP_METH_ZONE.tsv in the table folder (02_*), log V109_POP_CN_METH.txt in the log folder (04_*)
import io, os
import numpy as np, pandas as pd
from scipy import stats
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"; T = OUT + "/02_tables"
LOG = io.open(OUT + "/04_logs/V109_POP_CN_METH.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
rng = np.random.default_rng(20260923)
SP5 = ["AFR", "AMR", "EAS", "EUR", "SAS"]

P("[0] 자료·라벨")
L = pd.read_parquet(T + "/L2_HAP_UNIT.parquet")
L["slot"] = L.slot.astype(str)
P("   L2 행 %d · 하플로타입 %d · 개체 %d · 자료원 %s"
  % (len(L), L.groupby(["sample", "slot"]).ngroups, L["sample"].nunique(), sorted(L.data_source.astype(str).unique())))
PED = pd.read_csv(OUT + "/05_reference/1kGP_3202_ped_population_2026-09-23.txt", sep=r"\s+")
m_sp = dict(zip(PED.SampleID, PED.Superpopulation)); m_pop = dict(zip(PED.SampleID, PED.Population))
kid = {s for s, f, mo in zip(PED.SampleID, PED.FatherID, PED.MotherID) if str(f) != "0" or str(mo) != "0"}
par = {s: (str(f), str(mo)) for s, f, mo in zip(PED.SampleID, PED.FatherID, PED.MotherID)}

S = L.groupby("sample").agg(pop_old=("pop", "first"), sp_old=("super_pop", "first")).reset_index()
chk = S[S.sp_old.notna()]
agree = int(sum(m_sp.get(s) == v for s, v in zip(chk["sample"], chk.sp_old)))
P("   기존 라벨 %d 명 ↔ 3202 명부 일치 %d (불일치 %d) — 불일치면 멈춤" % (len(chk), agree, len(chk) - agree))
assert agree == len(chk)
S["sp"] = [v if pd.notna(v) else m_sp.get(s) for s, v in zip(S["sample"], S.sp_old)]
S["pop"] = [v if pd.notna(v) else m_pop.get(s) for s, v in zip(S["sample"], S.pop_old)]
S["label_src"] = np.where(S.sp_old.notna(), "panel2504", np.where(S.sp.notna(), "ped3202", "없음"))
P("   라벨 원천 : %s" % S.label_src.value_counts().to_dict())
P("   라벨 못 붙인 개체 : %s" % sorted(S[S.sp.isna()]["sample"]))

# per-haplotype table
H = L.groupby(["sample", "slot"]).agg(
    N=("N", "first"), depth=("depth", "first"), gate_rate=("gate_rate", "first"), cov_frac=("cov_frac", "first"),
    n_unit=("ad", "size"), n_read_slot=("n_read_slot", "first"), f_hap=("f_hap", "first"),
    Bp=("Bp", "first"), Bd=("Bd", "first"), n_low=("n_low", "first"), N_other=("N_other", "first")).reset_index()
H = H.merge(S[["sample", "pop", "sp", "label_src"]], on="sample", how="left")
inco = set(H["sample"])
H["kid_in_cohort"] = [s in kid and any(p in inco for p in par.get(s, ())) for s in H["sample"]]
H.to_csv(T + "/POP_COHORT.tsv", sep="\t", index=False, float_format="%.4f")
P("   코호트 : 하플로 %d · 개체 %d · 라벨 있는 하플로 %d" % (len(H), H["sample"].nunique(), H.sp.notna().sum()))
P("   친족(부모가 코호트 안에 있는 자식) 개체 %d · 하플로 %d" % (H[H.kid_in_cohort]["sample"].nunique(), H.kid_in_cohort.sum()))
G = H[H.sp.isin(SP5)].copy()
P("   집단별 하플로/개체 :")
for sp, g in G.groupby("sp"):
    P("     %-4s 하플로 %3d · 개체 %3d · 깊이 중앙 %5.1f · 친족 하플로 %d" % (sp, len(g), g["sample"].nunique(), g.depth.median(), g.kid_in_cohort.sum()))

def kw(d, col, by="sp"):
    gs = [g[col].dropna().values for _, g in d.groupby(by)]
    gs = [x for x in gs if len(x) >= 3]
    if len(gs) < 2: return np.nan, np.nan
    h, p = stats.kruskal(*gs)
    return h, p

def perm_p(d, col, by="sp", n=2000):
    obs = kw(d, col, by)[0]
    if not np.isfinite(obs): return np.nan
    v = d[[col, by]].dropna(); x = v[col].values; lab = v[by].values.copy()
    c = 0
    for _ in range(n):
        rng.shuffle(lab)
        h = kw(pd.DataFrame({col: x, by: lab}), col, by)[0]
        c += (h >= obs)
    return (c + 1) / (n + 1)

P()
P("[1] CN (array 유닛 수 N) — 하플로타입 단위")
rows = []
for sp, g in G.groupby("sp"):
    rows.append(dict(집단=sp, 하플로=len(g), 중앙=g.N.median(), 평균=g.N.mean(), sd=g.N.std(), q10=g.N.quantile(.1), q90=g.N.quantile(.9), 최소=g.N.min(), 최대=g.N.max()))
CN = pd.DataFrame(rows)
P(CN.round(2).to_string(index=False))
h, p = kw(G, "N"); pp = perm_p(G, "N")
P("   집단 간 차 (Kruskal-Wallis) H %.2f · p %.4f · 라벨 뒤섞기 p %.4f" % (h, p, pp))
pr = []
for i, a in enumerate(SP5):
    for b in SP5[i + 1:]:
        x, y = G[G.sp == a].N.values, G[G.sp == b].N.values
        if len(x) >= 5 and len(y) >= 5:
            u, pv = stats.mannwhitneyu(x, y)
            pr.append(dict(쌍="%s vs %s" % (a, b), n1=len(x), n2=len(y), 중앙1=np.median(x), 중앙2=np.median(y), 차=np.median(x) - np.median(y), p=pv))
PR = pd.DataFrame(pr)
if len(PR):
    PR["q"] = np.minimum(1, PR.p * len(PR))
    P("   쌍별 (Mann-Whitney · Bonferroni q):")
    P(PR.round(4).to_string(index=False))
P("   깊이와 N 의 관계 (스피어만) rho %+.3f p %.3g — 깊이가 N 을 만들지 않는지 확인"
  % stats.spearmanr(G.depth, G.N)[:2])
P("   깊이 자체가 집단 간 다른가 : KW p %.4f" % kw(G, "depth")[1])
G1 = G[~G.kid_in_cohort]
P("   친족 제외(자식 뺌 · 하플로 %d) : KW p %.4f · 중앙 %s" % (len(G1), kw(G1, "N")[1], G1.groupby("sp").N.median().round(1).to_dict()))
I = G.groupby(["sample", "sp"]).N.mean().reset_index()
P("   개체 단위(두 하플로 평균 · 개체 %d) : KW p %.4f · 중앙 %s" % (len(I), kw(I, "N")[1], I.groupby("sp").N.median().round(1).to_dict()))
CN.to_csv(T + "/POP_CN.tsv", sep="\t", index=False, float_format="%.3f")
if len(PR): PR.to_csv(T + "/POP_CN_PAIRS.tsv", sep="\t", index=False, float_format="%.5f")

P()
P("[2] 메틸 — 구역별 (6 구역 · 하플로타입 단위)")
ZON = [("M1 근1-3", "ap", 1, 3), ("M2 근4-10", "ap", 4, 10), ("M3 근11-25", "ap", 11, 25),
       ("D2 원11-22", "ad", 11, 22), ("D1 원1-10", "ad", 1, 10)]
L2 = L.merge(S[["sample", "sp"]], on="sample", how="left")
L2 = L2[L2.sp.isin(SP5)]
recs = []
for nm, col, a, b in ZON:
    z = L2[(L2[col] >= a) & (L2[col] <= b)]
    g = z.groupby(["sample", "slot", "sp"]).agg(f_hi=("f_hi", "mean"), b_mean=("b_mean", "mean"),
                                                b_lo=("b_lo_mean", "mean"), b_hi=("b_hi_mean", "mean"),
                                                n_unit=("ad", "size")).reset_index()
    g["구역"] = nm
    recs.append(g)
Z = pd.concat(recs, ignore_index=True)
Z = Z.merge(H[["sample", "slot", "depth", "N", "kid_in_cohort"]], on=["sample", "slot"], how="left")
out = []
for (nm, met) in [(z[0], m) for z in ZON for m in ("f_hi", "b_mean", "b_lo", "b_hi")]:
    d = Z[(Z["구역"] == nm)][[met, "sp", "depth", "N"]].dropna()
    if d[met].nunique() < 5: continue
    h, p = kw(d, met)
    med = d.groupby("sp")[met].median()
    out.append(dict(구역=nm, 지표=met, n=len(d), **{s: round(float(med.get(s, np.nan)), 4) for s in SP5},
                    최대차=round(float(med.max() - med.min()), 4), H=round(h, 2), p=p))
ZR = pd.DataFrame(out)
ZR["q"] = ZR.p * len(ZR)
ZR["q"] = ZR["q"].clip(upper=1)
P("   (지표 : f_hi = 고메틸 유닛 비율 · b_mean 유닛 β · b_lo/b_hi = 저/고 상태 안의 β)")
P(ZR.sort_values("p").to_string(index=False))
sig = ZR[ZR.q < 0.05]
P("   q<0.05 인 칸 %d / %d" % (len(sig), len(ZR)))
# haplotype-level metrics (boundary, number of switched-off units)
for met, lab in (("Bp", "근위 경계 위치"), ("Bd", "원위 경계 위치"), ("n_low", "저메틸 유닛 수"), ("f_hap", "하플로 고메틸 비율")):
    d = G[[met, "sp"]].dropna()
    if len(d) < 30: continue
    h, p = kw(d, met)
    P("   %-12s(%s) : 중앙 %s · KW p %.4f" % (met, lab, d.groupby("sp")[met].median().round(2).to_dict(), p))
# depth control: repeated within depth quintiles
P("   깊이 5 분위 안에서 다시 (f_hi · 각 구역) :")
Z["dq"] = pd.qcut(Z.depth, 5, labels=False, duplicates="drop")
for nm, _, _, _ in ZON:
    ps = []
    for q, g in Z[Z["구역"] == nm].groupby("dq"):
        if g.sp.nunique() >= 3: ps.append(kw(g, "f_hi")[1])
    if ps: P("     %-10s 분위별 p %s" % (nm, " ".join("%.3f" % x for x in ps)))
Z.to_csv(T + "/POP_METH_ZONE.tsv", sep="\t", index=False, float_format="%.4f")
ZR.to_csv(T + "/POP_METH_ZONE_TEST.tsv", sep="\t", index=False, float_format="%.5g")
P()
P("저장 : 02_tables/POP_COHORT.tsv · POP_CN.tsv · POP_CN_PAIRS.tsv · POP_METH_ZONE.tsv · POP_METH_ZONE_TEST.tsv")
LOG.close()
