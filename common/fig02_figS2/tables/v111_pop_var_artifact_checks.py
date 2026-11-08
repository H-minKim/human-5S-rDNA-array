# -*- coding: utf-8 -*-
# v111 (2026-09-23) two artifact checks on the per-population variant results
#   (1) sample size: AFR has the most haplotypes (153), so 'tokens found only in that population' become many by themselves
#      -> downsample the five populations to the same number (smallest population size) 200 times and recount 'private token count' and 'q<0.05 token count'.
#   (2) lineage: the top tokens have almost the same number of carrier haplotypes -> group them by Jaccard to see whether they are one lineage carried together on the same haplotypes.
#      (note: copies form cohesive haplotype lineages, and variant signatures are transmitted as blocks)
#   output POP_VAR_DOWNSAMPLE.tsv and POP_VAR_BLOCKS.tsv in the table folder (02_*, see T), log V111_POP_VAR_CHECK.txt in the log folder (04_*, see LOG)
import io
import numpy as np, pandas as pd
from scipy import stats
OUT = (__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"); T = OUT + "/02_tables"
LOG = io.open(OUT + "/04_logs/V111_POP_VAR_CHECK.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
rng = np.random.default_rng(20260923)
SP5 = ["AFR", "AMR", "EAS", "EUR", "SAS"]

H = pd.read_csv(T + "/POP_COHORT.tsv", sep="\t", dtype={"slot": str})
H = H[H.sp.isin(SP5)].copy()
H["hap"] = H["sample"] + "|" + H.slot
L = pd.read_parquet(T + "/L2_HAP_UNIT.parquet", columns=["sample", "slot", "varlist"])
L["slot"] = L.slot.astype(str)
L = L.merge(H[["sample", "slot", "hap", "sp"]], on=["sample", "slot"], how="inner")
E = L[["hap", "sp", "varlist"]].copy()
E["tok"] = E.varlist.fillna("").str.split(",")
E = E.explode("tok"); E = E[E.tok.astype(str).str.len() > 0]
CAR = E.groupby(["tok", "hap"]).size().reset_index()[["tok", "hap"]]     # haplotype carriage (regardless of unit count)
hap_sp = dict(zip(H.hap, H.sp))
CAR["sp"] = CAR.hap.map(hap_sp)
spn = H.groupby("sp").size().to_dict()
P("[0] 하플로 %d · 집단별 %s · 토큰 %d" % (len(H), spn, CAR.tok.nunique()))

P()
P("[1] 표본 크기를 맞추면 '그 집단에만 있는 토큰' 이 몇 개나 남나")
nmin = min(spn.values())
P("   맞추는 수 = 최소 집단 %s 의 %d 하플로 · 200 회 내려뽑기" % (min(spn, key=spn.get), nmin))
hap_by_sp = {s: H[H.sp == s].hap.values for s in SP5}
tok_haps = {t: set(g.hap) for t, g in CAR.groupby("tok")}
rows = []
for it in range(200):
    pick = {s: set(rng.choice(hap_by_sp[s], nmin, replace=False)) for s in SP5}
    allp = set().union(*pick.values())
    cnt = {s: 0 for s in SP5}; ntok = 0
    for t, hs in tok_haps.items():
        h2 = hs & allp
        if len(h2) < 3: continue
        ins = {s for s in SP5 if h2 & pick[s]}
        if len(ins) == 1:
            s = next(iter(ins))
            if len(h2) >= 3: cnt[s] += 1
        ntok += 1
    rows.append(dict(반복=it, 검정가능토큰=ntok, **cnt))
D = pd.DataFrame(rows)
D.to_csv(T + "/POP_VAR_DOWNSAMPLE.tsv", sep="\t", index=False)
P("   전용 토큰 수(보유 하플로 3+ · 200 회 중앙값) : %s" % {s: float(D[s].median()) for s in SP5})
P("   같은 수로 맞춘 뒤에도 AFR 이 가장 많은 반복 : %d / 200" % int((D[["AFR"] + [s for s in SP5 if s != "AFR"]].idxmax(axis=1) == "AFR").sum()))
P("   AFR 전용 토큰 수 범위 %d–%d · 2 위 집단 중앙 %.0f"
  % (D.AFR.min(), D.AFR.max(), float(np.median([sorted([r[s] for s in SP5])[-2] for _, r in D.iterrows()]))))

P()
P("[2] 상위 AFR 토큰들이 같은 하플로에 함께 실려 있나 (자카드)")
V = pd.read_csv(T + "/POP_VAR_TOKEN.tsv", sep="\t")
top = V[(V.보유_q < 0.05) & (V.보유최대차집단 == "AFR")].sort_values("보유_p").head(60).토큰.astype(str).tolist()
S = {t: tok_haps[t] for t in top if t in tok_haps}
ks = list(S)
J = np.zeros((len(ks), len(ks)))
for i, a in enumerate(ks):
    for j, b in enumerate(ks):
        J[i, j] = len(S[a] & S[b]) / max(1, len(S[a] | S[b]))
# single-linkage clustering (Jaccard 0.5 or more)
lab = list(range(len(ks)))
def find(x):
    while lab[x] != x: lab[x] = lab[lab[x]]; x = lab[x]
    return x
for i in range(len(ks)):
    for j in range(i + 1, len(ks)):
        if J[i, j] >= 0.5: lab[find(i)] = find(j)
grp = {}
for i, t in enumerate(ks): grp.setdefault(find(i), []).append(t)
P("   상위 %d 토큰 → 자카드 0.5 단일연결 군집 %d 개" % (len(ks), len(grp)))
rows = []
for g in sorted(grp.values(), key=len, reverse=True):
    hs = set.intersection(*[S[t] for t in g]) if len(g) > 1 else S[g[0]]
    un = set.union(*[S[t] for t in g])
    sp_cnt = pd.Series([hap_sp[h] for h in un]).value_counts().to_dict()
    rows.append(dict(토큰수=len(g), 예시토큰=" ".join(g[:6]), 합집합하플로=len(un), 교집합하플로=len(hs),
                     개체=len({h.split("|")[0] for h in un}), 집단분포=str(sp_cnt)))
    if len(g) >= 2 or len(rows) <= 8:
        P("     %2d 토큰 · 하플로 합 %3d 교 %3d · 개체 %3d · %s · %s"
          % (len(g), len(un), len(hs), len({h.split("|")[0] for h in un}), sp_cnt, " ".join(g[:8])))
BL = pd.DataFrame(rows)
BL.to_csv(T + "/POP_VAR_BLOCKS.tsv", sep="\t", index=False)
P("   해석 : 군집 하나 = 한 array 계통이 통째로 특정 집단에 몰린 것 · 토큰 하나하나가 독립 사건이 아니다")

P()
P("[3] 친족 제외하고 상위 결과가 남나 (자식 하플로 뺌)")
H2 = H[~H.kid_in_cohort]
c2 = CAR[CAR.hap.isin(set(H2.hap))]
spn2 = H2.groupby("sp").size().to_dict()
res = []
for t in top[:20]:
    g = c2[c2.tok == t]
    car = g.groupby("sp").size().to_dict()
    a = car.get("AFR", 0); b = spn2["AFR"] - a
    c = sum(car.get(x, 0) for x in SP5 if x != "AFR"); d = sum(spn2[x] for x in SP5 if x != "AFR") - c
    orr, pv = stats.fisher_exact([[a, b], [c, d]])
    res.append((t, a, c, pv))
P("   친족 제외 뒤 상위 20 토큰 : 여전히 p<0.001 인 것 %d / 20" % sum(p < 1e-3 for _, _, _, p in res))
P("     예 : %s" % " · ".join("%s AFR %d/타 %d p %.1g" % r for r in res[:5]))
LOG.close()

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('POP_VAR_DOWNSAMPLE.tsv',): _taben.to_en_file(_o.path.join(T, _f))
