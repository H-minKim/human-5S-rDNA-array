# -*- coding: utf-8 -*-
r"""Fig 2 N: permutation test of the excess of AFR-private tokens (2026-10-28 user question 'can it not be tested? check it' and 2026-10-29 'put it in the figure and legend')
   (1) Reproduction: rerun [1] of v111_pop_var_artifact_checks.py with the same random seed (20260923) and compare row by row with X12_POP_VAR_DOWNSAMPLE that the figure reads
   (2) Test: null = population labels are unrelated to variants. Labels are shuffled per individual (the two haplotypes of a person stay together, the number of individuals per population is preserved)
      statistics = median AFR-private token count over K(50) equal-size downsamples, and the median of (AFR - max of the other four populations)
      P = (1 + #{null >= observed}) / (N + 1), N = 1,000, random seed 20261028
   Environment variables (no defaults): RDNA_DATA_ROOT (root that holds v11_readvar), X12_OUT (output folder), X12_PUB (folder with X12_POP_VAR_DOWNSAMPLE.tsv)
   Output: X12_POP_VAR_PERM.tsv"""
import io, os, sys
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8")
T = os.environ["RDNA_DATA_ROOT"].rstrip("/\\") + "/v11_readvar/02_tables/"
OUTD = os.environ["X12_OUT"]; PUB = os.environ["X12_PUB"]
os.makedirs(OUTD, exist_ok=True)
SP5 = ["AFR", "AMR", "EAS", "EUR", "SAS"]
H = pd.read_csv(T + "POP_COHORT.tsv", sep="\t", dtype={"slot": str})
H = H[H.sp.isin(SP5)].copy(); H["hap"] = H["sample"] + "|" + H.slot
L = pd.read_parquet(T + "L2_HAP_UNIT.parquet", columns=["sample", "slot", "varlist"]); L["slot"] = L.slot.astype(str)
L = L.merge(H[["sample", "slot", "hap", "sp"]], on=["sample", "slot"], how="inner")
E = L[["hap", "varlist"]].copy(); E["tok"] = E.varlist.fillna("").str.split(","); E = E.explode("tok"); E = E[E.tok.astype(str).str.len() > 0]
CAR = E.groupby(["tok", "hap"]).size().reset_index()[["tok", "hap"]]
spn = H.groupby("sp").size().to_dict(); nmin = min(spn.values())
NHAP, NIND = len(H), H["sample"].nunique()
print("하플로 %d · 개인 %d · 집단별 하플로 %s · 토큰 %d · 맞추는 수 %d" % (NHAP, NIND, spn, CAR.tok.nunique(), nmin))
# -- (1) production code unchanged (same seed, same order)
rng = np.random.default_rng(20260923)
hap_by_sp = {s: H[H.sp == s].hap.values for s in SP5}
tok_haps = {t: set(g.hap) for t, g in CAR.groupby("tok")}
rows = []
for it in range(200):
    pick = {s: set(rng.choice(hap_by_sp[s], nmin, replace=False)) for s in SP5}
    allp = set().union(*pick.values()); cnt = {s: 0 for s in SP5}; ntok = 0
    for t, hs in tok_haps.items():
        h2 = hs & allp
        if len(h2) < 3: continue
        ins = {s for s in SP5 if h2 & pick[s]}
        if len(ins) == 1: cnt[next(iter(ins))] += 1
        ntok += 1
    rows.append(dict(반복=it, 검정가능토큰=ntok, **cnt))
D = pd.DataFrame(rows); sys.path.insert(0, os.environ.get("V11_CODE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))); import taben  # 
P0 = taben.read_csv_ko(os.path.join(PUB, "X12_POP_VAR_DOWNSAMPLE.tsv"), sep="\t")
same = all((D[c].values == P0[c].values).all() for c in ["반복", "검정가능토큰"] + SP5)
print("① 재현 : 그림 값표와 200 줄 × 7 열 모두 같음 %s" % same)
assert same, "재현 실패 — 멈춘다"
# -- (2) permutation shuffling labels per individual
haps = H.hap.values; hidx = {h: i for i, h in enumerate(haps)}
toks = sorted(tok_haps); M = np.zeros((len(haps), len(toks)), dtype=np.float32)
for j, t in enumerate(toks):
    for h in tok_haps[t]:
        if h in hidx: M[hidx[h], j] = 1.0
ind = H["sample"].values; inds = np.unique(ind); ind_sp = H.groupby("sample").sp.first()
assert (H.groupby("sample").sp.nunique() == 1).all()


def stat(sp_of_hap, rng_, K):
    by = {s: np.where(sp_of_hap == s)[0] for s in SP5}
    m = min(len(v) for v in by.values())
    afr, diff = [], []
    for _ in range(K):
        G = np.zeros((len(haps), 5), dtype=np.float32)
        for k, s in enumerate(SP5): G[rng_.choice(by[s], m, replace=False), k] = 1.0
        c = G.T @ M; tot = c.sum(0); ok = tot >= 3
        ex = [(ok & (c[k] == tot)).sum() for k in range(5)]
        afr.append(ex[0]); diff.append(ex[0] - max(ex[1:]))
    return float(np.median(afr)), float(np.median(diff))


K, N, SEED = 50, 1000, 20261028
rng2 = np.random.default_rng(SEED)
obs = stat(H.sp.values, rng2, K)
labs = ind_sp.loc[inds].values
nul = []
for i in range(N):
    perm = dict(zip(inds, rng2.permutation(labs)))
    nul.append(stat(np.array([perm[x] for x in ind]), rng2, K))
nul = np.array(nul)
out = []
for k, nm in enumerate(["AFR 고유 토큰 중앙값", "AFR − 다른 집단 최대 (중앙값)"]):
    p = (1 + (nul[:, k] >= obs[k]).sum()) / (N + 1)
    out.append(dict(통계량=nm, 관측=obs[k], 귀무중앙=float(np.median(nul[:, k])), 귀무최소=float(nul[:, k].min()),
                    귀무최대=float(nul[:, k].max()), 치환수=N, 내려뽑기수=K, P=p, 가장작은P=1 / (N + 1),
                    하플로타입=NHAP, 개인=NIND, 난수씨=SEED))
    print("② %s : 관측 %.1f · 귀무 %d 회 중앙 %.1f (범위 %.1f–%.1f) · P = %.4g" % (nm, obs[k], N, np.median(nul[:, k]), nul[:, k].min(), nul[:, k].max(), p))
pd.DataFrame(out).to_csv(os.path.join(OUTD, "X12_POP_VAR_PERM.tsv"), sep="\t", index=False, float_format="%.6g")

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('X12_POP_VAR_PERM.tsv',): _taben.to_en_file(_o.path.join(OUTD, _f))
