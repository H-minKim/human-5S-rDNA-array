# -*- coding: utf-8 -*-
r"""P08 XB2: fixes the distance definition error of XB
   XB error: divergent units were removed and the units were renumbered before counting distance => 'distance 1' mixed in pairs that were originally 2-3 apart
   Fix (same convention as X5): keep the original numbering, and use only pairs whose two ends are both non-divergent
   The null uses the same convention: shuffle the order with labels following the units, and apply the same mask after shuffling
   Output: for each of the two sets (all / divergent excluded), the observed curve and its own null band"""
import os, sys, io, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
SEED, NSIM, K = 20261010, 2499, 15
rng = np.random.default_rng(SEED)
import v11_io as G
from scipy import stats

U, _, _ = G.load("h1")
VPU = sum(len(t[4]) for v in U.values() for t in v) / sum(len(v) for v in U.values())
print("§0 자 : var_per_unit = %.6f" % VPU)
assert abs(VPU - 6.866) < 0.01
NCOPY = {("%s" % k[0], "%s" % k[1]): len(v) for k, v in U.items()}
with open(os.path.join(HERE, "OUT_P04", "S3_matrices.pkl"), "rb") as f:
    S3 = pickle.load(f)
MATS, APS = S3["mats"], S3["aps"]
KS = np.arange(1, K + 1)
PACK = []
for s in sorted(MATS):
    aps = APS[s]
    N = NCOPY.get(s, 0)
    sel = [i for i, ap in enumerate(aps) if 0.2 * N < ap <= 0.8 * N]
    if len(sel) < 22:
        continue
    ii = np.array(sel)
    M = MATS[s][np.ix_(ii, ii)].astype(np.int64)
    n = len(sel)
    med = int(np.argmin(M.sum(1)))
    div = M[med].astype(float)
    PACK.append((M, n, div > np.percentile(div, 75)))
print("haplotype %d · 유닛 %d · 발산 %d"
      % (len(PACK), sum(x[1] for x in PACK), sum(int(x[2].sum()) for x in PACK)))


def erl_env(obs, sims, alpha=0.05):
    A = np.vstack([np.asarray(obs)[None, :], sims])
    m, d = A.shape
    R = np.empty((m, d))
    for i in range(d):
        c = A[:, i]
        R[:, i] = np.minimum(stats.rankdata(c, method="min"), stats.rankdata(-c, method="min"))
    Rs = np.sort(R, axis=1)
    kk = [tuple(Rs[j]) for j in range(m)]
    p = sum(1 for z in kk if z <= kk[0]) / m
    order = sorted(range(m), key=lambda j: kk[j])
    keep = order[int(np.floor(alpha * m)):]
    return p, A[keep].min(0), A[keep].max(0)


def run(use_all):
    num = np.zeros(K)
    den = np.zeros(K)
    SN = np.zeros((NSIM, K))
    SD = np.zeros((NSIM, K))
    allsum = allcnt = 0.0
    for M, n, odd in PACK:
        ok = np.ones(n, bool) if use_all else ~odd
        # observed, original numbering kept, only pairs with both ends ok
        for ki, k in enumerate(KS):
            if n - k <= 0:
                continue
            i = np.arange(n - k)
            m = ok[i] & ok[i + k]
            if m.any():
                num[ki] += M[i[m], i[m] + k].sum()
                den[ki] += int(m.sum())
        tri = M[np.triu_indices(n, 1)]
        if use_all:
            allsum += tri.sum()
            allcnt += len(tri)
        else:
            w = np.where(ok)[0]
            t2 = M[np.ix_(w, w)][np.triu_indices(len(w), 1)]
            allsum += t2.sum()
            allcnt += len(t2)
        # null, shuffle the order with labels following the units
        O = np.argsort(rng.random((NSIM, n)), axis=1)
        okp = ok[O]                                  # (NSIM, n)
        for ki, k in enumerate(KS):
            if n - k <= 0:
                continue
            a = O[:, :-k]
            b = O[:, k:]
            m = okp[:, :-k] & okp[:, k:]
            v = np.where(m, M[a, b], 0)
            SN[:, ki] += v.sum(1)
            SD[:, ki] += m.sum(1)
    obs = num / den
    sims = SN / SD
    p, lo, hi = erl_env(obs, sims)
    return obs, sims, p, lo, hi, allsum / allcnt, den


print("\n" + "=" * 78)
rows = []
for tag, ua, key in (("전체", True, "all"), ("발산 사분위 제외", False, "noodd")):
    obs, sims, p, lo, hi, allmean, den = run(ua)
    out = (np.where((obs < lo) | (obs > hi))[0] + 1).tolist()
    nm = sims.mean(0)
    print("[%s]  쌍 수 : 거리 1 %d · 거리 15 %d" % (tag, den[0], den[-1]))
    print("    그 집합 안 모든 쌍의 평균 차이 = %.4f bp" % allmean)
    print("    귀무 곡선 평균 : 거리 1 %.4f · 거리 15 %.4f (차 %.4f — 평평)" % (nm[0], nm[-1], abs(nm[0] - nm[-1])))
    print("    귀무 띠 폭 평균 %.4f bp · 거리 1 [%.4f, %.4f]" % (np.mean(hi - lo), lo[0], hi[0]))
    print("    관측 : 거리 1 %.4f · 거리 15 %.4f" % (obs[0], obs[-1]))
    print("    거리 1 이 귀무보다 %.4f bp 낮다 = %.2f%%" % (nm[0] - obs[0], 100 * (nm[0] - obs[0]) / nm[0]))
    print("    전역 p = %.5g · 띠 밖 거리 %s\n" % (p, out))
    rows.append((tag, allmean, float(nm[0]), float(lo[0]), float(hi[0]), float(obs[0]), float(obs[-1]),
                 float(np.mean(hi - lo)), p, out))
    np.savez_compressed(os.path.join(OUT, "XB2_band_%s.npz" % key),
                        KS=KS, obs=obs, lo=lo, hi=hi, null_mean=nm, p=p, n_pair=den)
print("=" * 78)
print("두 집합의 귀무 띠 높이 : %.4f bp ↔ %.4f bp — 같은 띠를 쓰면 안 된다" % (rows[0][2], rows[1][2]))
g = io.open(os.path.join(OUT, "XB2_band_explained.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 XB2 · Fig2a 귀무 띠의 정체 · 거리 정의 고침 (원래 번호 유지) · nsim %d · seed %d\n" % (NSIM, SEED))
g.write("set\tmean_all_pairs\tnull_mean_lag1\tlo_lag1\thi_lag1\tobs_lag1\tobs_lag15\tband_width\tp_erl\tout_lags\n")
for r in rows:
    g.write("%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6g\t%s\n"
            % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], ",".join(map(str, r[9]))))
g.close()
print("\n산출 :", os.path.join(OUT, "XB2_band_explained.tsv"))
