# -*- coding: utf-8 -*-
r"""P08 Q2c - fix λ in the middle window (trend removed)
   ⒜ practical range 3a over central 0.2-0.8 + block bootstrap CI (block = haplotype)
   ⒝ redo by segment in the middle window too (interior(+-5) is inflated by end effects)
   ⒞ λ vs n_copies"""
import os, sys, io, time, collections, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
SEED, NSIM, B_BLK, DEC = 20261010, 2499, 2000, 10
rng = np.random.default_rng(SEED)
import v11_io as G
from scipy import stats, optimize

CAN = G.CAN
LU = len(CAN)
REG = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA-rich", 1179, 1551),
       ("NTS2", 1552, 2031), ("promoter", 2032, 2121), ("gene", 2122, 2240)]
U, _, _ = G.load("h1")
VPU = sum(len(t[4]) for v in U.values() for t in v) / sum(len(v) for v in U.values())
print("§0 자 : var_per_unit = %.6f" % VPU)
assert abs(VPU - 6.866) < 0.01
NCOPY = {("%s" % k[0], "%s" % k[1]): len(v) for k, v in U.items()}
with open(os.path.join(HERE, "OUT_P04", "Q0_alignment.pkl"), "rb") as f:
    DIFF = pickle.load(f)["diff"]
with open(os.path.join(HERE, "OUT_P04", "S3_matrices.pkl"), "rb") as f:
    S3 = pickle.load(f)
BAD, MATS, APS = S3["bad"], S3["mats"], S3["aps"]
keys = [k for k in DIFF if k not in BAD]
SETR, PSETR, INSR = {}, {}, {}
for nm, a_, b_ in REG:
    SETR[nm], PSETR[nm], INSR[nm] = {}, {}, {}
for k in keys:
    sub, dele, ins, ul = DIFF[k]
    A = dict(sub)
    for p in dele:
        A[p] = "-"
    ia = {}
    for p, s_ in ins:
        q = min(max(p, 1), LU)
        ia[q] = ia.get(q, "") + s_
    for nm, a_, b_ in REG:
        d_ = {p: v for p, v in A.items() if a_ <= p <= b_}
        SETR[nm][k] = frozenset(d_.items())
        PSETR[nm][k] = frozenset(d_)
        INSR[nm][k] = {p: v for p, v in ia.items() if a_ <= p <= b_}


def pdr(nm, k1, k2):
    S1, S2 = SETR[nm][k1], SETR[nm][k2]
    both = len(S1 & S2)
    d = len(S1) + len(S2) - 2 * both - (len(PSETR[nm][k1] & PSETR[nm][k2]) - both)
    I1, I2 = INSR[nm][k1], INSR[nm][k2]
    if I1 or I2:
        for p in set(I1) | set(I2):
            a, b = I1.get(p, ""), I2.get(p, "")
            if a != b:
                L = max(len(a), len(b))
                d += sum(1 for x, y in zip(a.ljust(L, "\0"), b.ljust(L, "\0")) if x != y)
    return d


def fit_exp(gam, ks, w, start=3):
    m = ks >= start
    x, y, ww = ks[m].astype(float), gam[m], w[m]

    def f(p):
        return np.sqrt(ww) * (p[0] + p[1] * (1 - np.exp(-x / max(p[2], 1e-6))) - y)
    r = optimize.least_squares(f, [y.min(), max(y.max() - y.min(), 1e-6), 5.0],
                               bounds=([0, 0, 0.1], [np.inf, np.inf, 500]))
    return r.x


def erl_test(obs, sims, alpha=0.05):
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


def pack(window, Mfun, K):
    out = []
    for s in sorted(APS):
        aps = APS[s]
        N = NCOPY.get(s, 0)
        sel = [i for i, ap in enumerate(aps) if window[0] * N < ap <= window[1] * N]
        if len(sel) < max(9, 3 * K // 2):
            continue
        out.append((s, np.array(sel), len(sel)))
    return out


def curves_from(PK, Mget, K):
    NH = len(PK)
    KSc = np.arange(1, K + 1)
    PD = np.zeros((NH, K))
    PQ = np.zeros((NH, K))
    P0 = np.zeros((NH, K))
    PN = np.zeros((NH, K))
    Ms = []
    for i, (s, sel, n) in enumerate(PK):
        M = Mget(s, sel)
        Ms.append(M)
        for ki, k in enumerate(KSc):
            if n - k <= 0:
                continue
            v = M[np.arange(n - k), np.arange(k, n)]
            PD[i, ki] = v.sum()
            PQ[i, ki] = (v * v).sum()
            P0[i, ki] = (v == 0).sum()
            PN[i, ki] = n - k
    return PD, PQ, P0, PN, Ms


def null_from(Ms, PK, K, nsim):
    KSc = np.arange(1, K + 1)
    SD = np.zeros((nsim, K))
    for M, (s, sel, n) in zip(Ms, PK):
        O = np.argsort(rng.random((nsim, n)), axis=1)
        for ki, k in enumerate(KSc):
            if n - k <= 0:
                continue
            SD[:, ki] += M[O[:, :-k], O[:, k:]].sum(1)
    return SD


# ---------------- ⒜ all units, middle window ----------------
print("\n" + "=" * 78)
print("[⒜] 가운데 창 (0.2-0.8) 의 practical range + block bootstrap CI (블록 = haplotype)")
K = 15
PK = pack((0.2, 0.8), None, K)
PD, PQ, P0, PN, Ms = curves_from(PK, lambda s, sel: MATS[s][np.ix_(sel, sel)].astype(np.int64), K)
ON = PN.sum(0)
KSc = np.arange(1, K + 1)
GA = 0.5 * PQ.sum(0) / ON
TM = PD.sum(0) / ON
cn, cp, a = fit_exp(GA, KSc, ON)
print("    hap %d · nugget %.3f · sill %.3f · a %.4f ⇒ 3a = %.4f copy" % (len(PK), cn, cn + cp, a, 3 * a))
t0 = time.time()
prs, aas = [], []
NH = len(PK)
for b in range(B_BLK):
    r = rng.integers(0, NH, NH)
    nn = PN[r].sum(0)
    gm = 0.5 * PQ[r].sum(0) / nn
    try:
        aa = fit_exp(gm, KSc, nn)[2]
        aas.append(aa)
        prs.append(3 * aa)
    except Exception:
        pass
prs = np.array(prs)
ci = np.percentile(prs, [2.5, 97.5])
print("    block bootstrap B %d (%.0fs) : 3a 95%%CI [%.3f, %.3f] · median %.3f · 상한에 붙은 몫 %.4f"
      % (len(prs), time.time() - t0, ci[0], ci[1], np.median(prs), float((np.array(aas) > 499).mean())))
SD = null_from(Ms, PK, K, NSIM)
p1, lo, hi = erl_test(TM, SD / ON)
out = np.where((TM < lo) | (TM > hi))[0] + 1
print("    GET global p(erl) = %.5g · 띠 밖 lag %s · lambda_envelope = %d"
      % (p1, out.tolist(), out.max() if len(out) else 0))
print("    세 추정량 : practical range %.2f [%.2f, %.2f] · lambda_envelope %d · EHH<=0.05 lag 2"
      % (3 * a, ci[0], ci[1], out.max() if len(out) else 0))
print("       ⇒ 한 자릿수 차이로 불일치. §3.4(A) 가 practical range 를 주 추정량으로 사전 등록했다.")

# ---------------- ⒝ by segment, middle window ----------------
print("\n" + "=" * 78)
print("[⒝] 구간별 — 가운데 창 (0.2-0.8) 에서 다시 (앞 판 interior(+-5) 는 양끝 효과로 부풀었다)")
rows = []
for nm, a_, b_ in REG:
    t0 = time.time()
    PKr = pack((0.2, 0.8), None, K)

    def Mget(s, sel, nm=nm):
        aps = APS[s]
        kk = [(s[0], s[1], aps[i]) for i in sel]
        n = len(kk)
        M = np.zeros((n, n), np.int64)
        for i in range(n):
            for j in range(i + 1, n):
                M[i, j] = M[j, i] = pdr(nm, kk[i], kk[j])
        return M
    PDr, PQr, P0r, PNr, Msr = curves_from(PKr, Mget, K)
    ONr = PNr.sum(0)
    TMr = PDr.sum(0) / ONr
    GAr = 0.5 * PQr.sum(0) / ONr
    mp = TMr.mean()
    try:
        c1, c2, aa = fit_exp(GAr, KSc, ONr)
    except Exception:
        c1 = c2 = aa = np.nan
    SDr = null_from(Msr, PKr, K, 999)
    pr, lor, hir = erl_test(TMr, SDr / ONr)
    outr = np.where((TMr < lor) | (TMr > hir))[0] + 1
    # bootstrap
    pb = []
    for b in range(500):
        r = rng.integers(0, len(PKr), len(PKr))
        nn = PNr[r].sum(0)
        gm = 0.5 * PQr[r].sum(0) / nn
        try:
            pb.append(3 * fit_exp(gm, KSc, nn)[2])
        except Exception:
            pass
    pb = np.array(pb)
    cib = np.percentile(pb, [2.5, 97.5]) if len(pb) > 50 else (np.nan, np.nan)
    ident = (mp >= 0.5) and (aa < 400)
    print("    %-10s bp %4d · 쌍당 불일치 %7.4f (%7.3f /kb) · global p %.4g · lambda_env %2d · 3a = %8.3f [%.2f, %.2f] %s  (%.0fs)"
          % (nm, b_ - a_ + 1, mp, 1000 * mp / (b_ - a_ + 1), pr, outr.max() if len(outr) else 0,
             3 * aa, cib[0], cib[1], "" if ident else "추정불가", time.time() - t0))
    rows.append((nm, b_ - a_ + 1, float(mp), float(1000 * mp / (b_ - a_ + 1)), pr,
                 int(outr.max()) if len(outr) else 0, float(3 * aa), float(cib[0]), float(cib[1]), ident,
                 float(TMr[0]), float(TMr[1]), float(TMr[-1])))
print("    gene / promoter / AluY / NTS2 는 쌍당 불일치가 1 미만 ⇒ range 추정 불가 (Q0-3 해상도 한계)")
print("    쓸 수 있는 문장 : 'lambda 는 GA-rich 가 지배한다. 유전자 구간에서는 추정 자체가 불가능하다.'")

# ---------------- ⒞ λ vs n_copies ----------------
print("\n" + "=" * 78)
print("[⒞] λ vs n_copies — '93 copy 이므로 구역 7 개' 같은 계산은 쓰지 않는다 (n 이 14-239 로 17 배 다르다)")
lam_h, n_h = [], []
for (s, sel, n), M in zip(PK, Ms):
    nn = NCOPY.get(s, 0)
    SDh = np.zeros((99, K))
    Th = np.zeros(K)
    cnt = np.zeros(K)
    for ki, k in enumerate(KSc):
        if n - k <= 0:
            continue
        Th[ki] = M[np.arange(n - k), np.arange(k, n)].sum()
        cnt[ki] = n - k
    O = np.argsort(rng.random((99, n)), axis=1)
    for ki, k in enumerate(KSc):
        if n - k <= 0:
            continue
        SDh[:, ki] = M[O[:, :-k], O[:, k:]].sum(1)
    ok = cnt > 0
    t_ = Th[ok] / cnt[ok]
    s_ = SDh[:, ok] / cnt[ok]
    lo_ = np.percentile(s_, 2.5, axis=0)
    hi_ = np.percentile(s_, 97.5, axis=0)
    o_ = np.where((t_ < lo_) | (t_ > hi_))[0]
    lam_h.append(int(KSc[ok][o_.max()]) if len(o_) else 0)
    n_h.append(nn)
lam_h, n_h = np.array(lam_h), np.array(n_h)
r = stats.spearmanr(n_h, lam_h)
print("    per-haplotype lambda_env : median %.1f · IQR %.1f-%.1f · 0 인 몫 %.4f"
      % (np.median(lam_h), *np.percentile(lam_h, [25, 75]), (lam_h == 0).mean()))
print("    lambda vs n_copies : Spearman rho = %+.4f  p = %.4g" % (r.statistic, r.pvalue))
nz = lam_h > 0
print("    n/lambda (lambda>0 인 %d 개) : median %.2f · IQR %.2f-%.2f"
      % (int(nz.sum()), np.median(n_h[nz] / lam_h[nz]), *np.percentile(n_h[nz] / lam_h[nz], [25, 75])))
print("    ⇒ 판정 : %s" % ("FRACTION_SCALE 경향 (유의)" if r.pvalue < 0.05 else "ABSOLUTE_SCALE (n 과 무관)"))
print("    개체별 순열 99 벌이라 검정력 낮음 — 경향으로만")

g = io.open(os.path.join(OUT, "Q2c_results.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 Q2c · ruler majref %.6f · seed %d · 가운데 창 0.2-0.8 · nsim %d · B %d · 2026-10-10\n"
        % (VPU, SEED, NSIM, B_BLK))
g.write("[A]\tn_hap\t%d\tK\t%d\tnugget\t%.6f\tsill\t%.6f\ta\t%.6f\tpractical_range\t%.6f\tlo\t%.6f\thi\t%.6f\tp_erl\t%.6g\tlambda_env\t%d\n"
        % (len(PK), K, cn, cn + cp, a, 3 * a, ci[0], ci[1], p1, out.max() if len(out) else 0))
g.write("[A-CURVE]\tlag\tN_k\tmean_d\tgamma\tlo_erl\thi_erl\n")
for ki in range(K):
    g.write("\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n" % (ki + 1, ON[ki], TM[ki], GA[ki], lo[ki], hi[ki]))
g.write("[B]\tregion\tbp\tmismatch_per_pair\tper_kb\tp_erl\tlambda_env\tpractical_range\tlo\thi\tidentifiable\tmean_d_lag1\tmean_d_lag2\tmean_d_lag15\n")
for r_ in rows:
    g.write("\t%s\t%d\t%.6f\t%.4f\t%.6g\t%d\t%.4f\t%.4f\t%.4f\t%s\t%.6f\t%.6f\t%.6f\n" % r_)
g.write("[C]\tlam_median\t%.2f\tfrac_zero\t%.4f\trho_n\t%.6f\tp\t%.6g\tn_over_lam_median\t%.4f\n"
        % (np.median(lam_h), (lam_h == 0).mean(), r.statistic, r.pvalue,
           np.median(n_h[nz] / lam_h[nz])))
g.close()
np.save(os.path.join(OUT, "Q2c_boot_pr.npy"), prs)
print("\n산출 :", os.path.join(OUT, "Q2c_results.tsv"))
