# -*- coding: utf-8 -*-
r"""P08 Q2b, three additions to Q2
   (a) restrict to the middle window (normalized position 0.2-0.8) + free order permutation + GET
       Reason: with decile-stratified permutation each decile spans 8 copies, so the null mimics the short-lag signal itself
               (the stratifying variable is collinear with the predictor) => detrending must be done by window restriction
   (b) section 3.8 practical range / global p per segment (the 119 bp gene may lack enough diagnostic sites)
   (c) second output path = direct pairwise alignment with parasail (TRASH/ModDotPlot not installed)"""
import os, sys, io, re, time, collections, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
SEED, EDGE, NSIM, KMAX, DEC = 20261010, 5, 2499, 25, 10
rng = np.random.default_rng(SEED)
import v11_io as G
import parasail
from scipy import stats, optimize

CAN = G.CAN
LU = len(CAN)
CAN_NP = np.frombuffer(CAN.encode("ascii"), np.uint8)
REG = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA-rich", 1179, 1551),
       ("NTS2", 1552, 2031), ("promoter", 2032, 2121), ("gene", 2122, 2240)]
MAT = parasail.matrix_create("ACGT", 2, -3)
GO, GE = 10, 1
CIG = re.compile(r"(\d+)([=XIDM])")
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


def st_test(obs, sims, alpha=0.05):
    mu = sims.mean(0)
    sd = np.where(sims.std(0, ddof=1) > 0, sims.std(0, ddof=1), 1e-12)
    us = np.abs((sims - mu) / sd).max(1)
    uo = float(np.abs((np.asarray(obs) - mu) / sd).max())
    return (1 + int((us >= uo).sum())) / (len(us) + 1), uo, float(np.quantile(us, 1 - alpha))


def fit_exp(gam, ks, w, start=3):
    m = ks >= start
    x, y, ww = ks[m].astype(float), gam[m], w[m]

    def f(p):
        return np.sqrt(ww) * (p[0] + p[1] * (1 - np.exp(-x / max(p[2], 1e-6))) - y)
    r = optimize.least_squares(f, [y.min(), max(y.max() - y.min(), 1e-6), 5.0],
                               bounds=([0, 0, 0.1], [np.inf, np.inf, 500]))
    return r.x


def run_get(SLc, K, label, strat=False, nsim=NSIM):
    NH = len(SLc)
    KSc = np.arange(1, K + 1)
    PD = np.zeros((NH, K))
    PQ = np.zeros((NH, K))
    P0 = np.zeros((NH, K))
    PN = np.zeros((NH, K))
    for i, (M, dec, n) in enumerate(SLc):
        for ki, k in enumerate(KSc):
            if n - k <= 0:
                continue
            v = M[np.arange(n - k), np.arange(k, n)]
            PD[i, ki] = v.sum()
            PQ[i, ki] = (v * v).sum()
            P0[i, ki] = (v == 0).sum()
            PN[i, ki] = n - k
    ON = PN.sum(0)
    T = PD.sum(0) / ON
    GA = 0.5 * PQ.sum(0) / ON
    FK = P0.sum(0) / ON
    SD = np.zeros((nsim, K))
    SQ = np.zeros((nsim, K))
    for (M, dec, n) in SLc:
        key = (dec[None, :].astype(np.float64) * 1e6 + rng.random((nsim, n))) if strat else rng.random((nsim, n))
        O = np.argsort(key, axis=1)
        for ki, k in enumerate(KSc):
            if n - k <= 0:
                continue
            v = M[O[:, :-k], O[:, k:]]
            SD[:, ki] += v.sum(1)
            SQ[:, ki] += (v * v).sum(1)
    SD /= ON
    SQ = 0.5 * SQ / ON
    p1, lo, hi = erl_test(T, SD)
    p2, uo, uq = st_test(T, SD)
    out = np.where((T < lo) | (T > hi))[0] + 1
    pg, logg, higg = erl_test(GA, SQ)
    outg = np.where((GA < logg) | (GA > higg))[0] + 1
    try:
        cn, cp, a = fit_exp(GA, KSc, ON)
    except Exception:
        cn = cp = a = np.nan
    print("  ── %-22s hap %3d · lag 1-%d · mean d %.4f→%.4f" % (label, NH, K, T[0], T[-1]))
    print("     global p(erl, mean_d) = %.5g · p(st) = %.5g · 띠 밖 lag %s%s"
          % (p1, p2, out.tolist() if len(out) else "없음",
             ("  lambda_env = %d" % out.max()) if len(out) else ""))
    print("     gamma : p(erl) = %.5g · 띠 밖 %s · nugget %.3f sill %.3f a %.3f ⇒ 3a = %.3f"
          % (pg, outg.tolist() if len(outg) else "없음", cn, cn + cp, a, 3 * a))
    return dict(label=label, NH=NH, K=K, T=T, GA=GA, FK=FK, ON=ON,
                p_erl=p1, p_st=p2, out=out.tolist(), lam=int(out.max()) if len(out) else 0,
                p_erl_g=pg, out_g=outg.tolist(), nugget=cn, sill=cn + cp, a=a, pr=3 * a)


# ---------------- (a) middle window ----------------
print("\n" + "=" * 78)
print("[⒜] 추세 제거 — 가운데 창 제한 (정규화 위치 0.2-0.8) · 자유 순서 순열")
print("    (십분위 층화는 층화 변수가 예측변수와 공선 ⇒ 짧은 lag 신호를 귀무가 흉내낸다)")
RES = {}
for lo_, hi_, K in ((0.0, 1.0, KMAX), (0.2, 0.8, 15), (0.3, 0.7, 10)):
    SLc = []
    for s in sorted(MATS):
        aps = APS[s]
        N = NCOPY.get(s, 0)
        if lo_ == 0.0:
            sel = [i for i, ap in enumerate(aps) if EDGE < ap <= N - EDGE]
        else:
            sel = [i for i, ap in enumerate(aps) if lo_ * N < ap <= hi_ * N]
        if len(sel) < 3 * K // 2:
            continue
        ii = np.array(sel)
        M = MATS[s][np.ix_(ii, ii)].astype(np.int64)
        n = len(ii)
        SLc.append((M, np.minimum((np.arange(n) * DEC) // n, DEC - 1), n))
    tag = "interior(+-5)" if lo_ == 0.0 else ("central %.1f-%.1f" % (lo_, hi_))
    RES[tag] = run_get(SLc, K, tag, strat=False)

# ---------------- (b) per segment ----------------
print("\n" + "=" * 78)
print("[⒝] §3.8 구간별 — lambda = 11 은 unit 전체 거리이고 GA-rich 가 지배한다")
SETR, PSETR, INSR = {}, {}, {}
for nm, a_, b_ in REG:
    SETR[nm] = {}
    PSETR[nm] = {}
    INSR[nm] = {}
for k in keys:
    sub, dele, ins, ul = DIFF[k]
    A = dict(sub)
    for p in dele:
        A[p] = "-"
    ia = {}
    for p, s_ in ins:
        ia[min(max(p, 1), LU)] = ia.get(min(max(p, 1), LU), "") + s_
    for nm, a_, b_ in REG:
        d_ = {p: v for p, v in A.items() if a_ <= p <= b_}
        SETR[nm][k] = frozenset(d_.items())
        PSETR[nm][k] = frozenset(d_)
        INSR[nm][k] = {p: v for p, v in ia.items() if a_ <= p <= b_}
bySlot = collections.defaultdict(list)
for k in keys:
    bySlot[(k[0], k[1])].append(k[2])
for s in bySlot:
    bySlot[s].sort()


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


RESR = {}
for nm, a_, b_ in REG:
    t0 = time.time()
    SLc = []
    for s in sorted(bySlot):
        aps = bySlot[s]
        N = NCOPY.get(s, 0)
        sel = [ap for ap in aps if EDGE < ap <= N - EDGE]
        if len(sel) < 9:
            continue
        n = len(sel)
        kk = [(s[0], s[1], ap) for ap in sel]
        M = np.zeros((n, n), np.int64)
        for i in range(n):
            for j in range(i + 1, n):
                M[i, j] = M[j, i] = pdr(nm, kk[i], kk[j])
        SLc.append((M, np.minimum((np.arange(n) * DEC) // n, DEC - 1), n))
    mp = np.mean([Mx[np.triu_indices(nx, 1)].mean() for Mx, _, nx in SLc])
    print("\n  [%s] bp %d · 쌍당 평균 불일치 %.4f (%.3f /kb) · 행렬 %.0fs"
          % (nm, b_ - a_ + 1, mp, 1000 * mp / (b_ - a_ + 1), time.time() - t0))
    if mp < 0.5:
        print("     쌍당 불일치 %.3f < 0.5 ⇒ 진단 부위 부족으로 range 추정 불가 (Q0-3 해상도 한계)" % mp)
    RESR[nm] = run_get(SLc, 15, nm, strat=False, nsim=999)
    RESR[nm]["mismatch_per_pair"] = float(mp)
    RESR[nm]["per_kb"] = float(1000 * mp / (b_ - a_ + 1))

# ---------------- (c) second path: direct pairwise alignment ----------------
print("\n" + "=" * 78)
print("[⒞] 두 번째 산출 경로 — parasail 직접 쌍별 전역정렬 (consensus 경유 없음)")
print("    TRASH / ModDotPlot / StainedGlass 전부 미설치 ⇒ 도구 독립성은 정렬 방식 독립으로 대체")


def rebuild(k):
    sub, dele, ins, ul = DIFF[k]
    code = CAN_NP.copy()
    if sub:
        code[np.fromiter((p - 1 for p in sub), np.int64, len(sub))] = \
            np.fromiter((ord(b) for b in sub.values()), np.uint8, len(sub))
    keep = np.ones(LU, bool)
    if dele:
        keep[np.fromiter((p - 1 for p in dele), np.int64, len(dele))] = False
    kept = code[keep]
    if not ins:
        return kept.tobytes().decode("ascii")
    cs = np.cumsum(keep)
    ia = collections.defaultdict(str)
    for p, t in ins:
        ia[p] += t
    parts, prev = [], 0
    for p in sorted(ia):
        at = int(cs[p - 1]) if p >= 1 else 0
        parts.append(kept[prev:at].tobytes().decode("ascii"))
        parts.append(ia[p])
        prev = at
    parts.append(kept[prev:].tobytes().decode("ascii"))
    return "".join(parts)


def editd(s1, s2):
    r = parasail.nw_trace_striped_16(s1, s2, GO, GE, MAT)
    cg = r.cigar.decode
    if isinstance(cg, bytes):
        cg = cg.decode("ascii")
    X = ib = db = 0
    qi = ri = 0
    for n, op in CIG.findall(cg):
        n = int(n)
        if op in "=XM":
            for t in range(n):
                if s1[qi + t] != s2[ri + t]:
                    X += 1
            qi += n
            ri += n
        elif op == "I":
            ib += n
            qi += n
        elif op == "D":
            db += n
            ri += n
    return X + ib + db


INTS = {}
for s, aps in bySlot.items():
    N = NCOPY.get(s, 0)
    v = [ap for ap in aps if EDGE < ap <= N - EDGE]
    if len(v) >= 9:
        INTS[s] = v
SEQ = {}
t0 = time.time()
rows_c = []
NPL = {1: 3000, 2: 3000, 3: 3000, 4: 3000, 5: 3000}
for k in list(range(1, 11)) + [15, 20, 25]:
    pl = []
    for s, aps in INTS.items():
        for i in range(len(aps) - k):
            pl.append(((s[0], s[1], aps[i]), (s[0], s[1], aps[i + k])))
    rng.shuffle(pl)
    pl = pl[:NPL.get(k, 1500)]
    dv = []
    for k1, k2 in pl:
        for kx in (k1, k2):
            if kx not in SEQ:
                SEQ[kx] = rebuild(kx)
        dv.append(editd(SEQ[k1], SEQ[k2]))
    dv = np.array(dv)
    se = dv.std(ddof=1) / np.sqrt(len(dv))
    rows_c.append((k, len(dv), float(dv.mean()), float(se), float(np.median(dv))))
    print("    lag %-3d n %5d · 직접정렬 mean d = %.4f ± %.4f (SE) · median %.1f   [%.0fs]"
          % (k, len(dv), dv.mean(), se, np.median(dv), time.time() - t0), flush=True)
T0 = RES["interior(+-5)"]["T"]
kk = [r[0] for r in rows_c]
pa = np.array([T0[k - 1] for k in kk])
pb = np.array([r[2] for r in rows_c])
rho = stats.spearmanr(pa, pb)
print("\n    두 경로 곡선 일치 : Spearman rho = %+.4f (p %.3g) · Pearson r = %+.4f"
      % (rho.statistic, rho.pvalue, np.corrcoef(pa, pb)[0, 1]))
print("    lag2 최소 재현 : 프로파일 lag1 %.4f > lag2 %.4f ? %s · 직접정렬 lag1 %.4f > lag2 %.4f ? %s"
      % (T0[0], T0[1], T0[0] > T0[1], pb[0], pb[1], pb[0] > pb[1]))
dz = (pb[0] - pb[1]) / np.sqrt(rows_c[0][3] ** 2 + rows_c[1][3] ** 2)
print("    직접정렬의 lag1 − lag2 차 = %+.4f · z = %+.3f · 양측 p = %.4g"
      % (pb[0] - pb[1], dz, 2 * (1 - stats.norm.cdf(abs(dz)))))

g = io.open(os.path.join(OUT, "Q2b_results.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 Q2b · ruler majref %.6f · seed %d · nsim %d(구간 999) · 2026-10-10\n" % (VPU, SEED, NSIM))
g.write("[A]\twindow\tn_hap\tK\tp_erl_mean_d\tp_st\tlambda_env\tout_lags\tp_erl_gamma\tnugget\tsill\ta\tpractical_range\n")
for t_, v in RES.items():
    g.write("\t%s\t%d\t%d\t%.6g\t%.6g\t%d\t%s\t%.6g\t%.4f\t%.4f\t%.4f\t%.4f\n"
            % (t_, v["NH"], v["K"], v["p_erl"], v["p_st"], v["lam"], ",".join(map(str, v["out"])),
               v["p_erl_g"], v["nugget"], v["sill"], v["a"], v["pr"]))
g.write("[A-CURVE]\twindow\tlag\tN_k\tmean_d\tgamma\tf_k\n")
for t_, v in RES.items():
    for ki in range(v["K"]):
        g.write("\t%s\t%d\t%d\t%.6f\t%.6f\t%.8f\n" % (t_, ki + 1, v["ON"][ki], v["T"][ki], v["GA"][ki], v["FK"][ki]))
g.write("[B]\tregion\tbp\tmismatch_per_pair\tper_kb\tp_erl_mean_d\tlambda_env\tp_erl_gamma\tnugget\tsill\ta\tpractical_range\n")
for nm, a_, b_ in REG:
    v = RESR[nm]
    g.write("\t%s\t%d\t%.4f\t%.4f\t%.6g\t%d\t%.6g\t%.4f\t%.4f\t%.4f\t%.4f\n"
            % (nm, b_ - a_ + 1, v["mismatch_per_pair"], v["per_kb"], v["p_erl"], v["lam"],
               v["p_erl_g"], v["nugget"], v["sill"], v["a"], v["pr"]))
g.write("[B-CURVE]\tregion\tlag\tmean_d\tf_k\n")
for nm, a_, b_ in REG:
    v = RESR[nm]
    for ki in range(v["K"]):
        g.write("\t%s\t%d\t%.6f\t%.8f\n" % (nm, ki + 1, v["T"][ki], v["FK"][ki]))
g.write("[C]\tlag\tn_pair\tmean_d_direct\tSE\tmedian\tmean_d_profile\n")
for (k, n_, m_, se_, md_) in rows_c:
    g.write("\t%d\t%d\t%.6f\t%.6f\t%.1f\t%.6f\n" % (k, n_, m_, se_, md_, T0[k - 1]))
g.write("\tspearman_two_paths\t%.6f\tp\t%.6g\tpearson\t%.6f\n" % (rho.statistic, rho.pvalue, np.corrcoef(pa, pb)[0, 1]))
g.write("\tlag1_minus_lag2_direct\t%.6f\tz\t%.4f\tp\t%.6g\n" % (pb[0] - pb[1], dz, 2 * (1 - stats.norm.cdf(abs(dz)))))
g.close()
print("\n산출 :", os.path.join(OUT, "Q2b_results.tsv"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('Q2b_results.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
