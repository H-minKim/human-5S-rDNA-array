# -*- coding: utf-8 -*-
r"""P08 §2 - Q1 : how different and how alike are the units of two haplotypes
   ruler = number of profile column mismatches (bp) fixed in Q0-1, identical = fully identical sequence
   three layers T1(within array) / T2(homolog) / T3(unrelated) x four quantities ①median ②identity ③identical partner ④f/φ
   CI is a haplotype-level cluster bootstrap (Field & Welsh 2007). The pair-level CI of P04 was wrong
   Ohta C2 is tested by individual label permutation
   A9 : size distribution of fully identical classes (same form as the EHH formula)"""
import os, sys, io, time, collections, pickle, itertools
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
os.makedirs(OUT, exist_ok=True)
SEED, EDGE, B_BOOT, N_PERM = 20261010, 5, 2000, 999
rng = np.random.default_rng(SEED)
import v11_io as G
from scipy import stats

CAN = G.CAN
LU = len(CAN)
U, _, _ = G.load("h1")
VPU = sum(len(t[4]) for v in U.values() for t in v) / sum(len(v) for v in U.values())
print("§0 자 : var_per_unit = %.6f" % VPU)
assert abs(VPU - 6.866) < 0.01
NCOPY = {("%s" % k[0], "%s" % k[1]): len(v) for k, v in U.items()}
POP = {}
PEDF = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/v11_readvar/05_reference/1kGP_3202_ped_population_2026-09-23.txt")
for Ln in io.open(PEDF, encoding="utf8"):
    f = Ln.split()
    if len(f) >= 7 and f[1] != "SampleID":
        POP[f[1]] = f[6]
print("집단 라벨 적재 %d 개체 (super-pop %s)" % (len(POP), sorted(set(POP.values()))[:8]))
with open(os.path.join(HERE, "OUT_P04", "Q0_alignment.pkl"), "rb") as f:
    DIFF = pickle.load(f)["diff"]
with open(os.path.join(HERE, "OUT_P04", "S3_matrices.pkl"), "rb") as f:
    S3 = pickle.load(f)
BAD, MATS, APS = S3["bad"], S3["mats"], S3["aps"]
keys = [k for k in DIFF if k not in BAD]

# ---------------- representation: (p,allele) set, position set, insertion dict, signature ----------------
SET, PSET, INS, SIG = {}, {}, {}, {}
for k in keys:
    sub, dele, ins, ul = DIFF[k]
    A = dict(sub)
    for p in dele:
        A[p] = "-"
    SET[k] = frozenset(A.items())
    PSET[k] = frozenset(A)
    ia = {}
    for p, s in ins:
        ia[p] = ia.get(p, "") + s
    INS[k] = ia
    SIG[k] = (";".join(sorted("%d%s" % (p, b) for p, b in A.items())) + "|"
              + ";".join("%d:%s" % (p, ia[p]) for p in sorted(ia)))


def pd2(k1, k2):
    """Number of profile column mismatches (Q0-1 ruler), set operation version"""
    S1, S2 = SET[k1], SET[k2]
    both = len(S1 & S2)
    d = len(S1) + len(S2) - 2 * both - (len(PSET[k1] & PSET[k2]) - both)
    I1, I2 = INS[k1], INS[k2]
    if I1 or I2:
        for p in set(I1) | set(I2):
            a, b = I1.get(p, ""), I2.get(p, "")
            if a != b:
                L = max(len(a), len(b))
                d += sum(1 for x, y in zip(a.ljust(L, "\0"), b.ljust(L, "\0")) if x != y)
    return d


# check: does it match the S3 matrix
chk = 0
for s in itertools.islice(sorted(MATS), 5):
    aps = APS[s]
    M = MATS[s]
    for i in range(min(len(aps), 12)):
        for j in range(i + 1, min(len(aps), 12)):
            assert pd2((s[0], s[1], aps[i]), (s[0], s[1], aps[j])) == M[i, j], (s, i, j)
            chk += 1
print("검산 pd2 == S3 행렬 : %d / %d 쌍 일치" % (chk, chk))

bySlot = collections.defaultdict(list)
for k in keys:
    bySlot[(k[0], k[1])].append(k[2])
for s in bySlot:
    bySlot[s].sort()
INT = {}
for s, aps in bySlot.items():
    N = NCOPY.get(s, 0)
    v = [ap for ap in aps if EDGE < ap <= N - EDGE]
    if len(v) >= 3:
        INT[s] = v
HSL = sorted(INT)
IDX = {s: i for i, s in enumerate(HSL)}
bySamp = collections.defaultdict(list)
for s in HSL:
    bySamp[s[0]].append(s)
PAIRED = [tuple(v) for v in bySamp.values() if len(v) == 2]
print("interior 슬롯 %d · 유닛 %d · 두 하플로 다 있는 개체 %d"
      % (len(HSL), sum(len(v) for v in INT.values()), len(PAIRED)))

# ---------------- ④ f / φ : exhaustive, via signature inner products ----------------
sid = {}
rowsig = []
for s in HSL:
    c = collections.Counter(SIG[(s[0], s[1], ap)] for ap in INT[s])
    rowsig.append(c)
    for g in c:
        sid.setdefault(g, len(sid))
NS = len(sid)
print("서로 다른 서명 %d (유닛 %d)" % (NS, sum(len(v) for v in INT.values())))
X = np.zeros((len(HSL), NS), np.int32)
for i, c in enumerate(rowsig):
    for g, n in c.items():
        X[i, sid[g]] = n
SAME = X.astype(np.int64) @ X.T                     # number of identical pairs
nvec = X.sum(1).astype(np.int64)
diagsame = (X.astype(np.int64) * (X.astype(np.int64) - 1)).sum(1)   # identical pairs within the same haplotype x2


def f_phi(idxs=None, exclude_self_indiv=True):
    """idxs = list of haplotype indices to use (duplicates allowed). f = within the same haplotype, φ = between different individuals"""
    I = np.arange(len(HSL)) if idxs is None else np.asarray(idxs)
    nn = nvec[I]
    num_f = diagsame[I].sum()
    den_f = (nn * (nn - 1)).sum()
    S = SAME[np.ix_(I, I)]
    samp = np.array([HSL[i][0] for i in I])
    diff = samp[:, None] != samp[None, :]
    num_p = int(S[diff].sum())
    den_p = int((nn[:, None] * nn[None, :])[diff].sum())
    return (num_f / den_f if den_f else np.nan, num_f, den_f,
            num_p / den_p if den_p else np.nan, num_p, den_p)


f_, nf, df_, phiU, npU, dpU = f_phi()
f_ = nf / df_
# φ_homolog : only individuals with both haplotypes
nh = dh = 0
for a, b in PAIRED:
    i, j = IDX[a], IDX[b]
    nh += int(SAME[i, j])
    dh += int(nvec[i]) * int(nvec[j])
phiH = nh / dh
# by population
samp_pop = {s[0]: POP.get(s[0], "?") for s in HSL}
nw = dw = nb = db = 0
for i in range(len(HSL)):
    for j in range(i + 1, len(HSL)):
        if HSL[i][0] == HSL[j][0]:
            continue
        pi, pj = samp_pop[HSL[i][0]], samp_pop[HSL[j][0]]
        if pi == "?" or pj == "?":
            continue
        c = int(SAME[i, j])
        d = int(nvec[i]) * int(nvec[j])
        if pi == pj:
            nw += c
            dw += d
        else:
            nb += c
            db += d
print("\n" + "=" * 78)
print("[Q1 ④] 쌍 동일 확률 (전수 · 서명 내적)")
print("  f (T1 array 내 interior)   = %.6f   (%d / %d)" % (f_, nf, df_))
print("  φ_homolog (T2)             = %.6f   (%d / %d · 개체 %d)" % (phiH, nh, dh, len(PAIRED)))
print("  φ_unrelated (T3)           = %.6f   (%d / %d)" % (phiU, npU, dpU))
print("    φ_within_pop             = %.6f   (%d / %d)" % (nw / dw, nw, dw))
print("    φ_between_pop            = %.6f   (%d / %d)" % (nb / db, nb, db))
print("  f/φ : homolog %.4f · unrelated %.4f · within_pop %.4f · between_pop %.4f"
      % (f_ / phiH, f_ / phiU, f_ / (nw / dw), f_ / (nb / db)))

# ---------------- cluster bootstrap (haplotype / individual unit) ----------------
print("\n[Q1 CI] cluster bootstrap — 재표본 단위 = haplotype(개체) · B = %d" % B_BOOT)
pi_ = np.array([IDX[a] for a, b in PAIRED])
pj_ = np.array([IDX[b] for a, b in PAIRED])
sh_num = np.array([SAME[i, j] for i, j in zip(pi_, pj_)], np.int64)
sh_den = nvec[pi_] * nvec[pj_]
bs_h, bs_f, bs_u = [], [], []
nP = len(PAIRED)
nH = len(HSL)
for b in range(B_BOOT):
    r = rng.integers(0, nP, nP)
    bs_h.append(sh_num[r].sum() / sh_den[r].sum())
    r2 = rng.integers(0, nH, nH)
    nn = nvec[r2]
    bs_f.append(diagsame[r2].sum() / (nn * (nn - 1)).sum())
    S = SAME[np.ix_(r2, r2)]
    samp = np.array([HSL[i][0] for i in r2])
    diff = samp[:, None] != samp[None, :]
    bs_u.append(S[diff].sum() / (nn[:, None] * nn[None, :])[diff].sum())
bs_h, bs_f, bs_u = map(np.array, (bs_h, bs_f, bs_u))
ciH = np.percentile(bs_h, [2.5, 97.5])
ciF = np.percentile(bs_f, [2.5, 97.5])
ciU = np.percentile(bs_u, [2.5, 97.5])
# the (wrong) binomial CI of P04
from statsmodels.stats.proportion import proportion_confint
p04H = proportion_confint(nh, dh, 0.05, "wilson")
p04U = proportion_confint(npU, dpU, 0.05, "wilson")
print("  f           %.6f  cluster 95%%CI [%.6f, %.6f]  폭 %.2e" % (f_, ciF[0], ciF[1], ciF[1] - ciF[0]))
print("  φ_homolog   %.6f  cluster 95%%CI [%.6f, %.6f]  폭 %.2e" % (phiH, ciH[0], ciH[1], ciH[1] - ciH[0]))
print("              P04 의 쌍 단위 Wilson CI [%.6f, %.6f]  폭 %.2e  ⇒ 과소평가 %.2f 배"
      % (p04H[0], p04H[1], p04H[1] - p04H[0], (ciH[1] - ciH[0]) / (p04H[1] - p04H[0])))
print("  φ_unrelated %.6f  cluster 95%%CI [%.6f, %.6f]  폭 %.2e" % (phiU, ciU[0], ciU[1], ciU[1] - ciU[0]))
print("              P04 의 쌍 단위 Wilson CI [%.6f, %.6f]  폭 %.2e  ⇒ 과소평가 %.2f 배"
      % (p04U[0], p04U[1], p04U[1] - p04U[0], (ciU[1] - ciU[0]) / (p04U[1] - p04U[0])))
# bootstrap of the difference (conservative, because resampling is independent rather than paired)
d_bs = bs_h - bs_u
ciD = np.percentile(d_bs, [2.5, 97.5])
print("  φ_homolog − φ_unrelated = %+.6f  cluster 95%%CI [%+.6f, %+.6f]  (0 포함 %s)"
      % (phiH - phiU, ciD[0], ciD[1], "YES ⇒ 검출 못 함" if ciD[0] <= 0 <= ciD[1] else "NO ⇒ 검출"))
print("  비 φ_homolog/φ_unrelated = %.4f  95%%CI [%.4f, %.4f]"
      % (phiH / phiU, *np.percentile(bs_h / bs_u, [2.5, 97.5])))

# ---------------- individual label permutation ----------------
print("\n[Q1 PERM] 개인 라벨 순열 — '같은 개인의 두 하플로' 관계를 무작위화 (N = %d)" % N_PERM)
flat = np.concatenate([pi_, pj_])
null = []
for b in range(N_PERM):
    p = rng.permutation(flat)
    a_, b_ = p[:nP], p[nP:]
    ok = np.array([HSL[x][0] != HSL[y][0] for x, y in zip(a_, b_)])
    if ok.sum() < nP * 0.8:
        continue
    num = np.array([SAME[x, y] for x, y in zip(a_[ok], b_[ok])], np.int64).sum()
    den = (nvec[a_[ok]] * nvec[b_[ok]]).sum()
    null.append(num / den)
null = np.array(null)
pperm = (1 + int((null >= phiH).sum())) / (len(null) + 1)
print("  관측 φ_homolog %.6f ↔ 귀무 mean %.6f  SD %.2e  p95 %.6f  p = %.5g"
      % (phiH, null.mean(), null.std(ddof=1), np.percentile(null, 95), pperm))
VERD_C2 = "검출 (생존)" if (pperm < 0.05 and ciD[0] > 0) else ("소멸" if pperm >= 0.05 else "부분")
print("  VERDICT Ohta C2 : %s" % VERD_C2)

# ---------------- ① median/IQR, ③ identical-partner carrier rate ----------------
print("\n[Q1 ①③] 쌍별 거리 분포 · 동일짝 보유율 — 크기 맞춘 비교 집합 (하플로 하나)")
t0 = time.time()
# T1 : exhaustive from the S3 matrix (interior)
d_T1 = []
tw0_T1 = tw5_T1 = n_T1 = 0
for s in HSL:
    aps = APS[s]
    pos = {ap: i for i, ap in enumerate(aps)}
    ii = np.array([pos[ap] for ap in INT[s]])
    M = MATS[s][np.ix_(ii, ii)]
    n = len(ii)
    tri = M[np.triu_indices(n, 1)]
    d_T1.append(tri)
    np.fill_diagonal(M, 9999)
    tw0_T1 += int((M.min(1) == 0).sum())
    tw5_T1 += int((M.min(1) <= 5).sum())
    n_T1 += n
d_T1 = np.concatenate(d_T1)
# T2 : homolog, exhaustive
d_T2 = []
tw0_T2 = tw5_T2 = n_T2 = 0
for a, b in PAIRED:
    ka = [(a[0], a[1], ap) for ap in INT[a]]
    kb = [(b[0], b[1], ap) for ap in INT[b]]
    M = np.empty((len(ka), len(kb)), np.int16)
    for i, k1 in enumerate(ka):
        for j, k2 in enumerate(kb):
            M[i, j] = pd2(k1, k2)
    d_T2.append(M.ravel())
    tw0_T2 += int((M.min(1) == 0).sum()) + int((M.min(0) == 0).sum())
    tw5_T2 += int((M.min(1) <= 5).sum()) + int((M.min(0) <= 5).sum())
    n_T2 += len(ka) + len(kb)
d_T2 = np.concatenate(d_T2)
print("   T2 전수 %d 쌍 (%.0fs)" % (len(d_T2), time.time() - t0))
# T3 : unrelated, size matched (one haplotype), 5 repeats
t0 = time.time()
d_T3 = []
tw0_T3 = tw5_T3 = n_T3 = 0
REPL = 5
for rep in range(REPL):
    for a in HSL:
        while True:
            b = HSL[rng.integers(0, len(HSL))]
            if b[0] != a[0]:
                break
        ka = [(a[0], a[1], ap) for ap in INT[a]]
        kb = [(b[0], b[1], ap) for ap in INT[b]]
        M = np.empty((len(ka), len(kb)), np.int16)
        for i, k1 in enumerate(ka):
            for j, k2 in enumerate(kb):
                M[i, j] = pd2(k1, k2)
        if rep == 0:
            d_T3.append(M.ravel())
        tw0_T3 += int((M.min(1) == 0).sum())
        tw5_T3 += int((M.min(1) <= 5).sum())
        n_T3 += len(ka)
d_T3 = np.concatenate(d_T3)
print("   T3 크기맞춤 %d 반복 · 유닛 비교 %d (%.0fs)" % (REPL, n_T3, time.time() - t0))


def wils(k, n):
    return proportion_confint(k, n, 0.05, "wilson")


print("\n  %-14s %10s %8s %8s %12s %12s" % ("층", "n_pair", "median", "IQR", "동일짝=0", "동일짝<=5"))
rows13 = []
for nm, d, k0, k5, nU in (("T1 array내", d_T1, tw0_T1, tw5_T1, n_T1),
                          ("T2 상동", d_T2, tw0_T2, tw5_T2, n_T2),
                          ("T3 무관(맞춤)", d_T3, tw0_T3, tw5_T3, n_T3)):
    q = np.percentile(d, [25, 75])
    c0, c5 = wils(k0, nU), wils(k5, nU)
    print("  %-14s %10d %8.1f %4.0f-%-4.0f %6.2f%% [%.2f,%.2f] %6.2f%% [%.2f,%.2f]"
          % (nm, len(d), np.median(d), q[0], q[1],
             100 * k0 / nU, 100 * c0[0], 100 * c0[1], 100 * k5 / nU, 100 * c5[0], 100 * c5[1]))
    rows13.append((nm, len(d), float(np.median(d)), float(q[0]), float(q[1]), nU,
                   k0, k0 / nU, k5, k5 / nU, float(d.mean())))
print("  (T3 는 하플로 하나와만 비교 — T2 와 크기 맞춤 · %d 반복 평균)" % REPL)
print("  Naish 2021 CEN180 기준 : 같은 염색체 중복 보유 57.1-69.0%% ↔ 다른 염색체 동일 copy 0.3%%")
print("     (배수 비교 금지 — 종·단위길이 178 vs 2240 bp·copy 수가 전부 다르다)")

# ---------------- A9 : fully identical class ----------------
print("\n" + "=" * 78)
print("[Q1 A9] 완전동일 class 크기 분포 (= array 내 unit 하플로타입 빈도 스펙트럼)")
PK0, EK, VM, MX, NC, NU_, nuniq = [], [], [], [], [], [], []
Kall = []
runs_obs, runs_null = [], []
for s in HSL:
    aps = INT[s]
    gs = [SIG[(s[0], s[1], ap)] for ap in aps]
    c = collections.Counter(gs)
    n = len(aps)
    sizes = np.array(sorted(c.values(), reverse=True))
    K = np.array([c[g] - 1 for g in gs])
    Kall.append(K)
    PK0.append((K == 0).mean())
    EK.append(K.mean())
    VM.append(K.var() / K.mean() if K.mean() > 0 else np.nan)
    MX.append(int(sizes.max()))
    NC.append(len(c))
    NU_.append(n)
    nuniq.append(int((sizes == 1).sum()))
Kall = np.concatenate(Kall)
PK0, EK, VM, MX, NC, NU_ = map(np.array, (PK0, EK, VM, MX, NC, NU_))
print("  haplotype %d · 평균 n_copies %.2f" % (len(HSL), NU_.mean()))
print("  P(K=0)        mean %.4f ± %.4f   median %.4f   (전체 유닛 기준 %.4f)"
      % (PK0.mean(), PK0.std(ddof=1), np.median(PK0), (Kall == 0).mean()))
print("  E[K]          mean %.4f ± %.4f   median %.4f   (전체 %.4f)"
      % (EK.mean(), EK.std(ddof=1), np.median(EK), Kall.mean()))
print("  분산/평균     mean %.4f   median %.4f   (전체 %.4f)"
      % (np.nanmean(VM), np.nanmedian(VM), Kall.var() / Kall.mean()))
print("  class 총수    mean %.2f ± %.2f   median %.0f   range %d-%d" % (NC.mean(), NC.std(ddof=1), np.median(NC), NC.min(), NC.max()))
print("  최대 class    mean %.2f   median %.0f   max %d" % (MX.mean(), np.median(MX), MX.max()))
print("  단일 class 몫 mean %.4f" % (np.array(nuniq) / NC).mean())
print("  K 히스토그램 (전체 유닛 %d) : " % len(Kall)
      + " ".join("%d:%d" % (v, c) for v, c in sorted(collections.Counter(Kall.tolist()).items())[:14]))
print("  역산 예측과 대조 : E[K] 2.30 ↔ 실측 %.4f · P(K=0) 0.4732 ↔ 실측 %.4f · class 53-63 종 ↔ 실측 %.1f"
      % (Kall.mean(), (Kall == 0).mean(), NC.mean()))
# independence consistency
pred = 1 - (1 - f_) ** (NU_.mean() - 1)
obs_ = (Kall > 0).mean()
print("  독립성 정합 : 예측 P(K>=1) = 1-(1-f)^(N-1) = %.4f ↔ 관측 %.4f ⇒ 배수 %.4f"
      % (pred, obs_, obs_ / pred))
# mean number of twins for units that carry a twin
kk = Kall[Kall > 0]
print("  쌍둥이 보유 unit 의 평균 쌍둥이 수 %.4f (역산 4.37)" % kk.mean())
# Poisson comparison
print("  Poisson 이면 P(K=0) = exp(-E[K]) = %.4f ↔ 관측 %.4f" % (np.exp(-Kall.mean()), (Kall == 0).mean()))
r_nb = Kall.mean() ** 2 / (Kall.var() - Kall.mean()) if Kall.var() > Kall.mean() else np.nan
print("  음이항 적률 : r = %.4f (역산 0.385) · 분산/평균 %.4f (역산 ~7)" % (r_nb, Kall.var() / Kall.mean()))

# ---------------- A9-4 : is a class contiguous along the array ----------------
print("\n  class 가 array 위에서 연속(run)인가 산발인가 — 귀무 = 같은 크기 분할을 무작위 배치 (N=%d)" % N_PERM)


def adjpairs(labels):
    return int(sum(1 for a, b in zip(labels[:-1], labels[1:]) if a == b))


obs_adj = nul_adj = 0
nulls = np.zeros(N_PERM)
for s in HSL:
    aps = INT[s]
    gs = [SIG[(s[0], s[1], ap)] for ap in aps]
    obs_adj += adjpairs(gs)
arr = [[SIG[(s[0], s[1], ap)] for ap in INT[s]] for s in HSL]
for b in range(N_PERM):
    t = 0
    for gs in arr:
        p = rng.permutation(len(gs))
        t += adjpairs([gs[i] for i in p])
    nulls[b] = t
p_adj = (1 + int((nulls >= obs_adj).sum())) / (N_PERM + 1)
print("     인접 동일 쌍 관측 %d ↔ 귀무 %.1f ± %.1f  배율 %.4f  p = %.5g"
      % (obs_adj, nulls.mean(), nulls.std(ddof=1), obs_adj / nulls.mean(), p_adj))
# distribution of position gaps within a class
gaps_o, gaps_n = [], []
for s in HSL:
    aps = np.array(INT[s])
    by = collections.defaultdict(list)
    for i, ap in enumerate(aps):
        by[SIG[(s[0], s[1], ap)]].append(ap)
    for g, v in by.items():
        if len(v) >= 2:
            gaps_o.extend(np.diff(sorted(v)).tolist())
    gs = [SIG[(s[0], s[1], ap)] for ap in aps]
    p = rng.permutation(len(aps))
    by2 = collections.defaultdict(list)
    for i, ip in enumerate(p):
        by2[gs[ip]].append(int(aps[i]))
    for g, v in by2.items():
        if len(v) >= 2:
            gaps_n.extend(np.diff(sorted(v)).tolist())
gaps_o, gaps_n = np.array(gaps_o), np.array(gaps_n)
mw = stats.mannwhitneyu(gaps_o, gaps_n, alternative="less")
print("     class 내 자리 간격 : 관측 median %.1f (n %d) ↔ 무작위 배치 median %.1f (n %d) · MWU p = %.4g"
      % (np.median(gaps_o), len(gaps_o), np.median(gaps_n), len(gaps_n), mw.pvalue))
print("  AluY 가짜 class 배제 : 서명 = unit 전체(2,240 자리 + 삽입) 의 완전 일치이므로 구조적으로 불가")

g = io.open(os.path.join(OUT, "Q1_results.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 Q1 · ruler majref %.6f · seed %d · B_boot %d · N_perm %d · 2026-10-10\n" % (VPU, SEED, B_BOOT, N_PERM))
g.write("# interior 슬롯 %d · 유닛 %d · 상동 개체 %d · 서명 %d\n"
        % (len(HSL), sum(len(v) for v in INT.values()), len(PAIRED), NS))
g.write("[Q1-4]\titem\tvalue\tnum\tden\n")
for nm, v, a, b in (("f_T1", f_, nf, df_), ("phi_homolog", phiH, nh, dh), ("phi_unrelated", phiU, npU, dpU),
                    ("phi_within_pop", nw / dw, nw, dw), ("phi_between_pop", nb / db, nb, db)):
    g.write("\t%s\t%.8f\t%d\t%d\n" % (nm, v, a, b))
g.write("[Q1-CI]\titem\test\tlo\thi\twidth\tp04_lo\tp04_hi\tunderest_x\n")
g.write("\tf\t%.8f\t%.8f\t%.8f\t%.3e\t\t\t\n" % (f_, ciF[0], ciF[1], ciF[1] - ciF[0]))
g.write("\tphi_homolog\t%.8f\t%.8f\t%.8f\t%.3e\t%.8f\t%.8f\t%.3f\n"
        % (phiH, ciH[0], ciH[1], ciH[1] - ciH[0], p04H[0], p04H[1], (ciH[1] - ciH[0]) / (p04H[1] - p04H[0])))
g.write("\tphi_unrelated\t%.8f\t%.8f\t%.8f\t%.3e\t%.8f\t%.8f\t%.3f\n"
        % (phiU, ciU[0], ciU[1], ciU[1] - ciU[0], p04U[0], p04U[1], (ciU[1] - ciU[0]) / (p04U[1] - p04U[0])))
g.write("\tdiff_h_minus_u\t%.8f\t%.8f\t%.8f\t\t\t\t\n" % (phiH - phiU, ciD[0], ciD[1]))
g.write("\tperm_p_individual_label\t%.6g\tnull_mean\t%.8f\tnull_sd\t%.3e\tVERDICT\t%s\n"
        % (pperm, null.mean(), null.std(ddof=1), VERD_C2))
g.write("[Q1-13]\ttier\tn_pair\tmedian\tq25\tq75\tn_unit\ttwin0\tfrac0\ttwin5\tfrac5\tmean_d\n")
for r in rows13:
    g.write("\t%s\t%d\t%.1f\t%.0f\t%.0f\t%d\t%d\t%.6f\t%d\t%.6f\t%.4f\n" % r)
g.write("[A9]\titem\tmean\tsd\tmedian\tmin\tmax\n")
for nm, a in (("P_K0", PK0), ("E_K", EK), ("n_class", NC.astype(float)), ("max_class", MX.astype(float)), ("n_copies", NU_.astype(float))):
    g.write("\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" % (nm, a.mean(), a.std(ddof=1), np.median(a), a.min(), a.max()))
g.write("\tpooled_P_K0\t%.6f\tpooled_E_K\t%.6f\tvar_over_mean\t%.6f\tnb_r\t%.6f\n"
        % ((Kall == 0).mean(), Kall.mean(), Kall.var() / Kall.mean(), r_nb))
g.write("\tindep_pred_PK1\t%.6f\tobs_PK1\t%.6f\tratio\t%.6f\n" % (pred, obs_, obs_ / pred))
g.write("\tadj_same_obs\t%d\tnull_mean\t%.2f\tratio\t%.6f\tp\t%.6g\n" % (obs_adj, nulls.mean(), obs_adj / nulls.mean(), p_adj))
g.write("\tgap_median_obs\t%.2f\tgap_median_null\t%.2f\tMWU_p\t%.6g\n" % (np.median(gaps_o), np.median(gaps_n), mw.pvalue))
g.close()
np.save(os.path.join(OUT, "Q1_K_hist.npy"), Kall)
np.save(os.path.join(OUT, "Q1_boot.npy"), np.vstack([bs_f, bs_h, bs_u]))
np.save(os.path.join(OUT, "Q1_perm_null.npy"), null)
print("\n산출 :", os.path.join(OUT, "Q1_results.tsv"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('Q1_results.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
