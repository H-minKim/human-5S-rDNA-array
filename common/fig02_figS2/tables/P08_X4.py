# -*- coding: utf-8 -*-
r"""P08 X4 - two things
   (a) Cause of the lag1 > lag2 phenomenon, testing the 'single interleaved odd unit' hypothesis
      Hypothesis: if some unit sits as a single interleaved unit, lag1 pairs often catch (normal, odd)
             while lag2 pairs skip over that unit and catch (normal, normal) => lag1 differs more
      Tests: (1) sign per haplotype  (2) recompute after removing single units  (3) after removing units carrying the proximal marker (2205C)
             (4) only haplotypes that lack such units  (5) is it a period-2 alternation (lags 3 and 4 included)
   (b) Homogenization patch length = Mansai 2011 interval-censored ML per lag (within arrays only)
      This is a different quantity from the longest shared identical stretch (descriptive) and the two-component mixture (distribution description)"""
import os, sys, io, time, collections, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
SEED, EDGE = 20261010, 5
rng = np.random.default_rng(SEED)
import v11_io as G
from scipy import stats, optimize

CAN = G.CAN
LU = len(CAN)
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
A_OF, SIG = {}, {}
for k in keys:
    sub, dele, ins, ul = DIFF[k]
    d = dict(sub)
    for p in dele:
        d[p] = "-"
    A_OF[k] = d
    ia = {}
    for p, s in ins:
        ia[p] = ia.get(p, "") + s
    SIG[k] = (";".join(sorted("%d%s" % (q, b) for q, b in d.items())) + "|"
              + ";".join("%d:%s" % (q, ia[q]) for q in sorted(ia)))

# ---------------- middle window groups ----------------
PACK = []
for s in sorted(MATS):
    aps = APS[s]
    N = NCOPY.get(s, 0)
    sel = [i for i, ap in enumerate(aps) if 0.2 * N < ap <= 0.8 * N]
    if len(sel) < 23:
        continue
    ii = np.array(sel)
    M = MATS[s][np.ix_(ii, ii)].astype(np.int64)
    kk = [(s[0], s[1], aps[i]) for i in sel]
    PACK.append((s, M, kk, len(sel)))
print("가운데 창 haplotype %d · 유닛 %d" % (len(PACK), sum(x[3] for x in PACK)))
KM = 6


def curve(mask_fn=None):
    """mask_fn(slot, kk, i) -> use the unit if True"""
    num = np.zeros(KM)
    den = np.zeros(KM)
    for s, M, kk, n in PACK:
        ok = np.ones(n, bool) if mask_fn is None else np.array([mask_fn(s, kk, i) for i in range(n)])
        for k in range(1, KM + 1):
            if n - k <= 0:
                continue
            i = np.arange(n - k)
            m = ok[i] & ok[i + k]
            if m.any():
                num[k - 1] += M[i[m], i[m] + k].sum()
                den[k - 1] += m.sum()
    return num / den, den


print("\n" + "=" * 78)
print("[⒜] lag1 > lag2 현상")
base, bn = curve()
print("    전체        : " + " · ".join("lag%d %.4f" % (k + 1, base[k]) for k in range(KM))
      + "   (lag1-lag2 = %+.4f)" % (base[0] - base[1]))

# (1) sign per haplotype
sgn = []
for s, M, kk, n in PACK:
    a1 = M[np.arange(n - 1), np.arange(1, n)].mean()
    a2 = M[np.arange(n - 2), np.arange(2, n)].mean()
    sgn.append(a1 - a2)
sgn = np.array(sgn)
bt = stats.binomtest(int((sgn > 0).sum()), len(sgn), 0.5)
print("    ① haplotype 별 lag1-lag2 : 양수 %d / %d = %.4f · median %+.4f · 부호검정 p = %.4g"
      % (int((sgn > 0).sum()), len(sgn), (sgn > 0).mean(), np.median(sgn), bt.pvalue))

# (2) single interleaved units
ODD = {}
nodd = 0
for s, M, kk, n in PACK:
    o = np.zeros(n, bool)
    for i in range(1, n - 1):
        if M[i - 1, i + 1] < min(M[i - 1, i], M[i, i + 1]):
            o[i] = True
    ODD[s] = o
    nodd += int(o.sum())
print("    ② 낱개 끼어든 유닛 = 양 이웃보다 '이웃끼리' 가 더 닮은 자리 : %d / %d = %.4f"
      % (nodd, sum(x[3] for x in PACK), nodd / sum(x[3] for x in PACK)))
c2, d2 = curve(lambda s, kk, i: not ODD[s][i])
print("       뺀 뒤   : " + " · ".join("lag%d %.4f" % (k + 1, c2[k]) for k in range(KM))
      + "   (lag1-lag2 = %+.4f)" % (c2[0] - c2[1]))
print("       ⇒ lag1-lag2 가 %+.4f → %+.4f  (%.1f%% 설명)"
      % (base[0] - base[1], c2[0] - c2[1],
         100 * (1 - (c2[0] - c2[1]) / (base[0] - base[1])) if (base[0] - base[1]) != 0 else np.nan))

# (3) units carrying the proximal marker 2205C
c3, d3 = curve(lambda s, kk, i: A_OF[kk[i]].get(2205) != "C")
n2205 = sum(1 for s, M, kk, n in PACK for i in range(n) if A_OF[kk[i]].get(2205) == "C")
print("    ③ 근위 표지 2205C 보유 유닛 %d (%.4f) 를 뺀 뒤 : lag1 %.4f · lag2 %.4f (차 %+.4f)"
      % (n2205, n2205 / sum(x[3] for x in PACK), c3[0], c3[1], c3[0] - c3[1]))

# (4) only haplotypes with no single units at all
keep = [s for s, M, kk, n in PACK if ODD[s].sum() == 0]
print("    ④ 낱개 유닛이 0 인 haplotype : %d 개" % len(keep))

# (5) is it a period-2 alternation, spacing distribution of single units
gaps = []
for s, M, kk, n in PACK:
    w = np.flatnonzero(ODD[s])
    if len(w) >= 2:
        gaps.extend(np.diff(w).tolist())
gaps = np.array(gaps)
if len(gaps):
    print("    ⑤ 낱개 유닛 사이 간격 : median %.1f · 간격 1 인 몫 %.4f · 간격 2 인 몫 %.4f"
          % (np.median(gaps), (gaps == 1).mean(), (gaps == 2).mean()))
    print("       (주기 2 교대라면 간격 2 가 압도해야 한다)")

# pairs that include a single unit vs pairs that do not
inn, out = [], []
for s, M, kk, n in PACK:
    o = ODD[s]
    i = np.arange(n - 1)
    v = M[i, i + 1]
    m = o[i] | o[i + 1]
    inn.extend(v[m].tolist())
    out.extend(v[~m].tolist())
print("    lag1 쌍 : 낱개 유닛 포함 %d 쌍 평균 %.4f ↔ 미포함 %d 쌍 평균 %.4f (차 %+.4f)"
      % (len(inn), np.mean(inn), len(out), np.mean(out), np.mean(inn) - np.mean(out)))

# ---------------- (b) Mansai patch length ----------------
print("\n" + "=" * 78)
print("[⒝] 동질화 패치 길이 = Mansai 2011 구간검열 ML (array 내부만 · Model 3 지수 · Model 2 gamma2)")
bySlot = collections.defaultdict(list)
for k in keys:
    bySlot[(k[0], k[1])].append(k[2])
for s in bySlot:
    bySlot[s].sort()
INT = {}
for s, aps in bySlot.items():
    N = NCOPY.get(s, 0)
    v = [ap for ap in aps if EDGE < ap <= N - EDGE]
    if len(v) >= 9:
        INT[s] = v


def events(k1, k2):
    A1, A2 = A_OF[k1], A_OF[k2]
    pos = sorted(set(A1) | set(A2))
    lab = [1 if (p in A1 and A2.get(p) == A1[p]) else 0 for p in pos]
    out = []
    i, Mn = 0, len(pos)
    while i < Mn:
        if lab[i] == 0:
            i += 1
            continue
        j = i
        while j + 1 < Mn and lab[j + 1] == 1:
            j += 1
        lo = pos[j] - pos[i]
        hi = (pos[j + 1] if j + 1 < Mn else LU) - (pos[i - 1] if i - 1 >= 0 else 1)
        out.append((lo, hi, 0 if (i - 1 >= 0 and j + 1 < Mn) else 1))
        i = j + 1
    return out


def S3f(t, T):
    return np.exp(-t / T)


def S2f(t, T):
    return (1 + 2 * t / T) * np.exp(-2 * t / T)


def nll(T, lo, hi, cen, S):
    T = max(float(T), 1e-3)
    pr = np.where(cen == 1, S(lo, T), S(lo, T) - S(hi, T))
    return -np.log(np.clip(pr, 1e-300, 1.0)).sum()


def fit(lo, hi, cen, S):
    r = optimize.minimize_scalar(lambda T: nll(T, lo, hi, cen, S), bounds=(1, 20000), method="bounded")
    base_ = r.fun

    def g(T):
        return nll(T, lo, hi, cen, S) - base_ - 1.920729
    try:
        a_ = optimize.brentq(g, 1, r.x)
    except Exception:
        a_ = np.nan
    try:
        b_ = optimize.brentq(g, r.x, 20000)
    except Exception:
        b_ = np.nan
    return r.x, a_, b_


LAGS = [1, 2, 3, 5, 8, 12, 18, 25]
NP_ = 4000
rows = []
print("    %5s %9s %9s %12s %22s %12s" % ("lag", "쌍", "사건", "T (지수)", "95%CI", "T (gamma2)"))
t0 = time.time()
for lag in LAGS:
    pl = []
    for s, aps in INT.items():
        for i in range(len(aps) - lag):
            pl.append(((s[0], s[1], aps[i]), (s[0], s[1], aps[i + lag])))
    rng.shuffle(pl)
    pl = pl[:NP_]
    ev = [e for k1, k2 in pl for e in events(k1, k2)]
    lo = np.array([e[0] for e in ev], float)
    hi = np.array([e[1] for e in ev], float)
    cen = np.array([e[2] for e in ev])
    T3, l3, h3 = fit(lo, hi, cen, S3f)
    T2, l2, h2 = fit(lo, hi, cen, S2f)
    print("    %5d %9d %9d %12.2f %22s %12.2f"
          % (lag, len(pl), len(ev), T3, "[%.1f, %.1f]" % (l3, h3), T2))
    rows.append((lag, len(pl), len(ev), float(T3), float(l3), float(h3), float(T2), float(l2), float(h2),
                 float(cen.mean())))
r_ = stats.spearmanr(LAGS, [x[3] for x in rows])
print("    lag ↔ 패치 길이 : Spearman rho = %+.4f  p = %.4g   (%.0fs)" % (r_.statistic, r_.pvalue, time.time() - t0))
print("    패치 길이 / unit 길이 : lag1 %.4f · lag25 %.4f" % (rows[0][3] / LU, rows[-1][3] / LU))

g = io.open(os.path.join(OUT, "X4_lag2_and_patch.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 X4 · ruler majref %.6f · seed %d · 2026-10-10\n" % (VPU, SEED))
g.write("[A]\tset\t%s\n" % "\t".join("lag%d" % (k + 1) for k in range(KM)))
g.write("\tall\t%s\n" % "\t".join("%.6f" % v for v in base))
g.write("\tdrop_odd\t%s\n" % "\t".join("%.6f" % v for v in c2))
g.write("\tdrop_2205C\t%s\n" % "\t".join("%.6f" % v for v in c3))
g.write("\tn_odd\t%d\tfrac\t%.6f\tn_2205C\t%d\n" % (nodd, nodd / sum(x[3] for x in PACK), n2205))
g.write("\thap_sign_pos\t%d\tn\t%d\tp\t%.6g\tmedian_diff\t%.6f\n"
        % (int((sgn > 0).sum()), len(sgn), bt.pvalue, np.median(sgn)))
g.write("\tlag1_with_odd\t%.6f\tn\t%d\tlag1_without\t%.6f\tn\t%d\n"
        % (np.mean(inn), len(inn), np.mean(out), len(out)))
if len(gaps):
    g.write("\toddgap_median\t%.2f\tfrac1\t%.6f\tfrac2\t%.6f\n" % (np.median(gaps), (gaps == 1).mean(), (gaps == 2).mean()))
g.write("[B]\tlag\tn_pair\tn_event\tT_exp\tlo\thi\tT_gamma2\tlo\thi\tcensored\n")
for r in rows:
    g.write("\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\n" % r)
g.write("\tspearman\t%.6f\tp\t%.6g\n" % (r_.statistic, r_.pvalue))
g.close()
print("\n산출 :", os.path.join(OUT, "X4_lag2_and_patch.tsv"))
