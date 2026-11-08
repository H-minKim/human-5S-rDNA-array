# -*- coding: utf-8 -*-
r"""P08 X2 - within arrays only: what does the homogenization unit length depend on
   ⒜ array orientation gate (U84C) -> fix proximal/distal
   ⒝ longest shared identical stretch by array zone (proximal -> distal, 5 bins), haplotype cluster bootstrap
   ⒞ probability that a whole zone is identical, by unit zone (AluY/NTS1/GA-rich/NTS2/promoter/gene) x lag
      => fit a separate exponential decay length per zone (comparable because it is a 0-1 probability even when zone lengths differ)
   ⒟ two-component mixture per lag: is the distance dependence only in the long component
   ⒠ per-pair mismatch by zone + mean shared-identical length within zone (truncated at zone length)"""
import os, sys, io, time, collections, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
os.makedirs(OUT, exist_ok=True)
SEED, EDGE, B_BOOT = 20261010, 5, 1000
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
RC = {}
for f, Hx in G._rows(G.SUMF):
    RC[(f[Hx["sample"]], f[Hx["slot"]])] = f[Hx["revcomp"]]
with open(os.path.join(HERE, "OUT_P04", "Q0_alignment.pkl"), "rb") as f:
    DIFF = pickle.load(f)["diff"]
with open(os.path.join(HERE, "OUT_P04", "S3_matrices.pkl"), "rb") as f:
    BAD = pickle.load(f)["bad"]
keys = [k for k in DIFF if k not in BAD]
A_OF = {}
for k in keys:
    sub, dele, ins, ul = DIFF[k]
    d = dict(sub)
    for p in dele:
        d[p] = "-"
    A_OF[k] = d
bySlot = collections.defaultdict(list)
for k in keys:
    bySlot[(k[0], k[1])].append(k[2])
for s in bySlot:
    bySlot[s].sort()

# ---------------- ⒜ orientation gate ----------------
print("\n" + "=" * 78)
print("[⒜] array 방향 관문 — U84C (v11_io 기준 : 근위 0.194 · 원위 0.000)")
ORI = {}
for tag in ("Y", "N"):
    lo_d = lo_n = hi_d = hi_n = 0
    for s, aps in bySlot.items():
        if RC.get(s) != tag:
            continue
        N = NCOPY.get(s, 0)
        for ap in aps:
            has = (A_OF[(s[0], s[1], ap)].get(84) == "C")
            if ap <= 3:
                lo_n += 1
                lo_d += has
            if ap >= N - 2:
                hi_n += 1
                hi_d += has
    r_lo = lo_d / max(lo_n, 1)
    r_hi = hi_d / max(hi_n, 1)
    # the end closer to 0.194 is proximal
    ORI[tag] = "LOW_IS_PROX" if abs(r_lo - 0.194) < abs(r_hi - 0.194) else "HIGH_IS_PROX"
    print("    revcomp=%s : array_pos<=3 %.4f (n %d) · >=N-2 %.4f (n %d)  ⇒ %s"
          % (tag, r_lo, lo_n, r_hi, hi_n, "array_pos 1 = 근위" if ORI[tag] == "LOW_IS_PROX" else "array_pos 1 = 원위"))


def prox_frac(s, ap):
    """Position normalized from proximal 0 to distal 1"""
    N = NCOPY.get(s, 0)
    x = (ap - 0.5) / N
    return x if ORI.get(RC.get(s), "LOW_IS_PROX") == "LOW_IS_PROX" else 1.0 - x


INT = {}
for s, aps in bySlot.items():
    N = NCOPY.get(s, 0)
    v = [ap for ap in aps if EDGE < ap <= N - EDGE]
    if len(v) >= 9:
        INT[s] = v
ALL = {}
for s, aps in bySlot.items():
    if len(aps) >= 12:
        ALL[s] = aps
HSL = sorted(INT)
print("    interior haplotype %d · 전체(양끝 포함) haplotype %d" % (len(HSL), len(ALL)))


def mmpos(k1, k2):
    A1, A2 = A_OF[k1], A_OF[k2]
    o = [p for p, b in A1.items() if A2.get(p) != b] + [p for p in A2 if p not in A1]
    o.sort()
    return o


def runs_of(mm):
    if not mm:
        return [(LU, 0)]
    out = [(mm[0] - 1, 0)]
    for a, b in zip(mm[:-1], mm[1:]):
        out.append((b - a - 1, 1))
    out.append((LU - mm[-1], 0))
    return [(L, c) for L, c in out if L >= 0]


def cboot(vals, haps, B=B_BOOT, fn=np.mean):
    uh = np.unique(haps)
    idx = {h: np.where(haps == h)[0] for h in uh}
    out = np.empty(B)
    for b in range(B):
        r = uh[rng.integers(0, len(uh), len(uh))]
        sel = np.concatenate([idx[h] for h in r])
        out[b] = fn(vals[sel])
    return out


# ---------------- ⒝ by array zone ----------------
print("\n" + "=" * 78)
print("[⒝] 최장 공유동일 구간 by array 구역 (근위 0 → 원위 1 · 양끝 포함 · lag 1 쌍)")
NB = 5
bins = np.linspace(0, 1, NB + 1)
acc = [([], []) for _ in range(NB)]
t0 = time.time()
for s, aps in ALL.items():
    for a, b in zip(aps[:-1], aps[1:]):
        if b != a + 1:
            continue
        k1 = (s[0], s[1], a)
        k2 = (s[0], s[1], b)
        x = 0.5 * (prox_frac(s, a) + prox_frac(s, b))
        bi = min(int(x * NB), NB - 1)
        mm = mmpos(k1, k2)
        L = np.array([v for v, c in runs_of(mm)], float)
        acc[bi][0].append(L.max())
        acc[bi][1].append("%s|%s" % (s[0], s[1]))
rows_b = []
print("    %-14s %8s %12s %22s %10s" % ("구역", "쌍", "최장동일", "95%CI (cluster boot)", "불일치/쌍"))
for bi in range(NB):
    v = np.array(acc[bi][0])
    h = np.array(acc[bi][1])
    bs = cboot(v, h)
    ci = np.percentile(bs, [2.5, 97.5])
    lab = "%.1f-%.1f" % (bins[bi], bins[bi + 1])
    nm = ("근위단" if bi == 0 else ("원위단" if bi == NB - 1 else "중앙"))
    print("    %-14s %8d %12.3f %22s %10.3f"
          % (lab + " " + nm, len(v), v.mean(), "[%.1f, %.1f]" % (ci[0], ci[1]), np.nan))
    rows_b.append((lab, nm, len(v), float(v.mean()), ci[0], ci[1], float(np.median(v))))
print("    (%.0fs)" % (time.time() - t0))
mid = np.array([r[3] for r in rows_b])
print("    근위단 %.1f ↔ 중앙 %.1f ↔ 원위단 %.1f bp · 최대/최소 %.4f 배"
      % (mid[0], mid[NB // 2], mid[-1], mid.max() / mid.min()))
# both ends vs center test (haplotype unit)
e_, m_ = [], []
byh = collections.defaultdict(lambda: ([], []))
for bi in range(NB):
    for v, h in zip(acc[bi][0], acc[bi][1]):
        (byh[h][0] if bi in (0, NB - 1) else byh[h][1]).append(v)
for h, (a_, b_) in byh.items():
    if len(a_) >= 3 and len(b_) >= 3:
        e_.append(np.mean(a_))
        m_.append(np.mean(b_))
w = stats.wilcoxon(np.array(e_), np.array(m_))
print("    양끝 두 칸 vs 중앙 세 칸 (haplotype 단위 Wilcoxon · n=%d) : 끝 %.3f ↔ 중앙 %.3f · 차 %+.3f · p = %.4g"
      % (len(e_), np.mean(e_), np.mean(m_), np.mean(e_) - np.mean(m_), w.pvalue))

# ---------------- ⒞ identity probability by unit zone ----------------
print("\n" + "=" * 78)
print("[⒞] 구역 전체가 동일할 확률 by unit 구역 × lag  (0-1 확률이라 구역 길이가 달라도 비교 가능)")
LAGS = list(range(1, 26))
NP_ = 2000
PID = {nm: np.zeros(len(LAGS)) for nm, a, b in REG}
PIDn = np.zeros(len(LAGS))
MMR = {nm: np.zeros(len(LAGS)) for nm, a, b in REG}
t0 = time.time()
for li, lag in enumerate(LAGS):
    pl = []
    for s, aps in INT.items():
        for i in range(len(aps) - lag):
            pl.append(((s[0], s[1], aps[i]), (s[0], s[1], aps[i + lag])))
    rng.shuffle(pl)
    pl = pl[:NP_]
    cnt = {nm: 0 for nm, a, b in REG}
    mmc = {nm: 0 for nm, a, b in REG}
    for k1, k2 in pl:
        mm = mmpos(k1, k2)
        for nm, a_, b_ in REG:
            c = sum(1 for p in mm if a_ <= p <= b_)
            mmc[nm] += c
            if c == 0:
                cnt[nm] += 1
    for nm, a_, b_ in REG:
        PID[nm][li] = cnt[nm] / len(pl)
        MMR[nm][li] = mmc[nm] / len(pl)
    PIDn[li] = len(pl)
print("    lag 1-%d · 갈래당 쌍 %d (%.0fs)" % (max(LAGS), NP_, time.time() - t0))
print("    %-10s %5s %8s %8s %8s %8s %12s %14s" % ("구역", "bp", "P(1)", "P(2)", "P(10)", "P(25)", "불일치/쌍(1)", "감쇠길이 a"))
rows_c = []


def fit_dec(y, x):
    """P(k) = p_inf + (p1 - p_inf) exp(-(k-1)/a)"""
    def f(p):
        return p[0] + (y[0] - p[0]) * np.exp(-(x - 1) / max(p[1], 1e-6)) - y
    r = optimize.least_squares(f, [y[-1], 5.0], bounds=([0, 0.1], [1, 500]))
    return r.x


X = np.array(LAGS, float)
for nm, a_, b_ in REG:
    y = PID[nm]
    try:
        pinf, a_dec = fit_dec(y, X)
    except Exception:
        pinf, a_dec = np.nan, np.nan
    ok = (y[0] - y[-1]) > 0.005
    print("    %-10s %5d %8.4f %8.4f %8.4f %8.4f %12.4f %14s"
          % (nm, b_ - a_ + 1, y[0], y[1], y[9], y[24], MMR[nm][0],
             ("%.3f" % a_dec) if ok else "추정불가"))
    rows_c.append((nm, b_ - a_ + 1, y[0], y[1], y[9], y[24], float(MMR[nm][0]),
                   float(pinf), float(a_dec), ok))
print("    읽는 법 : P(1) 이 높을수록 그 구역이 이미 동일하다 · 감쇠길이 a 가 짧을수록 국소적으로만 동질")

# ---------------- ⒟ two-component mixture by lag ----------------
print("\n" + "=" * 78)
print("[⒟] 두 성분 혼합을 lag 별로 — 거리 의존성이 긴 성분에만 있는가")
SUBL = [1, 2, 3, 5, 8, 12, 18, 25]
rows_d = []
print("    %5s %8s %10s %10s %12s %12s" % ("lag", "n", "w(짧은)", "T1(짧은)", "T2(긴)", "평균 = wT1+(1-w)T2"))
for lag in SUBL:
    pl = []
    for s, aps in INT.items():
        for i in range(len(aps) - lag):
            pl.append(((s[0], s[1], aps[i]), (s[0], s[1], aps[i + lag])))
    rng.shuffle(pl)
    pl = pl[:2500]
    T, C = [], []
    for k1, k2 in pl:
        for L, c in runs_of(mmpos(k1, k2)):
            if L >= 2:
                T.append(L)
                C.append(c)
    t_ = np.array(T, float)
    c_ = np.array(C, int)

    def nll(p):
        w = 1 / (1 + np.exp(-p[0]))
        T1, T2 = max(p[1], 1e-3), max(p[2], 1e-3)
        pdf = w * np.exp(-t_ / T1) / T1 + (1 - w) * np.exp(-t_ / T2) / T2
        sf = w * np.exp(-t_ / T1) + (1 - w) * np.exp(-t_ / T2)
        return -(np.where(c_ == 1, np.log(np.clip(pdf, 1e-300, None)), np.log(np.clip(sf, 1e-300, None)))).sum()
    r = optimize.minimize(nll, [0.0, t_.mean() / 4, t_.mean() * 3], method="Nelder-Mead",
                          options=dict(maxiter=6000, xatol=1e-5, fatol=1e-5))
    w_ = 1 / (1 + np.exp(-r.x[0]))
    T1, T2 = sorted([r.x[1], r.x[2]])
    if r.x[1] > r.x[2]:
        w_ = 1 - w_
    mean_ = w_ * T1 + (1 - w_) * T2
    print("    %5d %8d %10.4f %10.3f %12.3f %12.3f" % (lag, len(t_), w_, T1, T2, mean_))
    rows_d.append((lag, len(t_), float(w_), float(T1), float(T2), float(mean_)))
t1v = np.array([r[3] for r in rows_d])
t2v = np.array([r[4] for r in rows_d])
wv = np.array([r[2] for r in rows_d])
for nm, v in (("T1 (짧은 성분)", t1v), ("T2 (긴 성분)", t2v), ("w (짧은 성분 비중)", wv)):
    rr = stats.spearmanr(SUBL, v)
    print("    lag ↔ %-18s : Spearman rho = %+.4f  p = %.4g" % (nm, rr.statistic, rr.pvalue))

# ---------------- ⒠ mean shared-identical length within zone ----------------
print("\n" + "=" * 78)
print("[⒠] 구역 안 평균 공유동일 길이 = 구역 bp / (구역 불일치 + 1) · 구역 길이로 절단")
print("    %-10s %5s %14s %16s %10s" % ("구역", "bp", "불일치/쌍(lag1)", "평균 동일 길이", "절단?"))
rows_e = []
for nm, a_, b_ in REG:
    bp = b_ - a_ + 1
    m1 = MMR[nm][0]
    L = bp / (m1 + 1)
    cens = L >= bp * 0.9
    print("    %-10s %5d %14.4f %16.2f %10s" % (nm, bp, m1, min(L, bp), "절단" if cens else ""))
    rows_e.append((nm, bp, float(m1), float(min(L, bp)), cens))
print("    절단 = 쌍당 불일치가 1 미만이라 구역 전체가 보통 동일하다 ⇒ 그 구역의 '동질화 길이' 는 구역보다 길다")

g = io.open(os.path.join(OUT, "X2_zone_region.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 X2 · ruler majref %.6f · seed %d · B %d · 2026-10-10\n" % (VPU, SEED, B_BOOT))
g.write("# 방향 : revcomp Y -> %s · N -> %s (U84C 실측)\n" % (ORI.get("Y"), ORI.get("N")))
g.write("[B]\tbin\tzone\tn_pair\tmaxrun_mean\tlo\thi\tmaxrun_median\n")
for r in rows_b:
    g.write("\t%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.1f\n" % r)
g.write("\tends_vs_mid\tn\t%d\tends\t%.4f\tmid\t%.4f\tp\t%.6g\n" % (len(e_), np.mean(e_), np.mean(m_), w.pvalue))
g.write("[C]\tregion\tbp\tP_lag1\tP_lag2\tP_lag10\tP_lag25\tmismatch_lag1\tp_inf\tdecay_a\tidentifiable\n")
for r in rows_c:
    g.write("\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\n" % r)
g.write("[C-CURVE]\tlag\t%s\n" % "\t".join(nm for nm, a, b in REG))
for li, lag in enumerate(LAGS):
    g.write("\t%d\t%s\n" % (lag, "\t".join("%.6f" % PID[nm][li] for nm, a, b in REG)))
g.write("[D]\tlag\tn\tw_short\tT1\tT2\tmean\n")
for r in rows_d:
    g.write("\t%d\t%d\t%.6f\t%.4f\t%.4f\t%.4f\n" % r)
g.write("[E]\tregion\tbp\tmismatch_lag1\tmean_identical_len\tcensored\n")
for r in rows_e:
    g.write("\t%s\t%d\t%.6f\t%.4f\t%s\n" % r)
g.close()
np.savez_compressed(os.path.join(OUT, "X2_curves.npz"),
                    LAGS=np.array(LAGS), **{("PID_" + nm): PID[nm] for nm, a, b in REG},
                    **{("MMR_" + nm): MMR[nm] for nm, a, b in REG})
print("\n산출 :", os.path.join(OUT, "X2_zone_region.tsv"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('X2_zone_region.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
