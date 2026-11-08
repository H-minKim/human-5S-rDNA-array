# -*- coding: utf-8 -*-
r"""P08 X5 - test lag1 > lag2 with a non-circular definition
   (2) of X4 was invalid: 'single unit' was defined as lag2 distance < lag1 distance, and then those units were removed
     => the result was built into the definition (the clue was the sign flipping by 435%, from +0.51 to -1.70)
   Corrected definition: whether a unit is odd is decided from its own content only (independent of neighbors and position)
     div_i = distance to the medoid unit of that haplotype, odd = top quartile within the haplotype
   Tests (a) are odd units scattered (number of adjacent odd-odd pairs vs a position-permutation null)
        (b) do lag1 pairs include odd units more often than lag2 pairs  <- decisive
        (c) lag1 - lag2 after removing odd units
        (d) does the same hold with the median instead of the mean (is it due to a few extreme pairs)"""
import os, sys, io, collections, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
SEED, N_PERM, KM = 20261010, 999, 15
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
PACK = []
for s in sorted(MATS):
    aps = APS[s]
    N = NCOPY.get(s, 0)
    sel = [i for i, ap in enumerate(aps) if 0.2 * N < ap <= 0.8 * N]
    if len(sel) < 23:
        continue
    ii = np.array(sel)
    M = MATS[s][np.ix_(ii, ii)].astype(np.int64)
    n = len(sel)
    med = int(np.argmin(M.sum(1)))           # medoid, independent of position
    div = M[med].astype(float)
    thr = np.percentile(div, 75)
    PACK.append((s, M, n, div, div > thr))
NU = sum(x[2] for x in PACK)
NOD = sum(int(x[4].sum()) for x in PACK)
print("가운데 창 haplotype %d · 유닛 %d · 특이(상위 사분위) %d (%.4f)" % (len(PACK), NU, NOD, NOD / NU))
print("    정의 : medoid 와의 거리 상위 사분위 — 이웃 거리를 전혀 쓰지 않는다")


def curve(use):
    num = np.zeros(KM)
    den = np.zeros(KM)
    med = [[] for _ in range(KM)]
    for (s, M, n, div, odd), ok in zip(PACK, use):
        for k in range(1, KM + 1):
            if n - k <= 0:
                continue
            i = np.arange(n - k)
            m = ok[i] & ok[i + k]
            if m.any():
                v = M[i[m], i[m] + k]
                num[k - 1] += v.sum()
                den[k - 1] += m.sum()
                med[k - 1].append(v)
    return num / den, den, [np.concatenate(x) if x else np.array([]) for x in med]


base, bn, bmed = curve([np.ones(x[2], bool) for x in PACK])
print("\n" + "=" * 78)
print("[기준] " + " · ".join("lag%d %.4f" % (k + 1, base[k]) for k in range(KM))
      + "   (lag1 − lag2 = %+.4f)" % (base[0] - base[1]))
print("       중앙값 : " + " · ".join("lag%d %.1f" % (k + 1, np.median(bmed[k])) for k in range(KM)))

# (a) are odd units scattered
print("\n[⒜] 특이 유닛의 배치 — 인접 '특이-특이' 쌍 수 ↔ haplotype 안 위치 순열 (N = %d)" % N_PERM)
obs = 0
for s, M, n, div, odd in PACK:
    obs += int((odd[:-1] & odd[1:]).sum())
nul = np.zeros(N_PERM)
for b in range(N_PERM):
    t = 0
    for s, M, n, div, odd in PACK:
        p = odd[rng.permutation(n)]
        t += int((p[:-1] & p[1:]).sum())
    nul[b] = t
p_iso = (1 + int((nul <= obs).sum())) / (N_PERM + 1)
print("    관측 %d ↔ 귀무 %.1f ± %.1f · 배율 %.4f · p(흩어짐) = %.4g"
      % (obs, nul.mean(), nul.std(ddof=1), obs / nul.mean(), p_iso))
print("    ⇒ %s" % ("흩어져 있지 않다 — 오히려 뭉쳐 있다" if obs > nul.mean() else "흩어져 있다"))

# (b) do lag1 pairs include odd units more often  <- decisive
print("\n[⒝] 결정적 — lag1 쌍과 lag2 쌍이 특이 유닛을 무는 빈도가 다른가")
print("    %5s %10s %10s %10s %12s" % ("lag", "특이 0 개", "1 개", "2 개", "1 개 몫"))
rows_b = []
for k in (1, 2, 3):
    c = [0, 0, 0]
    dsum = [0.0, 0.0, 0.0]
    for s, M, n, div, odd in PACK:
        if n - k <= 0:
            continue
        i = np.arange(n - k)
        t = odd[i].astype(int) + odd[i + k].astype(int)
        v = M[i, i + k]
        for z in (0, 1, 2):
            m = t == z
            c[z] += int(m.sum())
            dsum[z] += float(v[m].sum())
    tot = sum(c)
    print("    %5d %10d %10d %10d %12.5f" % (k, c[0], c[1], c[2], c[1] / tot))
    rows_b.append((k, c[0], c[1], c[2], c[1] / tot,
                   dsum[0] / max(c[0], 1), dsum[1] / max(c[1], 1), dsum[2] / max(c[2], 1)))
print("    칸별 평균 거리 : " + " | ".join(
    "lag%d  0개 %.3f · 1개 %.3f · 2개 %.3f" % (r[0], r[5], r[6], r[7]) for r in rows_b))
d1 = rows_b[0][4] - rows_b[1][4]
print("    '특이 1 개' 몫 : lag1 %.5f ↔ lag2 %.5f · 차 %+.5f" % (rows_b[0][4], rows_b[1][4], d1))
# size of lag1-lag2 explained by that difference in share
exp_gap = sum((rows_b[0][4 + 1 + z] * 0) for z in range(0))  # placeholder
m1 = np.array([rows_b[0][5], rows_b[0][6], rows_b[0][7]])
f1 = np.array([rows_b[0][1], rows_b[0][2], rows_b[0][3]], float)
f1 = f1 / f1.sum()
f2 = np.array([rows_b[1][1], rows_b[1][2], rows_b[1][3]], float)
f2 = f2 / f2.sum()
print("    lag1 의 칸별 평균을 lag2 의 구성비로 다시 섞으면 : %.4f (실제 lag1 %.4f · lag2 %.4f)"
      % (float((m1 * f2).sum()), base[0], base[1]))
print("       ⇒ 구성비 차이가 설명하는 몫 = %.1f%%"
      % (100 * (base[0] - float((m1 * f2).sum())) / (base[0] - base[1]) if base[0] != base[1] else np.nan))

# (c) odd units excluded
use = [~x[4] for x in PACK]
c2, d2, m2 = curve(use)
print("\n[⒞] 특이 유닛을 뺀 뒤 : " + " · ".join("lag%d %.4f" % (k + 1, c2[k]) for k in range(KM))
      + "   (lag1 − lag2 = %+.4f)" % (c2[0] - c2[1]))
print("    ⇒ lag1 − lag2 가 %+.4f → %+.4f" % (base[0] - base[1], c2[0] - c2[1]))

# (d) median and haplotype-level sign
print("\n[⒟] 평균이 아니라 중앙값으로 · 그리고 haplotype 단위 부호")
for nm, arr in (("전체", bmed), ("특이 제외", m2)):
    print("    %-8s 중앙값 lag1 %.1f · lag2 %.1f (차 %+.1f) · lag3 %.1f"
          % (nm, np.median(arr[0]), np.median(arr[1]), np.median(arr[0]) - np.median(arr[1]), np.median(arr[2])))
sg = []
for s, M, n, div, odd in PACK:
    sg.append(M[np.arange(n - 1), np.arange(1, n)].mean() - M[np.arange(n - 2), np.arange(2, n)].mean())
sg = np.array(sg)
bt = stats.binomtest(int((sg > 0).sum()), len(sg), 0.5)
w = stats.wilcoxon(sg)
print("    haplotype 별 lag1 − lag2 : 양수 %d / %d (%.4f) · median %+.4f · 부호검정 p = %.4g · Wilcoxon p = %.4g"
      % (int((sg > 0).sum()), len(sg), (sg > 0).mean(), np.median(sg), bt.pvalue, w.pvalue))
q = np.percentile(sg, [10, 25, 50, 75, 90])
print("    분포 십분위 : " + " · ".join("%.0f%% %+.3f" % (p, v) for p, v in zip([10, 25, 50, 75, 90], q)))

g = io.open(os.path.join(OUT, "X5_lag2_noncircular.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 X5 · 순환하지 않는 정의 (medoid 거리 상위 사분위) · seed %d · N_perm %d\n" % (SEED, N_PERM))
g.write("n_hap\t%d\tn_unit\t%d\tn_odd\t%d\tfrac\t%.6f\n" % (len(PACK), NU, NOD, NOD / NU))
g.write("[curve]\tset\t%s\n" % "\t".join("lag%d" % (k + 1) for k in range(KM)))
g.write("\tall\t%s\n" % "\t".join("%.6f" % v for v in base))
g.write("\tdrop_odd\t%s\n" % "\t".join("%.6f" % v for v in c2))
g.write("\tall_median\t%s\n" % "\t".join("%.1f" % np.median(bmed[k]) for k in range(KM)))
g.write("[iso]\tobs\t%d\tnull_mean\t%.2f\tnull_sd\t%.2f\tratio\t%.6f\tp\t%.6g\n"
        % (obs, nul.mean(), nul.std(ddof=1), obs / nul.mean(), p_iso))
g.write("[B]\tlag\tn_odd0\tn_odd1\tn_odd2\tfrac_odd1\tmean_d_odd0\tmean_d_odd1\tmean_d_odd2\n")
for r in rows_b:
    g.write("\t%d\t%d\t%d\t%d\t%.6f\t%.4f\t%.4f\t%.4f\n" % r)
g.write("\treweighted_lag1_with_lag2_mix\t%.6f\texplained_pct\t%.4f\n"
        % (float((m1 * f2).sum()), 100 * (base[0] - float((m1 * f2).sum())) / (base[0] - base[1])))
g.write("[hap]\tn\t%d\tfrac_pos\t%.6f\tmedian\t%.6f\tsign_p\t%.6g\twilcoxon_p\t%.6g\n"
        % (len(sg), (sg > 0).mean(), np.median(sg), bt.pvalue, w.pvalue))
g.close()
print("\n산출 :", os.path.join(OUT, "X5_lag2_noncircular.tsv"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('X5_lag2_noncircular.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
