# -*- coding: utf-8 -*-
r"""P08 XA - how robust is the decay constant a
   a is the number obtained 'when an exponential model is fitted to the semivariance at distances 3-15'.
   Changing the model, the fitting window and the target curve shows how much that number moves.
   (a) fit start distance: 1 / 2 / 3 / 4 / 5
   (b) target curve: semivariance (squared) vs mean difference (bp)
   (c) fit end distance: 10 / 12 / 15
   (d) nugget fixed or not
   Verdict: if a moves by orders of magnitude, a three-digit number such as '7.4' must not be written in the main text"""
import os, sys, io, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
SEED, B = 20261010, 2000
rng = np.random.default_rng(SEED)
import v11_io as G
from scipy import optimize

U, _, _ = G.load("h1")
VPU = sum(len(t[4]) for v in U.values() for t in v) / sum(len(v) for v in U.values())
print("§0 자 : var_per_unit = %.6f" % VPU)
assert abs(VPU - 6.866) < 0.01
NCOPY = {("%s" % k[0], "%s" % k[1]): len(v) for k, v in U.items()}
with open(os.path.join(HERE, "OUT_P04", "S3_matrices.pkl"), "rb") as f:
    S3 = pickle.load(f)
MATS, APS = S3["mats"], S3["aps"]
K = 15
KS = np.arange(1, K + 1)
SL = []
for s in sorted(MATS):
    aps = APS[s]
    N = NCOPY.get(s, 0)
    sel = [i for i, ap in enumerate(aps) if 0.2 * N < ap <= 0.8 * N]
    if len(sel) < 22:
        continue
    ii = np.array(sel)
    SL.append((MATS[s][np.ix_(ii, ii)].astype(np.int64), len(sel)))
PD = np.zeros((len(SL), K))
PQ = np.zeros((len(SL), K))
PN = np.zeros((len(SL), K))
for i, (M, n) in enumerate(SL):
    for ki, k in enumerate(KS):
        if n - k <= 0:
            continue
        v = M[np.arange(n - k), np.arange(k, n)]
        PD[i, ki] = v.sum()
        PQ[i, ki] = (v * v).sum()
        PN[i, ki] = n - k
ON = PN.sum(0)
GAM = 0.5 * PQ.sum(0) / ON
MEA = PD.sum(0) / ON
print("haplotype %d · semivariance %.2f → %.2f · 평균 차이 %.3f → %.3f"
      % (len(SL), GAM[0], GAM[-1], MEA[0], MEA[-1]))


def fit(y, w, start, end, free_nugget=True):
    m = (KS >= start) & (KS <= end)
    x, yy, ww = KS[m].astype(float), y[m], w[m]

    def f(p):
        cn = p[0] if free_nugget else yy.min()
        return np.sqrt(ww) * (cn + p[1] * (1 - np.exp(-x / max(p[2], 1e-6))) - yy)
    r = optimize.least_squares(f, [yy.min(), max(yy.max() - yy.min(), 1e-6), 5.0],
                               bounds=([0, 0, 0.1], [np.inf, np.inf, 500]))
    return r.x[2]


def boot_a(y_src, start, end, nrep=B):
    out = []
    for b in range(nrep):
        r = rng.integers(0, len(SL), len(SL))
        nn = PN[r].sum(0)
        y = (0.5 * PQ[r].sum(0) / nn) if y_src == "gamma" else (PD[r].sum(0) / nn)
        try:
            out.append(fit(y, nn, start, end))
        except Exception:
            pass
    out = np.array(out)
    return np.percentile(out, [2.5, 97.5]), float((out > 499).mean())


print("\n" + "=" * 78)
print("[⒜⒞] 적합 창을 바꾸면 a 가 어떻게 되나 (대상 = semivariance · nugget 자유)")
print("    %12s %10s %24s %10s" % ("적합 창", "a", "95%CI (블록 부트)", "상한붙음"))
rows = []
for st, en in ((1, 15), (2, 15), (3, 15), (4, 15), (5, 15), (3, 10), (3, 12)):
    a_ = fit(GAM, ON, st, en)
    ci, hit = boot_a("gamma", st, en)
    print("    %12s %10.3f %24s %10.3f" % ("거리 %d-%d" % (st, en), a_, "[%.2f, %.2f]" % (ci[0], ci[1]), hit))
    rows.append(("gamma", st, en, True, float(a_), float(ci[0]), float(ci[1]), hit))
print("\n[⒝] 대상 곡선을 평균 차이(bp)로 바꾸면 (같은 적합 창)")
print("    %12s %10s %24s %10s" % ("적합 창", "a", "95%CI", "상한붙음"))
for st, en in ((1, 15), (3, 15), (5, 15)):
    a_ = fit(MEA, ON, st, en)
    ci, hit = boot_a("mean", st, en)
    print("    %12s %10.3f %24s %10.3f" % ("거리 %d-%d" % (st, en), a_, "[%.2f, %.2f]" % (ci[0], ci[1]), hit))
    rows.append(("mean", st, en, True, float(a_), float(ci[0]), float(ci[1]), hit))
print("\n[⒟] nugget 을 고정하면 (semivariance · 거리 3-15)")
a_ = fit(GAM, ON, 3, 15, free_nugget=False)
print("    nugget 고정 : a = %.3f" % a_)
rows.append(("gamma", 3, 15, False, float(a_), np.nan, np.nan, np.nan))
av = np.array([r[4] for r in rows if r[0] == "gamma" and r[3]])
print("\n" + "=" * 78)
print("semivariance 판 a 범위 : %.2f ~ %.2f (최대/최소 %.2f 배)" % (av.min(), av.max(), av.max() / av.min()))
print("평균차이 판 a 범위     : %.2f ~ %.2f"
      % (min(r[4] for r in rows if r[0] == "mean"), max(r[4] for r in rows if r[0] == "mean")))
print("부트스트랩 CI 는 거리 3-15 기준 [%.2f, %.2f] — 상한/하한 %.2f 배"
      % (rows[2][5], rows[2][6], rows[2][6] / rows[2][5]))
g = io.open(os.path.join(OUT, "XA_decay_constant_robustness.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 XA · 감쇠 상수 a 의 민감도 · seed %d · B %d · 2026-10-10\n" % (SEED, B))
g.write("curve\tstart\tend\tfree_nugget\ta\tlo\thi\tfrac_at_bound\n")
for r in rows:
    g.write("%s\t%d\t%d\t%s\t%.4f\t%s\t%s\t%s\n"
            % (r[0], r[1], r[2], r[3], r[4],
               "" if r[5] != r[5] else "%.4f" % r[5], "" if r[6] != r[6] else "%.4f" % r[6],
               "" if r[7] != r[7] else "%.4f" % r[7]))
g.close()
print("\n산출 :", os.path.join(OUT, "XA_decay_constant_robustness.tsv"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('XA_decay_constant_robustness.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
