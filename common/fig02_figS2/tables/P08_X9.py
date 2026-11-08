# -*- coding: utf-8 -*-
r"""P08 X9, the 'no exchange floor' of patch length with the same pipeline as X4
   X4 ran only within-array distances. To draw a floor line in Fig4b, T_hat for unrelated
   individual pairs must be produced with the same estimator and the same settings (Q4b values use a different sampling setup and cannot be reused as is)."""
import os, sys, io, collections, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
SEED, EDGE, NP_ = 20261010, 5, 4000
rng = np.random.default_rng(SEED)
import v11_io as G
from scipy import optimize

CAN = G.CAN
LU = len(CAN)
U, _, _ = G.load("h1")
VPU = sum(len(t[4]) for v in U.values() for t in v) / sum(len(v) for v in U.values())
print("§0 자 : var_per_unit = %.6f" % VPU)
assert abs(VPU - 6.866) < 0.01
NCOPY = {("%s" % k[0], "%s" % k[1]): len(v) for k, v in U.items()}
FINE = {}
for Ln in io.open((__import__("os").environ["RDNA_DATA_ROOT"] + r"/v11_readvar/05_reference/1kGP_3202_ped_population_2026-09-23.txt"), encoding="utf8"):
    f = Ln.split()
    if len(f) >= 7 and f[1] != "SampleID":
        FINE[f[1]] = f[5]
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
INT = {}
for s, aps in bySlot.items():
    N = NCOPY.get(s, 0)
    v = [ap for ap in aps if EDGE < ap <= N - EDGE]
    if len(v) >= 9:
        INT[s] = v
HSL = sorted(INT)
bySamp = collections.defaultdict(list)
for s in HSL:
    bySamp[s[0]].append(s)
PAIRED = [tuple(v) for v in bySamp.values() if len(v) == 2]
print("interior haplotype %d · 개체 %d · 상동 개체 %d" % (len(HSL), len(bySamp), len(PAIRED)))


def events(k1, k2):
    A1, A2 = A_OF[k1], A_OF[k2]
    pos = sorted(set(A1) | set(A2))
    lab = [1 if (p in A1 and A2.get(p) == A1[p]) else 0 for p in pos]
    out = []
    i, M = 0, len(pos)
    while i < M:
        if lab[i] == 0:
            i += 1
            continue
        j = i
        while j + 1 < M and lab[j + 1] == 1:
            j += 1
        lo = pos[j] - pos[i]
        hi = (pos[j + 1] if j + 1 < M else LU) - (pos[i - 1] if i - 1 >= 0 else 1)
        out.append((lo, hi, 0 if (i - 1 >= 0 and j + 1 < M) else 1))
        i = j + 1
    return out


def S3f(t, T):
    return np.exp(-t / T)


def S2f(t, T):
    return (1 + 2 * t / T) * np.exp(-2 * t / T)


def fit(lo, hi, cen, S):
    def nll(T):
        T = max(float(T), 1e-3)
        pr = np.where(cen == 1, S(lo, T), S(lo, T) - S(hi, T))
        return -np.log(np.clip(pr, 1e-300, 1.0)).sum()
    r = optimize.minimize_scalar(nll, bounds=(1, 20000), method="bounded")
    b = r.fun

    def g(T):
        return nll(T) - b - 1.920729
    try:
        a_ = optimize.brentq(g, 1, r.x)
    except Exception:
        a_ = np.nan
    try:
        c_ = optimize.brentq(g, r.x, 20000)
    except Exception:
        c_ = np.nan
    return r.x, a_, c_


def sample(kind):
    pl = []
    if kind == "homolog":
        while len(pl) < NP_:
            a, b = PAIRED[rng.integers(0, len(PAIRED))]
            pl.append(((a[0], a[1], INT[a][rng.integers(0, len(INT[a]))]),
                       (b[0], b[1], INT[b][rng.integers(0, len(INT[b]))])))
    else:
        want_same = (kind == "unrelated_same_pop")
        while len(pl) < NP_:
            i, j = rng.integers(0, len(HSL), 2)
            if i == j or HSL[i][0] == HSL[j][0]:
                continue
            fi, fj = FINE.get(HSL[i][0], "?"), FINE.get(HSL[j][0], "?")
            if kind != "unrelated_all":
                if fi == "?" or fj == "?" or ((fi == fj) != want_same):
                    continue
            a, b = HSL[i], HSL[j]
            pl.append(((a[0], a[1], INT[a][rng.integers(0, len(INT[a]))]),
                       (b[0], b[1], INT[b][rng.integers(0, len(INT[b]))])))
    return pl


print("\n%-24s %9s %9s %12s %22s %12s" % ("갈래", "쌍", "사건", "T (지수)", "95%CI", "T (gamma2)"))
rows = []
for kind, nm in (("homolog", "상동 (같은 개인)"), ("unrelated_same_pop", "무관 · 같은 세부집단"),
                 ("unrelated_diff_pop", "무관 · 다른 세부집단"), ("unrelated_all", "무관 · 전체")):
    pl = sample(kind)
    ev = [e for k1, k2 in pl for e in events(k1, k2)]
    lo = np.array([e[0] for e in ev], float)
    hi = np.array([e[1] for e in ev], float)
    cen = np.array([e[2] for e in ev])
    T3, l3, h3 = fit(lo, hi, cen, S3f)
    T2, l2, h2 = fit(lo, hi, cen, S2f)
    print("%-24s %9d %9d %12.2f %22s %12.2f" % (nm, len(pl), len(ev), T3, "[%.1f, %.1f]" % (l3, h3), T2))
    rows.append((kind, nm, len(pl), len(ev), float(T3), float(l3), float(h3), float(T2), float(cen.mean())))
a1 = 519.14
print("\narray 내 거리 1 = %.2f bp ↔ 교환 없는 바닥(무관 전체) = %.2f bp ⇒ 비 %.4f"
      % (a1, rows[-1][4], a1 / rows[-1][4]))
print("참고 : 최장 공유동일 구간은 1130.4 ↔ 664.7 bp ⇒ 비 1.70 — 패치 추정량이 바닥 대비 신호가 더 크다")
g = io.open(os.path.join(OUT, "X9_patch_floor.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 X9 · 패치 길이의 '교환 없는 바닥' · X4 와 같은 설정 (interior · 쌍 %d) · seed %d\n" % (NP_, SEED))
g.write("kind\tlabel\tn_pair\tn_event\tT_exp\tlo\thi\tT_gamma2\tcensored\n")
for r in rows:
    g.write("%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.6f\n" % r)
g.write("#\narray_lag1_T\t%.4f\tfloor_T\t%.4f\tratio\t%.4f\n" % (a1, rows[-1][4], a1 / rows[-1][4]))
g.close()
print("\n산출 :", os.path.join(OUT, "X9_patch_floor.tsv"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('X9_patch_floor.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
