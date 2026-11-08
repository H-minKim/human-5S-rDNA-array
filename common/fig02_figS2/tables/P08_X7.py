# -*- coding: utf-8 -*-
r"""P08 X7, a measurement that separates "is it really homogenization, or just similarity"
   What the order permutation null already controls: overall similarity and each unit's variant composition (units are shuffled whole)
   What the order permutation does not control:
     (a) recurrent mutation (hotspot), two units may have gained a variant independently at the same site
     (b) assembly artefact, the assembler may have placed similar units next to each other
     (c) local duplication (unequal crossing over), neighbours may be similar because they are recent duplicates (also homogenization, but a different mechanism)
   The test that separates (a) = sharing of rare variants. Recurrent mutations occur at common sites.
      If two units share a private variant seen only within one array, recurrence can hardly explain it.
      If the units carrying that private variant lie close to each other = direct evidence of local exchange.
   (d) sensitivity to the divergence quartile threshold (10/25/33/50%)"""
import os, sys, io, collections, pickle, time
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P08")
SEED, N_PERM, DEC = 20261010, 999, 10
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
with open(os.path.join(HERE, "OUT_P04", "Q0_alignment.pkl"), "rb") as f:
    DIFF = pickle.load(f)["diff"]
with open(os.path.join(HERE, "OUT_P04", "S3_matrices.pkl"), "rb") as f:
    S3 = pickle.load(f)
MATS, APS = S3["mats"], S3["aps"]
keys = [k for k in DIFF if k not in S3["bad"]]
A_OF = {}
for k in keys:
    sub, dele, ins, ul = DIFF[k]
    d = dict(sub)
    for p in dele:
        d[p] = "-"
    A_OF[k] = d
# Cohort-wide frequency + number of haplotypes carrying the variant
CF = collections.Counter()
HAPS_OF = collections.defaultdict(set)
for k in keys:
    for q, b in A_OF[k].items():
        CF[(q, b)] += 1
        HAPS_OF[(q, b)].add((k[0], k[1]))
NT = len(keys)
NH_ALL = len({(k[0], k[1]) for k in keys})
print("유닛 %d · haplotype %d · 서로 다른 (자리,대립) %d" % (NT, NH_ALL, len(CF)))

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
    kk = [(s[0], s[1], aps[i]) for i in sel]
    dec = np.minimum((np.arange(n) * DEC) // n, DEC - 1)
    PACK.append((s, M, kk, n, dec))
print("가운데 창 haplotype %d · 유닛 %d" % (len(PACK), sum(x[3] for x in PACK)))

# ---------------- (a)(b)(c) distance between variant-carrying units, by frequency stratum ----------------
print("\n" + "=" * 78)
print("[핵심] 변이를 공유한 유닛들이 array 위에서 얼마나 가까운가 — 코호트 희귀도 층별")
print("    귀무 = 같은 carrier 수를 위치 십분위 층화해 무작위 배치 (N = %d)" % N_PERM)
STRATA = [("사유 (그 haplotype 에만)", lambda q, b: len(HAPS_OF[(q, b)]) == 1),
          ("매우 희귀 (<0.1%)", lambda q, b: CF[(q, b)] / NT < 0.001),
          ("희귀 (0.1-1%)", lambda q, b: 0.001 <= CF[(q, b)] / NT < 0.01),
          ("중간 (1-10%)", lambda q, b: 0.01 <= CF[(q, b)] / NT < 0.10),
          ("흔함 (>=10%)", lambda q, b: CF[(q, b)] / NT >= 0.10)]
SETS = {nm: [] for nm, f in STRATA}
for s, M, kk, n, dec in PACK:
    carry = collections.defaultdict(list)
    for i in range(n):
        for q, b in A_OF[kk[i]].items():
            carry[(q, b)].append(i)
    for (q, b), v in carry.items():
        if len(v) < 2:
            continue
        for nm, f in STRATA:
            if f(q, b):
                SETS[nm].append((len(PACK) and PACK.index((s, M, kk, n, dec)) if False else s, np.array(v)))
                break
# Back to slot index
SIDX = {x[0]: i for i, x in enumerate(PACK)}


def nn_stats(groups, orders=None):
    """groups = [(slot, idx array)], if orders = {slot: permutation array} apply it"""
    tot, cnt = 0.0, 0
    for s, v in groups:
        w = v if orders is None else orders[s][v]
        pos = np.sort(w)
        if len(pos) < 2:
            continue
        d = np.diff(pos)
        dmin = np.minimum(np.concatenate([[1e9], d]), np.concatenate([d, [1e9]]))
        tot += dmin.sum()
        cnt += len(pos)
    return tot / cnt if cnt else np.nan, cnt


rows = []
t0 = time.time()
for nm, f in STRATA:
    gs = SETS[nm]
    if len(gs) < 20:
        print("    %-22s 사례 %d — 너무 적어 건너뜀" % (nm, len(gs)))
        continue
    obs, ncar = nn_stats(gs)
    nul = np.zeros(N_PERM)
    for b in range(N_PERM):
        orders = {}
        for s, M, kk, n, dec in PACK:
            o = np.arange(n)
            for d in range(DEC):
                w = np.where(dec == d)[0]
                if len(w) > 1:
                    o[w] = w[rng.permutation(len(w))]
            orders[s] = o
        nul[b] = nn_stats(gs, orders)[0]
    p = (1 + int((nul <= obs).sum())) / (N_PERM + 1)
    print("    %-22s 변이 %6d · carrier %7d · 최근접 거리 관측 %7.4f ↔ 귀무 %7.4f ± %.4f · 배율 %.4f · p = %.4g"
          % (nm, len(gs), ncar, obs, nul.mean(), nul.std(ddof=1), obs / nul.mean(), p))
    rows.append((nm, len(gs), ncar, float(obs), float(nul.mean()), float(nul.std(ddof=1)),
                 float(obs / nul.mean()), p))
print("    (%.0fs)  읽는 법 : 배율 < 1 = 관측이 귀무보다 가깝다 = 국소 교환" % (time.time() - t0))

# ---------------- distance distribution of private variant carrier pairs ----------------
print("\n[보강] 사유 변이 (그 haplotype 에만 있는 변이) carrier 쌍의 거리 분포")
pri = SETS["사유 (그 haplotype 에만)"]
dd = []
for s, v in pri:
    pos = np.sort(v)
    for a in range(len(pos)):
        for b2 in range(a + 1, len(pos)):
            dd.append(pos[b2] - pos[a])
dd = np.array(dd)
if len(dd):
    print("    carrier 쌍 %d · 거리 median %.1f · 1 copy 몫 %.4f · <=5 copy 몫 %.4f · <=10 몫 %.4f"
          % (len(dd), np.median(dd), (dd == 1).mean(), (dd <= 5).mean(), (dd <= 10).mean()))
    p2 = float("nan")
    nul2 = np.array([np.nan])
    print("    (쌍 거리의 별도 귀무는 생략 — 위 층별 최근접 귀무가 같은 것을 이미 검정한다)")

# ---------------- (d) sensitivity to the divergence quartile threshold ----------------
print("\n" + "=" * 78)
print("[민감도] 발산 임계값을 바꿔도 lag1 − lag2 가 0 으로 가는가")
print("    %8s %10s %12s %12s %12s" % ("임계", "제외 몫", "lag1", "lag2", "lag1 − lag2"))
rows_d = []
for pct in (0, 50, 67, 75, 90):
    num = np.zeros(2)
    den = np.zeros(2)
    ndrop = 0
    for s, M, kk, n, dec in PACK:
        med = int(np.argmin(M.sum(1)))
        div = M[med].astype(float)
        ok = np.ones(n, bool) if pct == 0 else (div <= np.percentile(div, pct))
        ndrop += n - int(ok.sum())
        for k in (1, 2):
            i = np.arange(n - k)
            m = ok[i] & ok[i + k]
            if m.any():
                num[k - 1] += M[i[m], i[m] + k].sum()
                den[k - 1] += m.sum()
    c = num / den
    print("    %7s%% %10.4f %12.4f %12.4f %12.4f"
          % ("없음" if pct == 0 else "상위 %d" % (100 - pct), ndrop / sum(x[3] for x in PACK),
             c[0], c[1], c[0] - c[1]))
    rows_d.append((pct, ndrop / sum(x[3] for x in PACK), float(c[0]), float(c[1]), float(c[0] - c[1])))

g = io.open(os.path.join(OUT, "X7_homogenisation_evidence.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 X7 · ruler majref %.6f · seed %d · N_perm %d · 2026-10-10\n" % (VPU, SEED, N_PERM))
g.write("[freq]\tstratum\tn_variant\tn_carrier\tobs_nn\tnull_nn\tnull_sd\tratio\tp\n")
for r in rows:
    g.write("\t%s\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6g\n" % r)
if len(dd):
    g.write("[private]\tn_pair\t%d\tmedian\t%.2f\tfrac_adj\t%.6f\tfrac_le5\t%.6f\tfrac_le10\t%.6f\n"
            % (len(dd), np.median(dd), (dd == 1).mean(), (dd <= 5).mean(), (dd <= 10).mean()))
    g.write("\tobs_mean\t%.6f\tnull_mean\t%.6f\tnull_sd\t%.6f\tratio\t%.6f\tp\t%.6g\n"
            % (dd.mean(), nul2.mean(), nul2.std(ddof=1), dd.mean() / nul2.mean(), p2))
g.write("[sens]\tpct_keep\tfrac_dropped\tlag1\tlag2\tdiff\n")
for r in rows_d:
    g.write("\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n" % r)
g.close()
np.save(os.path.join(OUT, "X7_private_pair_dist.npy"), dd)
print("\n산출 :", os.path.join(OUT, "X7_homogenisation_evidence.tsv"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('X7_homogenisation_evidence.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
