# -*- coding: utf-8 -*-
r"""P08 X1 - population control: is the homolog vs unrelated difference explained by shared ancestry
   logic: unrelated pairs are likely from different populations, while the two haplotypes of one person come from the same population.
          => homologs being more alike may reflect a population constraint rather than Ohta C2.
   method ⒜ φ in levels: homolog / unrelated, same subpopulation / unrelated, same superpopulation / unrelated, different superpopulation
        ⒝ paired cluster bootstrap (resampling unit = individual) for CIs of ratios and differences
        ⒞ individual label permutation only within population (re-pair haplotypes of the same subpopulation) => null with population fixed
        ⒟ twin carrier rate (=0, <=5) recomputed with size-matched unrelated pairs from the same population
   first: array orientation (proximal/distal) gate, confirmed empirically with U84C"""
import os, sys, io, time, collections, pickle
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
from statsmodels.stats.proportion import proportion_confint

CAN = G.CAN
LU = len(CAN)
U, _, _ = G.load("h1")
VPU = sum(len(t[4]) for v in U.values() for t in v) / sum(len(v) for v in U.values())
print("§0 자 : var_per_unit = %.6f" % VPU)
assert abs(VPU - 6.866) < 0.01
NCOPY = {("%s" % k[0], "%s" % k[1]): len(v) for k, v in U.items()}
FINE, SUP = {}, {}
for Ln in io.open((__import__("os").environ["RDNA_DATA_ROOT"] + r"/v11_readvar/05_reference/1kGP_3202_ped_population_2026-09-23.txt"), encoding="utf8"):
    f = Ln.split()
    if len(f) >= 7 and f[1] != "SampleID":
        FINE[f[1]] = f[5]
        SUP[f[1]] = f[6]
RC = {}
for f, Hx in G._rows(G.SUMF):
    RC[(f[Hx["sample"]], f[Hx["slot"]])] = f[Hx["revcomp"]]
with open(os.path.join(HERE, "OUT_P04", "Q0_alignment.pkl"), "rb") as f:
    DIFF = pickle.load(f)["diff"]
with open(os.path.join(HERE, "OUT_P04", "S3_matrices.pkl"), "rb") as f:
    S3 = pickle.load(f)
BAD, MATS, APS = S3["bad"], S3["mats"], S3["aps"]
keys = [k for k in DIFF if k not in BAD]
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


bySlot = collections.defaultdict(list)
for k in keys:
    bySlot[(k[0], k[1])].append(k[2])
for s in bySlot:
    bySlot[s].sort()

# ================= ⒪ orientation gate : U84C =================
print("\n" + "=" * 78)
print("[⒪] array 방향 관문 — v11_io 가 쓰는 U84C 로 실측 확인")
print("    v11_io.load 가 찍는 기준 : r[0](근위)<=3 보유율 0.194 · r[1](원위)<=3 보유율 0.000")
print("    CAN[84] = %s  ⇒ 보유 = A_OF[84] == 'C'" % CAN[83])
for tag in ("Y", "N", "ALL"):
    lo_n = lo_d = hi_n = hi_d = 0
    for s, aps in bySlot.items():
        if tag != "ALL" and RC.get(s) != tag:
            continue
        N = NCOPY.get(s, 0)
        for ap in aps:
            k = (s[0], s[1], ap)
            has = (dict(SET[k]).get(84) == "C")
            if ap <= 3:
                lo_n += 1
                lo_d += has
            if ap >= N - 2:
                hi_n += 1
                hi_d += has
    print("    revcomp=%-3s  array_pos<=3 : %.4f (n %5d)   ·   array_pos>=N-2 : %.4f (n %5d)"
          % (tag, lo_d / max(lo_n, 1), lo_n, hi_d / max(hi_n, 1), hi_n))
print("    ⇒ 판정 : 0.194 에 가까운 끝이 근위 · 0.000 쪽이 원위")

# ================= cohort =================
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
SAMPS = sorted(bySamp)
PAIRED_S = [sp for sp in SAMPS if len(bySamp[sp]) == 2]
print("\nhaplotype %d · 개체 %d · 두 하플로 다 있는 개체 %d" % (len(HSL), len(SAMPS), len(PAIRED_S)))
fine = np.array([FINE.get(s[0], "?") for s in HSL])
sup = np.array([SUP.get(s[0], "?") for s in HSL])
print("대집단 분포 (haplotype) : %s" % sorted(collections.Counter(sup).items(), key=lambda x: -x[1]))
cf = collections.Counter(FINE.get(sp, "?") for sp in PAIRED_S)
print("두 하플로 보유 개체의 세부집단 : %d 종 · 상위 %s" % (len(cf), cf.most_common(8)))

sid = {}
rows = []
for s in HSL:
    c = collections.Counter(SIG[(s[0], s[1], ap)] for ap in INT[s])
    rows.append(c)
    for g in c:
        sid.setdefault(g, len(sid))
X = np.zeros((len(HSL), len(sid)), np.int32)
for i, c in enumerate(rows):
    for g, n in c.items():
        X[i, sid[g]] = n
SAME = X.astype(np.int64) @ X.T
nvec = X.sum(1).astype(np.int64)
SAMPARR = np.array([s[0] for s in HSL])

# ================= ⒜ φ in four levels =================
print("\n" + "=" * 78)
print("[⒜] φ 를 네 단계로 — 집단을 좁혀가며 상동과 비교")


def phi_of(mask_fn, idxs=None):
    I = np.arange(len(HSL)) if idxs is None else np.asarray(idxs)
    num = den = 0
    for a in range(len(I)):
        for b in range(a + 1, len(I)):
            i, j = I[a], I[b]
            if SAMPARR[i] == SAMPARR[j]:
                continue
            if not mask_fn(i, j):
                continue
            num += int(SAME[i, j])
            den += int(nvec[i]) * int(nvec[j])
    return num, den


t0 = time.time()
nh = dh = 0
for sp in PAIRED_S:
    a, b = bySamp[sp]
    i, j = IDX[a], IDX[b]
    nh += int(SAME[i, j])
    dh += int(nvec[i]) * int(nvec[j])
res = collections.OrderedDict()
res["상동 (같은 개인)"] = (nh, dh)
res["무관 · 같은 세부집단"] = phi_of(lambda i, j: fine[i] == fine[j] and fine[i] != "?")
res["무관 · 같은 대집단"] = phi_of(lambda i, j: sup[i] == sup[j] and sup[i] != "?")
res["무관 · 다른 대집단"] = phi_of(lambda i, j: sup[i] != sup[j] and sup[i] != "?" and sup[j] != "?")
res["무관 · 전체"] = phi_of(lambda i, j: True)
print("    %-22s %12s %16s %10s" % ("층", "φ", "쌍 수 (분모)", "상동/이것"))
for nm, (n_, d_) in res.items():
    print("    %-22s %12.6f %16d %10s"
          % (nm, n_ / d_, d_, "-" if nm.startswith("상동") else "%.4f" % ((nh / dh) / (n_ / d_))))
print("    (%.0fs)" % (time.time() - t0))

# ================= ⒝ paired cluster bootstrap =================
print("\n[⒝] 짝지은 cluster bootstrap (재표본 단위 = 개체 · 같은 재표본에서 네 양 전부 · B = %d)" % B_BOOT)
IND_HAPS = [[IDX[s] for s in bySamp[sp]] for sp in SAMPS]
IND_FINE = [FINE.get(sp, "?") for sp in SAMPS]
IND_SUP = [SUP.get(sp, "?") for sp in SAMPS]
nI = len(SAMPS)


def stats_from(indlist):
    haps, owner = [], []
    for t, ii in enumerate(indlist):
        for h in IND_HAPS[ii]:
            haps.append(h)
            owner.append(t)
    haps = np.asarray(haps)
    owner = np.asarray(owner)
    fi = np.array([IND_FINE[ii] for ii in indlist])[owner]
    su = np.array([IND_SUP[ii] for ii in indlist])[owner]
    nn = nvec[haps]
    S = SAME[np.ix_(haps, haps)]
    den = nn[:, None] * nn[None, :]
    diff = owner[:, None] != owner[None, :]
    out = {}
    nh_ = dh_ = 0
    for t, ii in enumerate(indlist):
        hs = IND_HAPS[ii]
        if len(hs) == 2:
            nh_ += int(SAME[hs[0], hs[1]])
            dh_ += int(nvec[hs[0]]) * int(nvec[hs[1]])
    out["hom"] = nh_ / dh_ if dh_ else np.nan
    for nm, m in (("fine", diff & (fi[:, None] == fi[None, :]) & (fi[:, None] != "?")),
                  ("sup", diff & (su[:, None] == su[None, :]) & (su[:, None] != "?")),
                  ("all", diff)):
        a_, b_ = int(S[m].sum()), int(den[m].sum())
        out[nm] = a_ / b_ if b_ else np.nan
    return out


BS = {k: [] for k in ("hom", "fine", "sup", "all")}
t0 = time.time()
for b in range(B_BOOT):
    r = rng.integers(0, nI, nI).tolist()
    o = stats_from(r)
    if any(np.isnan(v) for v in o.values()):
        continue
    for k in BS:
        BS[k].append(o[k])
BS = {k: np.array(v) for k, v in BS.items()}
print("    B = %d (%.0fs)" % (len(BS["hom"]), time.time() - t0))
obs = stats_from(list(range(nI)))
rows_b = []
for nm, key in (("무관 · 같은 세부집단", "fine"), ("무관 · 같은 대집단", "sup"), ("무관 · 전체", "all")):
    r_ = BS["hom"] / BS[key]
    d_ = BS["hom"] - BS[key]
    ci_r = np.percentile(r_, [2.5, 97.5])
    ci_d = np.percentile(d_, [2.5, 97.5])
    p = 2 * min((d_ <= 0).mean(), (d_ >= 0).mean())
    print("    상동 / %-20s = %.4f  95%%CI [%.4f, %.4f]  (1 포함 %s) · 차이 p = %.4g"
          % (nm, obs["hom"] / obs[key], ci_r[0], ci_r[1], "YES" if ci_r[0] <= 1 <= ci_r[1] else "NO", p))
    rows_b.append((nm, obs["hom"] / obs[key], ci_r[0], ci_r[1], obs["hom"] - obs[key], ci_d[0], ci_d[1], p))

# ================= ⒞ individual label permutation within population =================
print("\n[⒞] 개인 라벨 순열을 집단 안에서만 — 집단을 고정하고 '같은 개인' 관계만 무작위화 (N = %d)" % N_PERM)
for level, lab in (("세부집단", IND_FINE), ("대집단", IND_SUP)):
    grp = collections.defaultdict(list)
    for t, sp in enumerate(SAMPS):
        if len(bySamp[sp]) == 2 and lab[t] != "?":
            grp[lab[t]].append(t)
    usable = {g: v for g, v in grp.items() if len(v) >= 2}
    nind = sum(len(v) for v in usable.values())
    print("    %s : 쓸 수 있는 집단 %d 종 · 개체 %d (집단당 median %d)"
          % (level, len(usable), nind, int(np.median([len(v) for v in usable.values()]))))
    on = od = 0
    for g, v in usable.items():
        for t in v:
            h = IND_HAPS[t]
            on += int(SAME[h[0], h[1]])
            od += int(nvec[h[0]]) * int(nvec[h[1]])
    o_ = on / od
    nul = np.zeros(N_PERM)
    for b in range(N_PERM):
        num = den = 0
        for g, v in usable.items():
            hs = [h for t in v for h in IND_HAPS[t]]
            p = rng.permutation(len(hs))
            for a in range(0, len(hs) - 1, 2):
                i, j = hs[p[a]], hs[p[a + 1]]
                if SAMPARR[i] == SAMPARR[j]:
                    continue
                num += int(SAME[i, j])
                den += int(nvec[i]) * int(nvec[j])
        nul[b] = num / den if den else np.nan
    nul = nul[~np.isnan(nul)]
    pp = (1 + int((nul >= o_).sum())) / (len(nul) + 1)
    print("       관측 φ_homolog(그 집단들) = %.6f ↔ 같은 집단 안 귀무 %.6f ± %.6f · 배율 %.4f · p = %.4g"
          % (o_, nul.mean(), nul.std(ddof=1), o_ / nul.mean(), pp))
    rows_b.append(("순열-" + level, o_, nul.mean(), nul.std(ddof=1), o_ / nul.mean(), pp, 0, 0))

# ================= ⒟ twin carrier rate with population control =================
print("\n" + "=" * 78)
print("[⒟] 쌍둥이 보유율(=0, <=5) — 무관을 같은 집단으로 묶어 크기 맞춰 재계산 (반복 5)")
REPL = 5


def twin_rates(pairs_fn, tag):
    t0 = 0
    t5 = 0
    n_ = 0
    for rep in range(REPL):
        for a in HSL:
            b = pairs_fn(a)
            if b is None:
                continue
            ka = [(a[0], a[1], ap) for ap in INT[a]]
            kb = [(b[0], b[1], ap) for ap in INT[b]]
            M = np.empty((len(ka), len(kb)), np.int16)
            for i, k1 in enumerate(ka):
                for j, k2 in enumerate(kb):
                    M[i, j] = pd2(k1, k2)
            t0 += int((M.min(1) == 0).sum())
            t5 += int((M.min(1) <= 5).sum())
            n_ += len(ka)
    c0 = proportion_confint(t0, n_, 0.05, "wilson")
    c5 = proportion_confint(t5, n_, 0.05, "wilson")
    print("    %-26s unit %7d · =0 %6.2f%% [%.2f, %.2f] · <=5 %6.2f%% [%.2f, %.2f]"
          % (tag, n_, 100 * t0 / n_, 100 * c0[0], 100 * c0[1], 100 * t5 / n_, 100 * c5[0], 100 * c5[1]))
    return (tag, n_, t0, t0 / n_, c0[0], c0[1], t5, t5 / n_, c5[0], c5[1])


byfine = collections.defaultdict(list)
bysup = collections.defaultdict(list)
for i, s in enumerate(HSL):
    if fine[i] != "?":
        byfine[fine[i]].append(s)
    if sup[i] != "?":
        bysup[sup[i]].append(s)


def pick(pool, a):
    cand = [x for x in pool if x[0] != a[0]]
    return cand[rng.integers(0, len(cand))] if cand else None


rows_d = []
t0_ = time.time()
rows_d.append(twin_rates(lambda a: (bySamp[a[0]][1] if bySamp[a[0]][0] == a else bySamp[a[0]][0])
                         if len(bySamp[a[0]]) == 2 else None, "상동 (같은 개인)"))
rows_d.append(twin_rates(lambda a: pick(byfine.get(FINE.get(a[0], "?"), []), a), "무관 · 같은 세부집단"))
rows_d.append(twin_rates(lambda a: pick(bysup.get(SUP.get(a[0], "?"), []), a), "무관 · 같은 대집단"))
rows_d.append(twin_rates(lambda a: pick([x for x in HSL if SUP.get(x[0], "?") != SUP.get(a[0], "?")], a),
                         "무관 · 다른 대집단"))
rows_d.append(twin_rates(lambda a: pick(HSL, a), "무관 · 전체"))
print("    (%.0fs · 반복 %d 평균 · 비교집합은 하플로타입 하나로 크기 맞춤)" % (time.time() - t0_, REPL))

g = io.open(os.path.join(OUT, "X1_population_control.tsv"), "w", encoding="utf8", newline="\n")
g.write("# P08 X1 집단 통제 · ruler majref %.6f · seed %d · B %d · N_perm %d · 2026-10-10\n"
        % (VPU, SEED, B_BOOT, N_PERM))
g.write("[A]\tstratum\tphi\tnum\tden\tratio_hom_over\n")
for nm, (n_, d_) in res.items():
    g.write("\t%s\t%.8f\t%d\t%d\t%.6f\n" % (nm, n_ / d_, n_, d_, (nh / dh) / (n_ / d_)))
g.write("[B]\tcomparison\tv1\tv2\tv3\tv4\tv5\tv6\tv7\n")
for r in rows_b:
    g.write("\t%s\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8g\n" % r)
g.write("[D]\ttier\tn_unit\ttwin0\tfrac0\tlo0\thi0\ttwin5\tfrac5\tlo5\thi5\n")
for r in rows_d:
    g.write("\t%s\t%d\t%d\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%.6f\t%.6f\n" % r)
g.close()
np.savez_compressed(os.path.join(OUT, "X1_boot.npz"), **BS)
print("\n산출 :", os.path.join(OUT, "X1_population_control.tsv"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('X1_population_control.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
