# -*- coding: utf-8 -*-
"""
EXPORT_origin_stability: stabilizes the inherited-type call so that it does not depend on a single lineage partition (KMeans K=8, seed 0).

Background (measured 2026-09-18)
  From v11 to v11-fix2, one more haplotype and one more background marker (45 -> 46) moved the KMeans 8-cluster boundaries,
  and the call changed for 68 of 838 variants (8%). The terminal-type constituent variants (635CT, DB43:2, 113GA, 380TA) and 2206GA
  had p=0.002, but their effect size (top-lineage occupancy - null) fell below 0.10, so inherited -> undetermined (small effect).
  -> The S25 strand-asymmetry main estimate flipped because of the single site 2206GA.

Method
  15 lineage partitions = KMeans K in {6, 8, 10} x random_state in {0..4} (n_init 20), on the 0/1 matrix of background markers
  (MID carrier rate 20-80%) over background haplotypes (MID >= 15 units). Self-reference guard: if the variant is itself a background marker, it is removed and the partition is rebuilt (same as the original call).
  Per partition: top-lineage occupancy, median of the null (500 random draws of as many haplotypes as carriers, the same permutation matrix for every partition), p = (null >= observed + 1)/501
             pass = p < 0.05 and (occupancy - null) >= 0.10
  frac_pass = fraction of partitions that pass, frac_sig = fraction of partitions with p<0.05, eff_median, p_median
Calls (stable version)
  ancestral            carrier haplotype fraction >= 0.90 (unchanged)
  inherited            frac_pass >= 0.5 (passes in a majority of partitions)
  undetermined (small effect)  frac_pass < 0.5 and frac_sig >= 0.5
  recurrence-suggestive  p_median >= 0.20, copies per array <= 1.5, recurrence context
  dispersed            p_median >= 0.20
  undetermined         otherwise (0.05 <= p_median < 0.20 and so on)
  stability    stable = frac_pass <= 0.2 or >= 0.8, borderline = in between
  Broad inherited set for sensitivity: frac_pass >= 0.2 (used by ORIGIN_MODE=broad in S25/S25b)

Usage
  python EXPORT_origin_stability.py F2   -> updates export\V11_ORIGIN_CALL.tsv in this folder to the stable version
                                           (originals kept: V11_ORIGIN_CALL.single_kmeans.tsv, V11_VARIANT_GROUPS.tsv.single_kmeans)
  python EXPORT_origin_stability.py v11  -> same computation on the original v11 data, writes only to export\stability_v11\ (for comparing agreement between the two versions)
"""
import sys, io, os, collections, shutil
MODE = sys.argv[1] if len(sys.argv) > 1 else "F2"
if MODE == "v11":
    sys.path.insert(0, __import__("os").environ["V11_0913_CODE"])
import numpy as np
from sklearn.cluster import KMeans
import unitdef_v11 as FC
import v11_io as G

F2DIR = G.OUT if MODE == "F2" else (__import__("os").environ["RDNA_DATA_ROOT"] + r"\5S_position_variant_work\2026-09-17_v11_fix2")   # correction: F2 mode output follows G.OUT (so the mirror guard applies). Previously the canonical path was hard-coded and a mirror run overwrote the canonical output.
SRC_EXPORT = os.path.join(G.OUT, "export")
OUTD = os.path.join(F2DIR, "export") if MODE == "F2" else os.path.join(F2DIR, "export", "stability_v11")
os.makedirs(OUTD, exist_ok=True)
print("판본 %s · 자료 %s · 출력 %s" % (MODE, G.OUT, OUTD))
KS = (6, 8, 10); SEEDS = (0, 1, 2, 3, 4); B = 500; EFF_MIN = 0.10; MINHAP = 8

U, A, M = G.load("primary")
rng = np.random.default_rng(20260918)

HAPS = []; MIDSET = {}
for k, rs in U.items():
    mid = [r for r in rs if r[0] > 30 and r[1] > 30]
    if len(mid) < 15: continue
    HAPS.append(k)
    c = collections.Counter()
    for r in mid: c.update(r[4])
    MIDSET[k] = (c, len(mid))
NH = len(HAPS)
have = collections.Counter()
for k in HAPS: have.update(MIDSET[k][0].keys())
MK = sorted([x for x in have if 0.20 * NH <= have[x] <= 0.80 * NH]); MKI = {x: i for i, x in enumerate(MK)}
X = np.zeros((NH, len(MK)), dtype=np.float32)
for i, k in enumerate(HAPS):
    for x in MIDSET[k][0]:
        j = MKI.get(x)
        if j is not None: X[i, j] = 1.0
print("배경 하플로타입 %d · 배경 마커 %d" % (NH, len(MK)))

PERM = np.array([rng.permutation(NH) for _ in range(B)])      # all partitions use the same permutations
def null_table(lab, K):
    """top lineage count when drawing n haplotypes: (B, NH) array, cum[:, n-1]"""
    oh = np.eye(K, dtype=np.int16)[lab]                        # NH×K
    cum = np.cumsum(oh[PERM], axis=1)                          # B×NH×K
    return cum.max(2)                                          # B×NH
PART = {}
for K in KS:
    for s in SEEDS:
        lab = KMeans(n_clusters=K, n_init=20, random_state=s).fit(X).labels_
        PART[(K, s)] = (lab, null_table(lab, K))
print("분할 %d 벌 준비" % len(PART))
EXCL_CACHE = {}
def part_excl(x, K, s):
    key = (x, K, s)
    if key not in EXCL_CACHE:
        keep = [j for j in range(len(MK)) if j != MKI[x]]
        lab = KMeans(n_clusters=K, n_init=20, random_state=s).fit(X[:, keep]).labels_
        EXCL_CACHE[key] = (lab, null_table(lab, K))
    return EXCL_CACHE[key]

carry_all = collections.defaultdict(set); nu = collections.Counter(); nh_all = collections.defaultdict(set)
for k, rs in U.items():
    s = set()
    for r in rs:
        for x in r[4]: nu[x] += 1; nh_all[x].add(k)
        s |= r[4]
    if k in MIDSET:
        for x in s: carry_all[x].add(k)
HI = {k: i for i, k in enumerate(HAPS)}
NHAPTOT = len(U)
cand = sorted([x for x in carry_all if len(carry_all[x]) >= MINHAP], key=lambda z: -nu[z])
print("대상 변이 %d" % len(cand))

rows = []
for i, x in enumerate(cand):
    idx = np.array(sorted(HI[k] for k in carry_all[x])); n = len(idx)
    hapfrac = len(nh_all[x]) / float(NHAPTOT); per = nu[x] / float(len(nh_all[x]))
    cl, pos, ref, alt, ln, blk, rep = FC.parse(x)
    recur = ((cl == "snv" and (ref, alt) in (("C", "T"), ("G", "A")) and FC.iscpg(pos)) or (cl in ("del", "ins") and bool(rep or blk)))
    tops = []; effs = []; ps = []; passes = []
    if hapfrac < 0.90:
        for (K, s) in PART:
            lab, NT = part_excl(x, K, s) if x in MKI else PART[(K, s)]
            top = np.bincount(lab[idx], minlength=K).max() / float(n)
            nul = NT[:, n - 1] / float(n)
            p = (np.sum(nul >= top) + 1) / (B + 1.0); eff = top - np.median(nul)
            tops.append(top); effs.append(eff); ps.append(p); passes.append(p < 0.05 and eff >= EFF_MIN)
    if hapfrac >= 0.90:
        call = "조상형"; fp = fs = em = pm = float("nan"); stab = ""
    else:
        fp = float(np.mean(passes)); fs = float(np.mean(np.array(ps) < 0.05)); em = float(np.median(effs)); pm = float(np.median(ps))
        if fp >= 0.5: call = "계승형"
        elif fs >= 0.5: call = "판정불가(효과작음)"
        elif pm >= 0.20 and per <= 1.5 and recur: call = "재발시사형"
        elif pm >= 0.20: call = "분산형"
        else: call = "판정불가"
        stab = "안정" if (fp <= 0.2 or fp >= 0.8) else "경계"
    rows.append((x, n, len(nh_all[x]), nu[x], hapfrac, per, recur, call, fp, fs, em, pm, stab,
                 (min(effs), max(effs)) if effs else (float("nan"), float("nan"))))
    if (i + 1) % 200 == 0: print("   ... %d / %d" % (i + 1, len(cand)))

def f4(v): return "" if v != v else "%.4g" % v
p = os.path.join(OUTD, "V11_ORIGIN_STABILITY.tsv")
g = io.open(p, "w", encoding="utf8", newline="\n")
g.write("variant\tn_carrier_bg_hap\tn_haplotypes\tn_units\thap_freq\tunits_per_carrier_array\trecurrence_context\t"
        "origin_call_stable\tfrac_pass\tfrac_sig\teff_median\teff_min\teff_max\tp_median\tstability\n")
for r in rows:
    g.write("%s\t%d\t%d\t%d\t%.4f\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        r[0], r[1], r[2], r[3], r[4], r[5], "Y" if r[6] else "N", r[7], f4(r[8]), f4(r[9]), f4(r[10]),
        f4(r[13][0]), f4(r[13][1]), f4(r[11]), r[12]))
g.close(); print("저장:", p)
ST = {r[0]: r for r in rows}

# ---------- comparison with the single-partition call ----------
single = {}
sp = os.path.join(SRC_EXPORT, "V11_ORIGIN_CALL.single_kmeans.tsv")
if not os.path.exists(sp): sp = os.path.join(SRC_EXPORT, "V11_ORIGIN_CALL.tsv")
for j, L in enumerate(io.open(sp, encoding="utf8")):
    f = L.rstrip("\n").split("\t")
    if j == 0: H = {c: q for q, c in enumerate(f)}; continue
    single[f[0]] = f[H["origin_call"]]
tab = collections.Counter((single.get(x, "(없음)"), ST[x][7]) for x in ST)
print("\n=== 한 벌(K8·seed0) 판정 → 안정판 ===")
for (a, b), n in sorted(tab.items(), key=lambda t: -t[1]): print("  %-14s → %-14s %4d" % (a, b, n))
cnt = collections.Counter(r[7] for r in rows)
print("\n=== 안정판 분포 ===  " + " · ".join("%s %d" % (c, n) for c, n in cnt.most_common()))
sc = collections.Counter(r[12] for r in rows if r[12])
print("안정도 : " + " · ".join("%s %d" % (c, n) for c, n in sc.items()) +
      " · 넓은 계승 집합(frac_pass≥0.2) %d" % sum(1 for r in rows if r[8] == r[8] and r[8] >= 0.2))
print("\n=== 확인 대상 변이 ===")
for x in ("2205TC", "2235TG", "2206GA", "635CT", "DB43:2", "113GA", "380TA", "1164GA", "740GA", "DB16:17", "D1096:1", "1250CT", "R1498TA", "2210GT"):
    if x in ST:
        r = ST[x]; print("  %-9s 보유 %4d · %-16s frac_pass %s · frac_sig %s · 효과 중앙 %s [%s, %s] · p중앙 %s · %s (한 벌 : %s)" % (
            x, r[2], r[7], f4(r[8]), f4(r[9]), f4(r[10]), f4(r[13][0]), f4(r[13][1]), f4(r[11]), r[12], single.get(x, "-")))

if MODE == "F2":
    # ---------- keep the originals, then replace origin_call with the stable version ----------
    for name in ("V11_ORIGIN_CALL.tsv", "V11_VARIANT_GROUPS.tsv"):
        src = os.path.join(OUTD, name)
        bak = os.path.join(OUTD, name.replace(".tsv", ".single_kmeans.tsv"))
        if not os.path.exists(bak): shutil.copy2(src, bak); print("원본 보존:", bak)
        hd = None; out = []
        for L in io.open(bak, encoding="utf8"):
            f = L.rstrip("\n").split("\t")
            if hd is None:
                hd = f; H = {c: q for q, c in enumerate(hd)}
                out.append(hd + ["origin_call_single_k8s0", "frac_pass", "frac_sig", "eff_median", "p_median", "stability"]); continue
            r = ST.get(f[0]); old = f[H["origin_call"]]
            if r is not None:
                f[H["origin_call"]] = r[7]
                out.append(f + [old, f4(r[8]), f4(r[9]), f4(r[10]), f4(r[11]), r[12]])
            else:
                out.append(f + [old, "", "", "", "", ""])
        g = io.open(src, "w", encoding="utf8", newline="\n")
        for f in out: g.write("\t".join(f) + "\n")
        g.close(); print("갱신:", src)
