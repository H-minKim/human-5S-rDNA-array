# -*- coding: utf-8 -*-
"""
S26b (new 3-way split) - split variants into 'inherited' and 'newly arisen', look at their distribution within the array,
      and measure which side the variants that make the boundaries found per array (S20 change points) belong to. (v11, primary)

(haplotype, variant) origin classes, called from V11_ORIGIN_CALL.tsv (lineage test) and carrier extent
  ancestral        origin_call = ancestral (all haplotypes, nearly all units)
  inherited        origin_call = inherited (passed down along the lineage background)
  other shared     several haplotypes, >=2 units in that array, not called ancestral/inherited (dispersed, recurrence-suggestive, undecidable, untested)
  private spread   only in this haplotype in the cohort and >=2 units in the array (arose in that lineage and then spread to neighbours)
  rare private (1 unit)   1 unit in the array, only in this haplotype in the cohort (proxy for the most recent events)
  rare shared (1 unit)    1 unit in the array, also present in other haplotypes, not ancestral/inherited (recurrence or a remnant of an old variant)
  'Rare (1 unit)' is not de novo. In 5 trio pairs on 2026-09-20, 97% of the rare variants of this class were also in the transmitting parent.
  How to read: not 'newly arisen variants' but 'where rare variants remain'.

(1) Position distribution: count per unit by class, proximal 1-30 | middle (>30 on both sides) | distal 30-1, two-stage aggregation, ratio to middle
(2) Variants that make a boundary: at each S20 change point (criterion A individual variants, criterion B variant character, PEN 1),
      per-variant carrier rate difference |fL-fR| in W=5 unit windows on both sides of the boundary (within the segment, >=3 units per side).
      Variants with >=0.6 = 'variants that split the boundary' -> class composition. Null = random positions in the same haplotype that are not change points.
      A variant in 1 unit has |fL-fR| <= 1/3 by definition and cannot be a 'splitting variant' -> rare (1 unit) variants are examined separately by density shift (2b).
    (2b) difference in count per unit by class between the two sides of the boundary |L-R|, observed vs null
(3) Finding boundaries from classes alone: the same binary segmentation (standardized within haplotype, penalty 1*d*logN, minimum segment 3)
      for inherited (carrier table), other shared and rare (1 unit) (counts of 12 characters) each, number of change points vs shuffle null, junction enrichment (observed/uniform expectation)
(4) Character by zone: best zones (proximal 1-3, 4-17, interior, distal 1-8) x 16 characters, multi-individual shared type vs single-unit type
      compute log2(zone/interior) per class and check whether the zone deviations of the two classes are alike (correlation).
Output: primary\S26_*.tsv, log
"""
import sys, io, os, collections, math
import numpy as np
import v11_io as G

PRESET = "primary"
OUT = os.path.join(G.OUT, PRESET)
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "S16_abs_vs_rel.py"), encoding="utf8").read()
exec(compile(src.split("STR=[(40,80")[0], "head", "exec"))   # units, ap, ad, Nn, hid, HK, U, parse, is_cpg, elem
rng = np.random.default_rng(20260916)

ORIGIN = {}
for i, L in enumerate(io.open(os.path.join(G.OUT, "export", "V11_ORIGIN_CALL.tsv"), encoding="utf8")):
    f = L.rstrip("\n").split("\t")
    if i == 0: H = {c: j for j, c in enumerate(f)}; continue
    ORIGIN[f[0]] = f[H["origin_call"]]

# ---------- classes ----------
hapcnt = collections.Counter()
for h in range(NH):
    s = set()
    for i in np.where(hid == h)[0]: s |= set(units[i])
    hapcnt.update(s)
arrcnt = [collections.Counter() for _ in range(NH)]
for i, vs in enumerate(units): arrcnt[hid[i]].update(vs)
# new 3-way split (definition of fork 4815ef09, reproduced with code/V2_apply.py)
#   cohort = haplotypes with >= 30 units, f = carrier fraction within that cohort, fT = 0.10
FT3 = 0.10; MINU3 = 30
_big = [h for h in range(NH) if int((hid == h).sum()) >= MINU3]
NS3 = len(_big)
hap3 = collections.Counter()
for h in _big:
    s = set()
    for i in np.where(hid == h)[0]: s |= set(units[i])
    hap3.update(s)
print("새 분류 코호트 : 하플로타입 %d / %d (유닛 >= %d) · fT %.2f" % (NS3, NH, MINU3, FT3))
CLS = ["다개체 공유형", "배열내 다중형", "단일유닛형"]
def klass(h, t):
    if hap3[t] / NS3 >= FT3: return "다개체 공유형"
    return "배열내 다중형" if arrcnt[h][t] >= 2 else "단일유닛형"
UC = []   # count of each class per unit
for i, vs in enumerate(units):
    c = collections.Counter(klass(hid[i], t) for t in vs)
    UC.append([c[k] for k in CLS])
UC = np.array(UC, float)
tot = UC.sum(0)
print("부류별 (하플로타입×유닛×변이) 건수 : " + " · ".join("%s %d (%.1f%%)" % (k, tot[j], 100 * tot[j] / tot.sum()) for j, k in enumerate(CLS)))
print("유닛당 평균 : " + " · ".join("%s %.3f" % (k, UC[:, j].mean()) for j, k in enumerate(CLS)))

# ---------- (1) position distribution ----------
print("\n=== (1) 부류별 유닛당 개수 — 가운데 대비 비 (2단계 · 하플로타입 부트스트랩 300) ===")
regs = [("P", d) for d in range(1, 31)] + [("MID", 0)] + [("D", d) for d in range(30, 0, -1)]
def regmask(side, d):
    if side == "P": return (ap == d)
    if side == "D": return (ad == d)
    return (ap > 30) & (ad > 30)
HM = {}
for side, d in regs:
    m = regmask(side, d)
    s = np.zeros((NH, len(CLS))); c = np.zeros(NH)
    np.add.at(s, hid[m], UC[m]); np.add.at(c, hid[m], 1)
    ok = c > 0
    HM[(side, d)] = s[ok] / c[ok][:, None]
mid = HM[("MID", 0)].mean(0)
g = io.open(os.path.join(OUT, "S26b_position_profile.tsv"), "w", encoding="utf8", newline="\n")
g.write("region\tclass\tn_hap\tper_unit\tratio_mid\tlo\thi\n")
for side, d in regs:
    Mx = HM[(side, d)]; n = len(Mx)
    bs = np.array([Mx[rng.integers(0, n, n)].mean(0) for _ in range(300)])
    for j, k in enumerate(CLS):
        v = Mx[:, j].mean()
        g.write("%s\t%s\t%d\t%.5f\t%.4f\t%.4f\t%.4f\n" % ("MID" if side == "MID" else "%s%d" % (side, d), k, n, v, v / mid[j],
                np.percentile(bs[:, j], 2.5) / mid[j], np.percentile(bs[:, j], 97.5) / mid[j]))
g.close()
print("  구역   " + "".join("%10s" % k for k in CLS))
for side, d in [("P", 1), ("P", 2), ("P", 3), ("P", 5), ("P", 10), ("P", 17), ("P", 30), ("MID", 0), ("D", 30), ("D", 17), ("D", 10), ("D", 8), ("D", 5), ("D", 3), ("D", 2), ("D", 1)]:
    v = HM[(side, d)].mean(0)
    print("  %-5s " % ("MID" if side == "MID" else "%s%d" % (side, d)) + "".join("%10.2f" % (v[j] / mid[j]) for j in range(len(CLS))))
print("  가운데 유닛당 : " + " · ".join("%s %.3f" % (k, mid[j]) for j, k in enumerate(CLS)))

# ---------- reading change points ----------
CP = collections.defaultdict(dict)
for i, L in enumerate(io.open(os.path.join(OUT, "S20_cp_unbiased.tsv"), encoding="utf8")):
    f = L.rstrip("\n").split("\t")
    if i == 0: continue
    if f[1] != "1.0": continue
    CP[f[0]][f[2]] = (int(f[3]), [int(x) for x in f[5].split(",") if x])
HKS = {"%s_%s" % k: h for h, k in enumerate(HK)}
ORD = {h: np.where(hid == h)[0][np.argsort(ap[np.where(hid == h)[0]])] for h in range(NH)}

# ---------- (2) variants that split a boundary ----------
W = 5; THR = 0.6
def window_stats(h, t, lo_b, hi_b):
    """t = boundary (sorted index t-1 | t). At most W units on each side, within the segment [lo_b, hi_b)"""
    idx = ORD[h]
    L = idx[max(lo_b, t - W):t]; R = idx[t:min(hi_b, t + W)]
    if len(L) < 3 or len(R) < 3: return None
    cL = collections.Counter(); cR = collections.Counter()
    for i in L: cL.update(units[i])
    for i in R: cR.update(units[i])
    comp = collections.Counter(); dens = UC[L].mean(0) - UC[R].mean(0)
    for tkn in set(cL) | set(cR):
        dd = abs(cL[tkn] / float(len(L)) - cR[tkn] / float(len(R)))
        if dd >= THR: comp[klass(h, tkn)] += 1
    return comp, np.abs(dens)
# split boundary positions into the 4 cohort zones (proximal 1-3, 4-17, middle, distal 1-10)
#   for boundary t, units on the proximal side = t and on the distal side = n - t. If it spans two zones, take the shallower one (same rule as VTZONE.py zone_of). Every boundary falls in one bin
#   (old: proximal junction 1-6 (t <= 6), distal junction 4-10 (4 <= n-t <= 10), interior (>=15 on both sides), the remaining 'in between' was not in the figure)
def loc_of(t, n):
    dp, dd = t, n - t
    pz = "근위 1-3" if dp <= 3 else ("근위 4-17" if dp <= 17 else None)
    dz = "원위 1-10" if dd <= 10 else None
    if pz and dz: return pz if dp <= dd else dz
    return pz or dz or "가운데"
g = io.open(os.path.join(OUT, "S26b_boundary_drivers.tsv"), "w", encoding="utf8", newline="\n")
g.write("basis\tset\tlocation\tn_boundaries\t" + "\t".join("n_%s" % k for k in CLS) + "\t" + "\t".join("dens_diff_%s" % k for k in CLS) + "\n")
print("\n=== (2) 경계를 가르는 변이의 기원 (|보유율 차| ≥%.1f · 창 %d 유닛) ===" % (THR, W))
for basis in ("A", "B"):
    RES = collections.defaultdict(lambda: [collections.Counter(), [], 0])
    for hk, (n, cps) in CP[basis].items():
        h = HKS[hk]
        bounds = [0] + cps + [n]
        for j, t in enumerate(cps):
            r = window_stats(h, t, bounds[j], bounds[j + 2])
            if r is None: continue
            for key in (("실제", loc_of(t, n)), ("실제", "전체")):
                RES[key][0].update(r[0]); RES[key][1].append(r[1]); RES[key][2] += 1
        # null: random positions at least ±4 away from change points, same number
        avoid = set()
        for t in cps:
            for q in range(t - 4, t + 5): avoid.add(q)
        cand = [q for q in range(3, n - 2) if q not in avoid]
        for t in (rng.choice(cand, min(len(cps), len(cand)), replace=False) if cand and cps else []):
            r = window_stats(h, int(t), 0, n)
            if r is None: continue
            for key in (("귀무", loc_of(int(t), n)), ("귀무", "전체")):
                RES[key][0].update(r[0]); RES[key][1].append(r[1]); RES[key][2] += 1
    print("\n  [기준 %s · %s]" % (basis, "개별 변이 160 종" if basis == "A" else "변이 성격 20 개"))
    print("   %-4s %-16s %6s | 경계당 가르는 변이 수 : %s | 경계 양쪽 유닛당 개수 차 : %s" % ("", "위치", "경계", " ".join("%8s" % k for k in CLS), " ".join("%8s" % k for k in CLS)))
    for loc in ("전체", "근위 1-3", "근위 4-17", "가운데", "원위 1-10"):
        for st in ("실제", "귀무"):
            c, dl, nb = RES[(st, loc)]
            if nb == 0: continue
            dm = np.mean(dl, 0)
            g.write("%s\t%s\t%s\t%d\t%s\t%s\n" % (basis, st, loc, nb, "\t".join("%.4f" % (c[k] / float(nb)) for k in CLS), "\t".join("%.4f" % x for x in dm)))
            print("   %-4s %-16s %6d | %s | %s" % (st, loc, nb, " ".join("%8.2f" % (c[k] / float(nb)) for k in CLS), " ".join("%8.3f" % x for x in dm)))
g.close()

# ---------- (3) finding boundaries from classes alone ----------
def best_split(S, min_size=3):
    m = len(S)
    if m < 2 * min_size: return None
    cs1 = np.cumsum(S, 0); cs2 = np.cumsum(S * S, 0); tot1 = cs1[-1]; tot2 = cs2[-1]
    mu = tot1 / m; c0 = float((tot2 - 2 * mu * tot1 + m * mu * mu).sum())
    best = None
    for t in range(min_size, m - min_size + 1):
        l1 = cs1[t - 1]; l2 = cs2[t - 1]
        cl = float((l2 - l1 * l1 / t).sum()); r1 = tot1 - l1; r2 = tot2 - l2; rr = m - t
        cr = float((r2 - r1 * r1 / rr).sum())
        if best is None or c0 - (cl + cr) > best[0]: best = (c0 - (cl + cr), t)
    return best
def binseg(X, pen, min_size=3):
    n = len(X); segs = [(0, n)]; cps = []
    while True:
        cand = None
        for si, (a, b) in enumerate(segs):
            r = best_split(X[a:b], min_size)
            if r is None: continue
            if cand is None or r[0] > cand[0]: cand = (r[0], si, r[1])
        if cand is None or cand[0] <= pen: break
        _, si, t = cand
        a, b = segs.pop(si); segs.append((a, a + t)); segs.append((a + t, b)); segs.sort(); cps.append(a + t)
    return sorted(cps)
CHAR = ["C>T", "C>A", "C>G", "T>C", "T>A", "T>G", "CpG", "AluY", "NTS1", "GArich", "NTS2", "coding", "del", "ins"]
def charvec(tkns):
    v = np.zeros(len(CHAR))
    for t in tkns:
        cl, p, r, a, ln, rep = parse(t)
        if cl == "snv":
            key = {("C", "T"): "C>T", ("G", "A"): "C>T", ("C", "A"): "C>A", ("G", "T"): "C>A", ("C", "G"): "C>G", ("G", "C"): "C>G",
                   ("T", "C"): "T>C", ("A", "G"): "T>C", ("T", "A"): "T>A", ("A", "T"): "T>A", ("T", "G"): "T>G", ("A", "C"): "T>G"}.get((r, a))
            if key: v[CHAR.index(key)] += 1
            if (r, a) in (("C", "T"), ("G", "A")) and is_cpg(p, r): v[CHAR.index("CpG")] += 1
        else:
            v[CHAR.index("del" if cl == "del" else "ins")] += 1
        e = elem(p)
        if e: v[CHAR.index(e)] += 1
    return v
def feats(h, cls):
    idx = ORD[h]; n = len(idx)
    if cls in ("다개체 공유형", "배열내 다중형"):
        tk = sorted(t for t in arrcnt[h] if klass(h, t) == cls and 0.05 * n <= arrcnt[h][t] <= 0.95 * n)
        if len(tk) < 3: return None
        TI = {t: j for j, t in enumerate(tk)}
        X = np.zeros((n, len(tk)))
        for r_, i in enumerate(idx):
            for t in units[i]:
                j = TI.get(t)
                if j is not None: X[r_, j] = 1
    else:
        X = np.array([charvec([t for t in units[i] if klass(h, t) == "단일유닛형"]) for i in idx])
    X = X - X.mean(0); sd = X.std(0); keep = sd > 0
    if keep.sum() < 2: return None
    X = X[:, keep] / sd[keep]
    if X.shape[1] > 5:
        u, s_, vt = np.linalg.svd(X, full_matrices=False); X = u[:, :5] * s_[:5]; X = X / np.maximum(X.std(0), 1e-9)
    return X
print("\n=== (3) 부류만으로 배열마다 경계 찾기 (벌점 1·d·logN · 순서섞기 귀무) ===")
bins = [(1, 3), (4, 6), (7, 10), (11, 15), (16, 20), (21, 30), (31, 60)]
g = io.open(os.path.join(OUT, "S26b_class_segmentation.tsv"), "w", encoding="utf8", newline="\n")
g.write("class\tn_hap\tmean_cp\tmean_cp_null\tside\tdist_bin\tobs\texp\tobs_over_exp\n")
for cls in ("다개체 공유형", "배열내 다중형", "단일유닛형"):
    ncp = []; nnull = []; so = {s: np.zeros(400) for s in "PD"}; se = {s: np.zeros(400) for s in "PD"}
    for h in range(NH):
        n = len(ORD[h])
        if n < 30: continue
        X = feats(h, cls)
        if X is None: continue
        pen = 1.0 * X.shape[1] * np.log(n)
        c = binseg(X, pen); ncp.append(len(c)); nnull.append(len(binseg(X[rng.permutation(n)], pen)))
        m = n - 5
        for x in c:
            if x <= n - x: so["P"][x] += 1
            else: so["D"][n - x] += 1
        if c:
            for x in range(3, n - 2):
                w = len(c) / float(m)
                if x <= n - x: se["P"][x] += w
                else: se["D"][n - x] += w
    print("  [%s] 하플로타입 %d · 변화점 평균 %.2f (순서섞기 귀무 %.2f)" % (cls, len(ncp), np.mean(ncp), np.mean(nnull)))
    for s, lab in (("P", "근위"), ("D", "원위")):
        row = []
        for a_, b_ in bins:
            o = so[s][a_:b_ + 1].sum(); e = se[s][a_:b_ + 1].sum()
            row.append("%d-%d %.2f" % (a_, b_, o / e if e > 0 else float("nan")))
            g.write("%s\t%d\t%.3f\t%.3f\t%s\t%d-%d\t%d\t%.2f\t%.3f\n" % (cls, len(ncp), np.mean(ncp), np.mean(nnull), lab, a_, b_, o, e, o / e if e > 0 else float("nan")))
        print("     %s 관측/기대 : %s" % (lab, " · ".join(row)))
g.close()

# ---------- (4) character by zone: multi-individual shared type vs single-unit type ----------
print("\n=== (4) 최적 구역별 성격 — 다개체 공유형 vs 단일유닛형 (유닛 1000 개당 · log2 구역/내부) ===")
def zone_of(i):
    if ad[i] <= 8: return "원위1-8"
    if ap[i] <= 3: return "근위1-3"
    if ap[i] <= 17: return "근위4-17"
    return "내부"
ZN = ["근위1-3", "근위4-17", "내부", "원위1-8"]
zn = np.array([zone_of(i) for i in range(NU)])
g = io.open(os.path.join(OUT, "S26b_zone_character.tsv"), "w", encoding="utf8", newline="\n")
g.write("class\tzone\tn_units\tfeature\tper_1000\tlog2_vs_interior\n")
DEV = {}
for cls, members in (("다개체 공유형", ("다개체 공유형",)), ("단일유닛형", ("단일유닛형",)), ("배열내 다중형", ("배열내 다중형",))):
    Z = {}
    for z in ZN:
        m = np.where(zn == z)[0]
        v = np.zeros(len(CHAR))
        for i in m: v += charvec([t for t in units[i] if klass(hid[i], t) in members])
        Z[z] = (v / len(m) * 1000, len(m))
    DEV[cls] = {}
    for z in ZN:
        dv = np.log2((Z[z][0] + 0.5) / (Z["내부"][0] + 0.5))
        DEV[cls][z] = dv
        for j, fname in enumerate(CHAR):
            g.write("%s\t%s\t%d\t%s\t%.2f\t%.3f\n" % (cls, z, Z[z][1], fname, Z[z][0][j], dv[j]))
    print("  [%s] 유닛 1000 개당 총 변이 : %s" % (cls, " · ".join("%s %.0f" % (z, Z[z][0][CHAR.index("AluY")] + Z[z][0][CHAR.index("NTS1")] + Z[z][0][CHAR.index("GArich")] + Z[z][0][CHAR.index("NTS2")] + Z[z][0][CHAR.index("coding")]) for z in ZN)))
    for z in ("근위1-3", "근위4-17", "원위1-8"):
        print("     %-8s log2(구역/내부) : %s" % (z, " ".join("%s %+.2f" % (CHAR[j], DEV[cls][z][j]) for j in range(len(CHAR)))))
g.close()
print("\n  계승형 편차와 드문(1유닛) 편차의 상관 (성격 14 개 · 구역별)")
for z in ("근위1-3", "근위4-17", "원위1-8"):
    a_ = DEV["다개체 공유형"][z]; b_ = DEV["단일유닛형"][z]
    r = np.corrcoef(a_, b_)[0, 1]
    rr = np.array([np.corrcoef(a_, rng.permutation(b_))[0, 1] for _ in range(2000)])
    print("   %-8s r = %+.2f (순열 p = %.3f)" % (z, r, (np.sum(np.abs(rr) >= abs(r)) + 1) / 2001.0))
print("\n저장 : S26_position_profile.tsv · S26_boundary_drivers.tsv · S26_class_segmentation.tsv · S26_zone_character.tsv")

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('S26b_boundary_drivers.tsv', 'S26b_position_profile.tsv'): _taben.to_en_file(_o.path.join(OUT, _f))
