# -*- coding: utf-8 -*-
# Stage table for all individuals, classification of 313 individuals, recomputed on G8v4 (2026-09-12, rules of stage_all3.py unchanged)
#   changed inputs: map CONSENSUS_MAP_G6C2 -> CONSENSUS_MAP_G8v4_full, merge READMAP_MERGED_v7 -> READMAP_G8V4_BYSOURCE
#   rules (same as the old version): cov1 = fraction of units with 1 or more reads, cov3 = fraction with depth 3 or more
#                       agreement = TP/(TP+FN+FP) of units with depth 3 or more, decision thresholds cov3 >= 90% and agreement >= 95%
#   The new merge has rows split by source (19 haplotypes with two sources). So two versions are produced together.
#     union : sources merged per position (sample, slot, array_pos), the main version comparable with the old table (source-merged)
#     rep   : only the source with the largest summed depth in that haplotype, the version that follows the user rule (do not mix sources)
#   The canonical source (canon_src) and slot list are taken from CONSENSUS_SUM_G6.tsv unchanged, since the canonical did not change.
#   Output is written to the work folder, not the canonical path.
#   Arguments: SUM_G6, MAP_G8v4, READMAP_G8V4, HAP_READ_NEW44, HAP_READ_FINAL, output folder
import io
import os
import sys
import collections

NL = chr(10)
TAB = chr(9)
SUMF, MAPF, RMF, NEW44F, HRSF, OUT = sys.argv[1:7]
os.makedirs(OUT, exist_ok=True)
PUB = set(["asmunits3", "hgsvcunits"])
SPAN = set(["tier1ref"])
C2 = set([("HG02683", "1"), ("NA19090", "1"), ("HG03126", "1"),
          ("HG03366", "1"), ("NA19120", "1"), ("HG03198", "1")])

HAP = {}
f = io.open(SUMF, encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    d = r[2].split("__")[0]
    HAP[(r[0], r[1])] = ("cons2" if (r[0], r[1]) in C2 else
                         "pub" if d in PUB else "span" if d in SPAN else "hifi")
SMP = sorted(set(k[0] for k in HAP))
NSL = collections.Counter(k[0] for k in HAP)
NEW44 = set(L.split(TAB)[0] for i, L in enumerate(io.open(NEW44F, encoding="utf8")) if i)
HRS = set(L.split(TAB)[0] for i, L in enumerate(io.open(HRSF, encoding="utf8")) if i)


def cohort(s):
    if s in ("hg001", "hg002", "hg003", "hg004"):
        return "GIAB"
    if s in ("COLO829-BL", "COLO829-T", "HG008-T"):
        return "암"
    return "new44" if s in NEW44 else "HPRC"


# map: unit count per haplotype (G8v4, seg_type UNIT)
NU = collections.Counter()
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] == "UNIT":
        NU[(r[0], r[1])] += 1

# merge: per position and per source
f = io.open(RMF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
SITE = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0]))
# SITE[(smp, slot)][(src, pos)] = [depth, TP, FN, FP, jx, ao]
for L in f:
    r = L.rstrip(NL).split(TAB)
    k = (r[0], r[ix["slot"]])
    src = r[ix["data_source"]]
    pos = int(r[ix["array_pos"]])
    d = int(r[ix["depth"]])
    e = SITE[k][(src, pos)]
    e[0] += d
    e[4] += int(r[ix["depth_jx"]]) > 0
    e[5] += int(r[ix["depth_ao"]]) > 0
    if d >= 3 and r[ix["agree"]] != "NOMAP":
        e[1] += int(r[ix["n_TP"]])
        e[2] += int(r[ix["n_FN"]])
        e[3] += int(r[ix["n_FP"]])


def tally(k, keys):
    """keys = list of (src,pos) -> merged per position into cov1 cov3 TP FN FP jx ao"""
    P = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0])
    for kk in keys:
        e = SITE[k][kk]
        p = P[kk[1]]
        for i, v in enumerate((e[0], e[1], e[2], e[3], e[4], e[5])):
            p[i] += v
    c1 = c3 = tp = fn = fp = jx = ao = 0
    for pos, p in P.items():
        c1 += p[0] >= 1
        c3 += p[0] >= 3
        tp += p[1]
        fn += p[2]
        fp += p[3]
        jx += p[4] > 0
        ao += p[5] > 0
    return c1, c3, tp, fn, fp, jx, ao


def stage(kind, jx, ao, cov3, idt):
    if kind == "span":
        return "B_정본무효(관통read)"
    if jx == 0 and ao == 0:
        return "C_배치산출0"
    if cov3 >= 90 and idt >= 95:
        return "OK"
    if cov3 < 90:
        return "D_피복부족"
    return "E_정본불일치"


rows = []
for k in sorted(HAP):
    smp = k[0]
    nu = NU.get(k, 0)
    allk = list(SITE.get(k, {}).keys())
    # rep: source with the largest summed depth
    bysrc = collections.Counter()
    for (src, pos) in allk:
        bysrc[src] += SITE[k][(src, pos)][0]
    rep = bysrc.most_common(1)[0][0] if bysrc else "-"
    res = {}
    for tag, keys in (("u", allk), ("r", [x for x in allk if x[0] == rep])):
        c1, c3, tp, fn, fp, jx, ao = tally(k, keys)
        cov1 = 100.0 * c1 / nu if nu else 0.0
        cov3 = 100.0 * c3 / nu if nu else 0.0
        idt = 100.0 * tp / (tp + fn + fp) if tp + fn + fp else 0.0
        res[tag] = (cov1, cov3, idt, fn, fp, jx, ao, stage(HAP[k], jx, ao, cov3, idt))
    u, rp = res["u"], res["r"]
    ok1 = HAP[k] != "span" and u[0] >= 90 and u[2] >= 95
    rows.append([k, cohort(smp), HAP[k], nu, u[5], u[6], u[0], u[1], u[2], u[3], u[4], u[7],
                 ok1, smp in HRS, len(set(x[0] for x in allk)), rep, rp[1], rp[2], rp[7]])

o = io.open(os.path.join(OUT, "STAGE_ALL_G8V4.tsv"), "w", encoding="utf8", newline=NL)
o.write(TAB.join(["sample", "slot", "cohort", "canon_src", "n_unit", "jx_units", "ao_units",
                  "cov1%", "cov3%", "var_ident%", "FN", "FP", "stage", "ok_cov1def",
                  "hapread_labeled", "n_source", "rep_source", "cov3_rep%", "var_ident_rep%",
                  "stage_rep"]) + NL)
for x in rows:
    o.write(TAB.join([x[0][0], x[0][1], x[1], x[2], str(x[3]), str(x[4]), str(x[5]),
                      "%.1f" % x[6], "%.1f" % x[7], "%.2f" % x[8], str(x[9]), str(x[10]), x[11],
                      "Y" if x[12] else "N", "Y" if x[13] else "N", str(x[14]), x[15],
                      "%.1f" % x[16], "%.2f" % x[17], x[18]]) + NL)
o.close()

CO = ["HPRC", "new44", "GIAB", "암"]
print("개체 %d · 하플로타입 %d · 지도 유닛 합 %d" % (len(SMP), len(rows), sum(x[3] for x in rows)))
print("== 하플로타입 판정 (union · 깊이3 정의)")
print("  %-22s %s    계" % ("", " ".join("%6s" % c for c in CO)))
for s in sorted(set(x[11] for x in rows)):
    c = collections.Counter(x[1] for x in rows if x[11] == s)
    print("  %-22s %s  %4d" % (s, " ".join("%6d" % c[q] for q in CO), sum(c.values())))
print("== 대표출처판(rep) 판정")
for s in sorted(set(x[18] for x in rows)):
    c = collections.Counter(x[1] for x in rows if x[18] == s)
    print("  %-22s %s  %4d" % (s, " ".join("%6d" % c[q] for q in CO), sum(c.values())))
d = [x for x in rows if x[14] > 1]
print("== 출처 2개인 하플로타입 %d : union 과 rep 판정이 다른 것 %d"
      % (len(d), sum(1 for x in d if x[11] != x[18])))
for x in d:
    if x[11] != x[18]:
        print("   %s_%s : union %s (피복 %.1f · 일치 %.2f) / rep %s %s (피복 %.1f · 일치 %.2f)"
              % (x[0][0], x[0][1], x[11], x[7], x[8], x[15], x[18], x[16], x[17]))

# -- CLASS313 (on union, only the decision part of class313.py, unchanged)
ST = dict(((x[0][0], x[0][1]), x) for x in rows)
cls = []
for s in SMP:
    sl = sorted(k[1] for k in ST if k[0] == s)
    ok = [ST[(s, x)][11] == "OK" for x in sl]
    ok1 = [ST[(s, x)][12] for x in sl]

    def cat(v):
        n = sum(1 for z in v if z)
        return "양쪽 신뢰" if n == 2 else ("한쪽만 신뢰" if n == 1 else "둘다 신뢰불가")
    stat = " / ".join("슬롯%s %s(피복%.1f%% 일치%.2f%%)"
                      % (x, ST[(s, x)][11], ST[(s, x)][7], ST[(s, x)][8]) for x in sl)
    cls.append([s, cohort(s), str(len(sl)), cat(ok), cat(ok1), stat])
o = io.open(os.path.join(OUT, "CLASS313_G8V4.tsv"), "w", encoding="utf8", newline=NL)
o.write(TAB.join(["sample", "cohort", "n_canon_slots", "reliability_depth3",
                  "reliability_depth1", "status_per_slot"]) + NL)
for r in cls:
    o.write(TAB.join(r) + NL)
o.close()
CT = ["양쪽 신뢰", "한쪽만 신뢰", "둘다 신뢰불가"]
print("== 개체 분류 (union · 깊이3)")
for c in CO:
    rs = [r for r in cls if r[1] == c]
    if not rs:
        continue
    a = collections.Counter(r[3] for r in rs)
    print("  %-8s %3d 개체 : %s" % (c, len(rs), " · ".join("%s %d" % (t, a[t]) for t in CT)))
a = collections.Counter(r[3] for r in cls)
print("  %-8s %3d 개체 : %s" % ("합계", len(cls), " · ".join("%s %d" % (t, a[t]) for t in CT)))
print("  정본 슬롯 1 개 개체 %d" % sum(1 for r in cls if r[2] == "1"))
