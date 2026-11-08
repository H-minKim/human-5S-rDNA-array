# -*- coding: utf-8 -*-
"""
stage_v11 - recompute the haplotype stage table on the v11 basis.
  The rules were verified against `01_code\\reference_G8v4_code\\stage_all4.py` and carried over **as is**:
    cov1 = fraction of units with 1 or more reads, cov3 = fraction with depth 3 or more
    agreement = TP/(TP+FN+FP) over units with depth 3 or more
    stage : B_정본무효 if canon_src is span(tier1ref), C_배치산출0 if jx=ao=0
            cov3>=90 & agreement>=95 -> OK, cov3<90 -> D_피복부족, otherwise E_정본불일치
    two versions: union (sources merged) and rep (the single source with the largest depth sum).
  Differences from the original (recorded from measurement):
    - canon_src is read from v11 `02_canonical\\hap_seq_node_v8.tsv` (no header) instead of SUM_G6.
      prefix distribution: asmunits3 400, arr 142, fill60 26, tier1ref 12, resub2 10, refix 8,
                  cons2seq 6, sweep31 5, nopurge 5, resub 1, hg002asm 1 (total 616)
    - cohort is not read from the old table, the rule of stage_all4.py is carried over into code as is
      (hg001-4=GIAB, COLO829-BL/T and HG008-T=cancer, the rest=HPRC).
      The NEW44 list is not in the release, so the new44 split cannot be produced (it merges into HPRC).
      new44 is only a reporting split and is not used in stage calls or filters. The hapread_labeled column is also dropped.
  Output is written to the work folder, not to the canonical folder.
"""
import io, os, sys, collections

V = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11-fix2")   # v11-fix2: the map is from fix1, the merged table from fix2
V1 = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11-fix1")
V0 = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11")
OUT = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\5S_position_variant_work\2026-09-17_v11_fix2")
os.makedirs(OUT, exist_ok=True)
MAPF = os.path.join(V1, "03_maps", "g8canonv11f1", "CONSENSUS_MAP_G8_full.tsv")
RMF = os.path.join(V, "07_merged", "READMAP_V11F2_BYSOURCE.tsv")
HSN = os.path.join(V0, "02_canonical", "hap_seq_node_v8.tsv")
def cohort(s):
    """Carries over the rules of stage_all4.py as is (does not read the old table)"""
    if s in ("hg001", "hg002", "hg003", "hg004"): return "GIAB"
    if s in ("COLO829-BL", "COLO829-T", "HG008-T"): return "암"
    return "HPRC"

PUB = set(["asmunits3", "hgsvcunits"])
SPAN = set(["tier1ref"])
C2 = set([("HG02683", "1"), ("NA19090", "1"), ("HG03126", "1"),
          ("HG03366", "1"), ("NA19120", "1"), ("HG03198", "1")])

HAP = {}
for L in io.open(HSN, encoding="utf8"):
    r = L.rstrip("\n").split("\t")
    if len(r) < 3:
        continue
    d = r[2].split("__")[0]
    HAP[(r[0], r[1])] = ("cons2" if (r[0], r[1]) in C2 else
                         "pub" if d in PUB else "span" if d in SPAN else "hifi")
print("정본 슬롯 %d · canon_src %s" % (len(HAP), dict(collections.Counter(HAP.values()))))

NU = collections.Counter()
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip("\n").split("\t")))
for L in f:
    r = L.rstrip("\n").split("\t")
    if r[ix["seg_type"]] == "UNIT":
        NU[(r[0], r[1])] += 1

f = io.open(RMF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip("\n").split("\t")))
SITE = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0]))
for L in f:
    r = L.rstrip("\n").split("\t")
    k = (r[0], r[ix["slot"]])
    e = SITE[k][(r[ix["data_source"]], int(r[ix["array_pos"]]))]
    d = int(r[ix["depth"]])
    e[0] += d
    e[4] += int(r[ix["depth_jx"]]) > 0
    e[5] += int(r[ix["depth_ao"]]) > 0
    if d >= 3 and r[ix["agree"]] != "NOMAP":
        e[1] += int(r[ix["n_TP"]]); e[2] += int(r[ix["n_FN"]]); e[3] += int(r[ix["n_FP"]])

def tally(k, keys):
    P = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0])
    for kk in keys:
        e = SITE[k][kk]; p = P[kk[1]]
        for i, v in enumerate(e):
            p[i] += v
    c1 = c3 = tp = fn = fp = jx = ao = 0
    for pos, p in P.items():
        c1 += p[0] >= 1; c3 += p[0] >= 3
        tp += p[1]; fn += p[2]; fp += p[3]
        jx += p[4] > 0; ao += p[5] > 0
    return c1, c3, tp, fn, fp, jx, ao

def stage(kind, jx, ao, cov3, idt):
    if kind == "span": return "B_정본무효(관통read)"
    if jx == 0 and ao == 0: return "C_배치산출0"
    if cov3 >= 90 and idt >= 95: return "OK"
    if cov3 < 90: return "D_피복부족"
    return "E_정본불일치"

rows = []
for k in sorted(HAP):
    nu = NU.get(k, 0)
    allk = list(SITE.get(k, {}).keys())
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
    rows.append([k, cohort(k[0]), HAP[k], nu, u[5], u[6], u[0], u[1], u[2], u[3], u[4], u[7],
                 len(set(x[0] for x in allk)), rep, rp[1], rp[2], rp[7]])

p = os.path.join(OUT, "STAGE_ALL_V11F2.tsv")
o = io.open(p, "w", encoding="utf8", newline="\n")
o.write("\t".join(["sample", "slot", "cohort", "canon_src", "n_unit", "jx_units", "ao_units",
                   "cov1%", "cov3%", "var_ident%", "FN", "FP", "stage",
                   "n_source", "rep_source", "cov3_rep%", "var_ident_rep%", "stage_rep"]) + "\n")
for x in rows:
    o.write("\t".join([x[0][0], x[0][1], x[1], x[2], str(x[3]), str(x[4]), str(x[5]),
                       "%.1f" % x[6], "%.1f" % x[7], "%.2f" % x[8], str(x[9]), str(x[10]), x[11],
                       str(x[12]), x[13], "%.1f" % x[14], "%.2f" % x[15], x[16]]) + "\n")
o.close()

CO = ["HPRC", "GIAB", "암"]
print("개체 %d · 하플로타입 %d · 지도 유닛 합 %d" % (len(set(x[0][0] for x in rows)), len(rows), sum(x[3] for x in rows)))
print("== union 판정")
for s in sorted(set(x[11] for x in rows)):
    c = collections.Counter(x[1] for x in rows if x[11] == s)
    print("  %-22s %s  계 %4d" % (s, " ".join("%6d" % c[q] for q in CO), sum(c.values())))
print("== rep 판정")
for s in sorted(set(x[16] for x in rows)):
    print("  %-22s %4d" % (s, sum(1 for x in rows if x[16] == s)))
d = [x for x in rows if x[12] > 1]
print("출처 2 개인 하플로타입 %d · union≠rep %d" % (len(d), sum(1 for x in d if x[11] != x[16])))
print("저장: %s" % p)
