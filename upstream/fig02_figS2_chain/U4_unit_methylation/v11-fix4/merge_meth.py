# -*- coding: utf-8 -*-
# methylation merge (2026-09-12), collects read x unit methylation from the 3 nodes into per-position tables
#   input: METH_READUNIT_<node>.tsv.gz (MM/ML of the placement BAM transferred to canonical coordinates)
#   duplicate rule: if the same read is in both junction and ARRAY_ONLY (DUAL), the junction one is kept
#                 (same rule as the merged map, the ARRAY_ONLY DUAL entries are removed)
#   source separation: even for the same sample, different data sources are written separately (same rule as the variant merge)
#   outputs
#     READMAP_METH_G8V4.tsv  : per sample, source, slot and unit number, read count, CpG call count, beta (call-based),
#                              mean/median read beta, junction/ARRAY_ONLY read counts
#     METH_SITE_G8V4.tsv.gz  : the above plus n, n_meth and beta per within-unit coordinate (1-2240)
#   Arguments: input folder, output folder
import io
import os
import sys
import gzip
import glob
import statistics
import collections

NL = chr(10)
TAB = chr(9)
IN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)

seen = {}                       # (smp, src, slot, ap, read) -> kind
RU = collections.defaultdict(list)      # (smp, src, slot, ap) -> [(beta, ncpg, nmeth, kind)]
SITE = collections.Counter()
SITEM = collections.Counter()
T = collections.Counter()

files = sorted(glob.glob(os.path.join(IN, "METH_READUNIT_*.tsv.gz")))
print("입력 파일 %d" % len(files))
for p in files:
    f = gzip.open(p, "rt", encoding="utf8")
    ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    for L in f:
        r = L.rstrip(NL).split(TAB)
        k = (r[0], r[ix["data_source"]], r[ix["slot"]], int(r[ix["array_pos"]]))
        rid = r[ix["read_id"]]
        kind = r[ix["kind"]]
        if kind == "junction" and r[0] in ("HG00384", "HG01150"):
            T["junction 제외 개체(HG00384·HG01150)"] += 1      # same rule as the variant merge (self-reference of the assignment table)
            continue
        kk = k + (rid,)
        if kk in seen:
            T["중복(DUAL 등) 버림"] += 1
            if seen[kk] == "array_only" and kind == "junction":
                pass                     # if junction comes later it is not replaced below (already counted)
            continue
        seen[kk] = kind
        n = int(r[ix["n_cpg"]])
        nm = int(r[ix["n_meth"]])
        RU[k].append((nm / float(n) if n else 0.0, n, nm, kind))
        T["read-유닛"] += 1
        for c in r[ix["calls"]].split(","):
            if not c:
                continue
            a, v = c.split(":")
            key = k + (int(a),)
            SITE[key] += 1
            if int(v) >= 128:
                SITEM[key] += 1
    f.close()

o = io.open(os.path.join(OUT, "READMAP_METH_G8V4.tsv"), "w", encoding="utf8", newline=NL)
o.write(TAB.join(["sample", "data_source", "slot", "array_pos", "n_reads", "n_calls", "n_meth",
                  "beta_call", "beta_read_mean", "beta_read_median", "reads_jx", "reads_ao"]) + NL)
for k in sorted(RU):
    v = RU[k]
    nc = sum(x[1] for x in v)
    nm = sum(x[2] for x in v)
    bs = [x[0] for x in v]
    o.write(TAB.join([k[0], k[1], k[2], str(k[3]), str(len(v)), str(nc), str(nm),
                      "%.4f" % (nm / float(nc) if nc else 0),
                      "%.4f" % statistics.mean(bs),
                      "%.4f" % statistics.median(bs),
                      str(sum(1 for x in v if x[3] == "junction")),
                      str(sum(1 for x in v if x[3] == "array_only"))]) + NL)
o.close()

g = gzip.open(os.path.join(OUT, "METH_SITE_G8V4.tsv.gz"), "wt", encoding="utf8", newline=NL)
g.write(TAB.join(["sample", "data_source", "slot", "array_pos", "unit_pos", "n", "n_meth", "beta"]) + NL)
for k in sorted(SITE):
    n = SITE[k]
    nm = SITEM.get(k, 0)
    g.write(TAB.join([k[0], k[1], k[2], str(k[3]), str(k[4]), str(n), str(nm),
                      "%.4f" % (nm / float(n))]) + NL)
g.close()

print("자리(유닛) %d · 자리(염기) %d" % (len(RU), len(SITE)))
for k, v in sorted(T.items()):
    print("  %-20s %d" % (k, v))
bs = [x[0] for v in RU.values() for x in v]
if bs:
    lo = sum(1 for b in bs if b < 0.2)
    hi = sum(1 for b in bs if b > 0.7)
    print("  read-유닛 beta : 저메틸<0.2 %d (%.1f%%) · 중간 %d · 고메틸>0.7 %d (%.1f%%)"
          % (lo, 100.0 * lo / len(bs), len(bs) - lo - hi, hi, 100.0 * hi / len(bs)))
