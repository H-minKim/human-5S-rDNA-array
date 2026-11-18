# -*- coding: utf-8 -*-
# Re-measure the bonobo 5S array methylation per basecaller release, to settle the default value.
#   Input: the bonobo read_unit table from 03_reads, holding read x unit with n_calls and n_5mC, plus the read-to-release mapping from the read ids in the flowcell selection SAMs.
#   Values: the call-weighted 5mC fraction, sum(n_5mC) / sum(n_calls), the median over reads, and the read count.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import os
import collections
import csv
import io
import statistics
import sys

RU = os.environ["RDNA_OTHERSP_DIR"] + "/data_v2/03_reads/bonobo_mPanPan1/bonobo_mPanPan1__bonobo_mPanPan1.read_unit.tsv"
MAP = sys.argv[1]
OUT = os.environ["RDNA_OTHERSP_DIR"] + "/data_v2/00_catalog/CMD01_BONOBO_BY_BASECALLER.tsv"
bc = {}
for L in io.open(MAP, encoding="utf8"):
    a, b = L.rstrip("\n").split("\t")
    bc[a] = b
per_read = collections.defaultdict(lambda: [0, 0])
for r in csv.DictReader(io.open(RU, encoding="utf8"), delimiter="\t"):
    if r["context"] != "CG":
        continue
    x = per_read[r["read_id"]]
    x[0] += int(r["n_calls"])
    x[1] += int(r["n_5mC"])
g = collections.defaultdict(list)
miss = 0
for rid, (c, m) in per_read.items():
    v = bc.get(rid)
    if v is None:
        miss += 1
        v = "unknown"
    g[v].append((c, m))
g["all(mixed)"] = [x for k in list(g) for x in g[k]]
o = io.open(OUT, "w", encoding="utf8", newline="\n")
o.write("basecaller\tn_read\tcalls\tcall_weighted_frac\tread_median_frac\n")
for k in sorted(g):
    v = [x for x in g[k] if x[0] > 0]
    c = sum(x[0] for x in v)
    m = sum(x[1] for x in v)
    med = statistics.median(x[1] / x[0] for x in v)
    o.write("%s\t%d\t%d\t%.3f\t%.3f\n" % (k, len(v), c, m / c, med))
    print("%-12s read %4d, calls %8d, call-weighted %.3f, read median %.3f" % (k, len(v), c, m / c, med))
o.close()
print("reads of unknown version", miss)
