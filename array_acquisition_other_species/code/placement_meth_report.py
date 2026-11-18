# -*- coding: utf-8 -*-
# Read placement and a methylation summary across every species, for the other species only.
#   Placement: counted from the array summaries in 03_reads and the fragment tables in 06_units_from_reads.
#   Methylation: the call count and fraction for CpG, or CG, CHG and CHH in the plant, from those same array summaries.
#   Outputs: PLACEMENT_METH_REPORT.tsv, by species, sample and context, with a per-species summary printed to the terminal.
#   Argument: the data_v2 folder.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import csv
import glob
import collections

TAB = chr(9)
NL = chr(10)
B = sys.argv[1]
AR = []
for f in glob.glob(os.path.join(B, "03_reads", "*", "5S_array_summary.tsv")):
    AR += list(csv.DictReader(io.open(f, encoding="utf8"), delimiter=TAB))
# the placement depth, from the read fragment table
DEP = {}
for f in glob.glob(os.path.join(B, "06_units_from_reads", "*.read_unit_seq.tsv")):
    base = os.path.basename(f)[:-len(".read_unit_seq.tsv")]
    sp, sam = base.split("__", 1)
    n = span = 0
    reads = set()
    ua = collections.Counter()
    for r in csv.DictReader(io.open(f, encoding="utf8"), delimiter=TAB):
        n += 1
        reads.add(r["read_id"])
        if r["spans_unit"] == "Y":
            span += 1
        ua[(r["array_id"], r["array_pos"])] += 1
    d = sorted(ua.values())
    DEP[(sp, sam.replace(".g12", ""))] = (len(reads), n, span, d[len(d) // 2] if d else 0)
HDR = ["species_key", "sample", "context", "n_arrays", "n_units", "reads_placed", "read_unit_frags",
       "frag_spanning_pct", "unit_depth_median", "calls", "frac_mod", "read_level_median",
       "frac_5hmC", "min_array_frac", "max_array_frac"]
o = io.open(os.path.join(B, "00_catalog", "PLACEMENT_METH_REPORT.tsv"), "w", encoding="utf8", newline=NL)
o.write(TAB.join(HDR) + NL)
g = collections.defaultdict(list)
for r in AR:
    g[(r["species_key"], r["sample"], r["context"])].append(r)
print("%-22s %-13s %-4s %4s %6s %7s %9s %7s %6s %9s %s"
      % ("종", "샘플", "문맥", "arr", "유닛", "read", "호출", "깊이", "비율", "read중앙", "array 범위"))
for k in sorted(g):
    sp, sam, ctx = k
    rs = g[k]
    nu = sum(int(x["N"]) for x in rs)
    ca = sum(int(x["calls"]) for x in rs if x["calls"].isdigit())
    fr = [float(x["frac_5mC"]) for x in rs if x["frac_5mC"] not in ("", "NA")]
    w = sum(float(x["frac_5mC"]) * int(x["calls"]) for x in rs
            if x["calls"].isdigit() and x["frac_5mC"] not in ("", "NA"))
    rm = [float(x["read_level_frac_5mC_median"]) for x in rs
          if x["read_level_frac_5mC_median"] not in ("", "NA")]
    hm = [float(x["frac_5hmC"]) for x in rs if x["frac_5hmC"] not in ("", "NA")]
    d = DEP.get((sp, sam), (0, 0, 0, 0))
    row = [sp, sam, ctx, str(len(rs)), str(nu), str(d[0]), str(d[1]),
           "%.1f" % (100.0 * d[2] / d[1]) if d[1] else "-", str(d[3]), str(ca),
           "%.3f" % (w / ca) if ca else "-",
           "%.3f" % (sorted(rm)[len(rm) // 2]) if rm else "-",
           "%.3f" % (sum(hm) / len(hm)) if hm else "-",
           "%.3f" % min(fr) if fr else "-", "%.3f" % max(fr) if fr else "-"]
    o.write(TAB.join(row) + NL)
    print("%-22s %-13s %-4s %4d %6d %7s %9s %7s %6s %9s %s~%s"
          % (sp[:22], sam[:13], ctx, len(rs), nu, row[5], row[9], row[8], row[10], row[11],
             row[13], row[14]))
o.close()
