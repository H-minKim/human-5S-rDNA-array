# -*- coding: utf-8 -*-
# 전 종 read 배치 결과 + 메틸 요약 (타 생물종 전용)
#   배치 : 03_reads 의 array 요약과 06_units_from_reads 의 조각 표에서 센다
#   메틸 : 같은 array 요약의 CpG(식물은 CG/CHG/CHH) 호출 수와 비율
#   산출 : PLACEMENT_METH_REPORT.tsv (종×샘플×문맥) · 화면에 종 단위 요약
#   인자 : data_v2 폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
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
# read 조각 표에서 배치 깊이
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
