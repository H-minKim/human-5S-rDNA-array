# -*- coding: utf-8 -*-
# 보노보 5S array 메틸을 호출기 판본별로 다시 잰다(기본값 확정용)
#   입력 : 03_reads/bonobo…read_unit.tsv(read×유닛 · n_calls · n_5mC) + read→판본 대응(flowcell 선별 SAM 의 read id)
#   값 : 호출 가중 5mC 비율(Σn_5mC/Σn_calls) · read 중앙값 · read 수
# MIT License · Copyright (c) 2026 Hungmin Kim
import collections
import csv
import io
import statistics
import sys

RU = "$RDNA_OTHERSP_DIR/data_v2/03_reads/bonobo_mPanPan1/bonobo_mPanPan1__bonobo_mPanPan1.read_unit.tsv"
MAP = sys.argv[1]
OUT = "$RDNA_OTHERSP_DIR/data_v2/00_catalog/CMD01_BONOBO_BY_BASECALLER.tsv"
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
    print("%-12s read %4d · 호출 %8d · 호출가중 %.3f · read 중앙 %.3f" % (k, len(v), c, m / c, med))
o.close()
print("reads of unknown version", miss)
