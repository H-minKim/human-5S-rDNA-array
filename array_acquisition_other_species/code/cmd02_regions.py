# -*- coding: utf-8 -*-
# array 구간 + 유전자(코딩부) 구간 표 — 금화조 · 제브라(ONT 배치에 쓴 array 만)
#   array  : 5S_array_structure.tsv(좌표) ∩ AVF 에서 잰 array
#   코딩부 : cms_tbl 의 cmsearch 히트(inc '!') 중 array 안 · 겹치면 합침
#   출력   : data_v2/13_array_level_diff/REGIONS_<species>.tsv  (array_id chrom start end strand gene_intervals)
# MIT License · Copyright (c) 2026 Hungmin Kim
import collections
import csv
import io
import os

R = "$RDNA_OTHERSP_DIR/data_v2/"
OUT = R + "13_array_level_diff/"
os.makedirs(OUT, exist_ok=True)
T = "\t"
used = collections.defaultdict(set)
for r in csv.DictReader(io.open(R + "00_catalog/CMD01_ARRAY_MINUS_FLANK.tsv", encoding="utf8"), delimiter=T):
    used[r["species"]].add(r["array_id"])
st = {r["array_id"]: r for r in csv.DictReader(io.open(R + "00_catalog/5S_array_structure.tsv", encoding="utf8"), delimiter=T)}
for sp in ("zfinch_v03_dip", "zebrafish_GRCz12tu"):
    hits = collections.defaultdict(list)
    for L in io.open(R + "cms_tbl/%s.tbl" % sp, encoding="utf8"):
        if L.startswith("#"):
            continue
        f = L.split()
        if len(f) < 17 or f[16] != "!":
            continue
        a, b = int(f[7]), int(f[8])
        hits[f[0]].append((min(a, b), max(a, b), f[9]))
    o = io.open(OUT + "REGIONS_%s.tsv" % sp, "w", encoding="utf8", newline="\n")
    o.write(T.join(["array_id", "chrom", "start", "end", "gene_strand", "n_gene", "coding_bp", "gene_intervals"]) + "\n")
    for aid in sorted(used[sp]):
        s = st[aid]
        ch, a0, a1 = s["chrom"], int(s["array_start"]), int(s["array_end"])
        g = sorted((x, y, z) for x, y, z in hits[ch] if x >= a0 and y <= a1)
        merged = []
        for x, y, z in g:
            if merged and x <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], y)
            else:
                merged.append([x, y, z])
        strand = collections.Counter(z for _, _, z in g).most_common(1)[0][0] if g else "."
        cod = sum(y - x + 1 for x, y, _ in merged)
        o.write(T.join([aid, ch, str(a0), str(a1), strand, str(len(merged)), str(cod),
                        ",".join("%d-%d" % (x, y) for x, y, _ in merged)]) + "\n")
    o.close()
    print(sp, "array", len(used[sp]))
