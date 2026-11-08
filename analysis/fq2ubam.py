# -*- coding: utf-8 -*-
"""hits fastq(헤더에 MM/ML 태그) -> 정렬 안 한 BAM (2026-10-14)
usage: python fq2ubam.py <in.fastq.gz> <out.bam>
SEQ 는 원본 read 방향(samtools fastq -T)이라 MM/ML 과 같은 틀이다."""
import array
import gzip, sys
import pysam
src, dst = sys.argv[1], sys.argv[2]
hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrFake", "LN": 2000000}]})
out = pysam.AlignmentFile(dst, "wb", header=hdr)
n = bad = 0
with gzip.open(src, "rt") as f:
    while True:
        h = f.readline()
        if not h:
            break
        seq = f.readline().rstrip("\n"); f.readline(); q = f.readline().rstrip("\n")
        toks = h[1:].rstrip("\n").split()
        name = toks[0]
        tags = {}
        for t in toks[1:]:
            k, ty, v = t.split(":", 2)
            tags[k] = (ty, v)
        if "MM" not in tags or "ML" not in tags:
            bad += 1
            continue
        a = pysam.AlignedSegment(hdr)
        a.query_name = name
        a.query_sequence = seq
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigartuples = [(0, len(seq))]
        a.query_qualities = pysam.qualitystring_to_array(q)
        a.set_tag("MM", tags["MM"][1], "Z")
        mlv = tags["ML"][1].split(",")
        if mlv[0] == "C":
            mlv = mlv[1:]
        a.set_tag("ML", array.array("B", [int(x) for x in mlv]))
        out.write(a)
        n += 1
out.close()
print(dst, "records", n, "skipped", bad)
