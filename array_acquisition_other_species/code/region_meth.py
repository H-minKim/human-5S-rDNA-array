# -*- coding: utf-8 -*-
# BAM 의 한 구역에서 CpG 메틸 비율을 낸다 (MM/ML 직접 · mapq 로 거르지 않음)
#   같은 BAM 안에서 array 와 옆 구역을 비교하는 내부 대조용이다.
#   조직·호출기·파이프라인이 같으므로 차이가 곧 구역의 차이다.
#   인자 : bam chrom start end [라벨]
# MIT License · Copyright (c) 2026 Hungmin Kim
import sys
import pysam

BAM, CHROM = sys.argv[1], sys.argv[2]
LO, HI = int(sys.argv[3]), int(sys.argv[4])
LAB = sys.argv[5] if len(sys.argv) > 5 else ""
b = pysam.AlignmentFile(BAM, "rb")
n = m = h = 0
nread = 0
for r in b.fetch(CHROM, LO - 1, HI):
    if r.is_secondary or r.is_supplementary or r.is_unmapped:
        continue
    mb = r.modified_bases
    if not mb:
        continue
    nread += 1
    ap = dict((q, rp) for q, rp in r.get_aligned_pairs(matches_only=True))
    best = {}
    for k, v in mb.items():
        if k[2] not in ("m", "h"):
            continue
        for pos, ql in v:
            p = best.get(pos, [0.0, 0.0])
            if k[2] == "m":
                p[0] = max(p[0], ql / 255.0)
            else:
                p[1] = max(p[1], ql / 255.0)
            best[pos] = p
    for pos, (pm, ph) in best.items():
        rp = ap.get(pos)
        if rp is None or not (LO <= rp + 1 <= HI):
            continue
        n += 1
        if pm >= 0.5:
            m += 1
        if ph >= 0.5:
            h += 1
print("%-28s %s:%d-%d  read %4d · CpG %8d · 5mC %.3f · 5hmC %.3f"
      % (LAB, CHROM, LO, HI, nread, n, m / float(n) if n else 0, h / float(n) if n else 0))
