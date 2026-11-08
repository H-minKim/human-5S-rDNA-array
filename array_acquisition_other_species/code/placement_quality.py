# -*- coding: utf-8 -*-
# array 마다 배치 건전성 · 깊이 · 관통 read 를 잰다
#
#   깊이      = array 구간에 정렬된 염기 합 / array 길이 (읽은 깊이)
#   관통 read = 한 분자가 array 전체(시작~끝)를 덮은 것. 양쪽 flank 까지 물면 flank-to-flank.
#   junction  = array 경계를 물고 바깥 flank 를 >= JFLANK bp 가진 read (한쪽 경계)
#   배치 건전성 축
#     · mapq0 비율        — 탠덤이라 높은 게 정상이지만, 얼마나 높은지는 적어 둔다
# - clip fraction: if many reads are heavily clipped at both ends, the placement is doubtful
#     · 정렬 일치도 중앙값  — NM 태그로 낸 read-참조 일치도
#   mapq 로 read 를 거르지 않는다.
#
#   인자 : 종키 샘플 bam chrom array시작 array끝 [JFLANK=1000]
# MIT License · Copyright (c) 2026 Hungmin Kim
import sys
import pysam

SP, SAMPLE, BAM, CHROM = sys.argv[1:5]
A0, A1 = int(sys.argv[5]), int(sys.argv[6])
JF = int(sys.argv[7]) if len(sys.argv) > 7 else 1000
AL = A1 - A0 + 1
b = pysam.AlignmentFile(BAM, "rb")
nread = nspan = njL = njR = nmq0 = 0
bases = 0
ident = []
clip = []
for r in b.fetch(CHROM, max(0, A0 - 1 - 50000), A1 + 50000):
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue
    s, e = r.reference_start + 1, r.reference_end
    if e < A0 or s > A1:
        continue
    nread += 1
    if r.mapping_quality == 0:
        nmq0 += 1
    bases += min(e, A1) - max(s, A0) + 1
    if s <= A0 and e >= A1:
        nspan += 1
    if s <= A0 - JF and e >= A0 + JF:
        njL += 1
    if s <= A1 - JF and e >= A1 + JF:
        njR += 1
    nm = r.get_tag("NM") if r.has_tag("NM") else None
    al = r.query_alignment_length or 1
    if nm is not None:
        ident.append(100.0 * (1 - nm / float(al)))
    cl = 0
    for op, ln in (r.cigartuples or []):
        if op in (4, 5):
            cl += ln
    clip.append(100.0 * cl / max(1, (r.query_length or 1)))
b.close()
ident.sort()
clip.sort()
print("\t".join(str(x) for x in [
    SP, SAMPLE, CHROM, A0, A1, AL, nread,
    "%.1f" % (bases / float(AL)),
    nspan, njL, njR,
    "%.1f" % (100.0 * nmq0 / nread) if nread else "-",
    "%.1f" % ident[len(ident) // 2] if ident else "-",
    "%.1f" % clip[len(clip) // 2] if clip else "-"]))
