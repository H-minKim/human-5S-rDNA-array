# -*- coding: utf-8 -*-
# Measure the soundness of the placement, the depth, and the array-spanning reads, per array.
#
#   Depth is the total aligned bases over the array interval divided by the array length, that is the read depth.
#   An array-spanning read is one molecule covering the array from start to end; reaching both flanks as well makes it flank-to-flank.
#   A junction read crosses the array boundary and carries at least JFLANK bp of outside flank, on one side.
#   The axes of placement soundness:
#     · the mapq 0 fraction, which is normally high in a tandem array, but how high is recorded;
# - clip fraction: if many reads are heavily clipped at both ends, the placement is doubtful
#     · the median alignment identity, the read-to-reference identity computed from the NM tag.
#   Reads are not filtered by mapq.
#
#   Arguments: species key, sample, bam, chrom, array start, array end, [JFLANK = 1000].
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
