# -*- coding: utf-8 -*-
# Can reads be assigned in an array with duplicate haplotypes?
#
#   In the zebra finch the maternal and paternal chr2 arrays are byte-identical in sequence,
#   so a read covering only the interior of the array matches both equally and cannot be assigned in principle, giving mapq 0.
#   A read reaching into the flank, however, should be assignable, because the flank sequences differ.
#   This script tests that prediction by splitting the reads into two groups and comparing their mapq distributions:
#     group A: the interior only, with less than ANCHOR bp of flank;
#     group B: at least ANCHOR bp of flank on at least one side.
#   It also counts whether the reads landing on the maternal and paternal copies split evenly, which would be the mark of the aligner dividing them arbitrarily.
#   Arguments: bam, chrom, array start, array end, label, [ANCHOR = 1000].
# MIT License · Copyright (c) 2026 Hungmin Kim
import sys
import collections
import pysam

BAM, CHROM = sys.argv[1], sys.argv[2]
A0, A1 = int(sys.argv[3]), int(sys.argv[4])
LAB = sys.argv[5]
ANCHOR = int(sys.argv[6]) if len(sys.argv) > 6 else 1000
b = pysam.AlignmentFile(BAM, "rb")
grp = collections.defaultdict(lambda: collections.Counter())
mq = collections.defaultdict(list)
for r in b.fetch(CHROM, max(0, A0 - 1 - 50000), A1 + 50000):
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue
    s, e = r.reference_start + 1, r.reference_end
    if e < A0 or s > A1:
        continue
    left = max(0, A0 - s)      # the length reaching outside the left edge of the array
    right = max(0, e - A1)
    k = "B_flank앵커" if (left >= ANCHOR or right >= ANCHOR) else "A_array내부"
    grp[k]["n"] += 1
    grp[k]["mq0"] += 1 if r.mapping_quality == 0 else 0
    grp[k]["mq60"] += 1 if r.mapping_quality >= 60 else 0
    grp[k]["mq1_59"] += 1 if 0 < r.mapping_quality < 60 else 0
    mq[k].append(r.mapping_quality)
b.close()
print("%-26s %-12s %6s %7s %8s %8s %9s" % (LAB, "무리", "read", "mapq0", "mapq1-59", "mapq>=60", "mapq중앙"))
for k in sorted(grp):
    g = grp[k]
    v = sorted(mq[k])
    print("%-26s %-12s %6d %6.1f%% %7.1f%% %7.1f%% %9d"
          % ("", k, g["n"], 100.0 * g["mq0"] / g["n"], 100.0 * g["mq1_59"] / g["n"],
             100.0 * g["mq60"] / g["n"], v[len(v) // 2]))
