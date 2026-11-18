# -*- coding: utf-8 -*-
# Preprocessing of single-end PBAT reads, in place of cutadapt, which is not installed on the node.
#   1) Trim 9 bp from the 5' end: PBAT is primed randomly, so the start of the read does not match the reference.
#   2) Trim the 3' adapter, the common Illumina AGATCGGAAGAGC; the adapter is ligated after bisulfite treatment and so is not converted.
#      If there is no full match, a prefix of at least 5 bp running off the end of the read is also trimmed.
#   3) Discard reads shorter than 30 bp.
#   Arguments: in.fastq.gz out.fastq.gz
# MIT License. Copyright (c) 2026 Hyungmin Kim
import gzip
import sys

ADP = "AGATCGGAAGAGC"
CLIP5 = 9
MINL = 30
src, dst = sys.argv[1], sys.argv[2]
n = kept = adp = 0
with gzip.open(src, "rt") as f, gzip.open(dst, "wt") as o:
    while True:
        h = f.readline()
        if not h:
            break
        s = f.readline().rstrip("\n")
        f.readline()
        q = f.readline().rstrip("\n")
        n += 1
        i = s.find(ADP)
        if i < 0:
            for k in range(len(ADP) - 1, 4, -1):
                if s.endswith(ADP[:k]):
                    i = len(s) - k
                    break
        if i >= 0:
            adp += 1
            s, q = s[:i], q[:i]
        s, q = s[CLIP5:], q[CLIP5:]
        if len(s) < MINL:
            continue
        kept += 1
        o.write(h + s + "\n+\n" + q + "\n")
sys.stderr.write("%s read %d · 어댑터 %d(%.1f%%) · 남김 %d(%.1f%%)\n" % (
    src.split("/")[-1], n, adp, 100.0 * adp / max(n, 1), kept, 100.0 * kept / max(n, 1)))
