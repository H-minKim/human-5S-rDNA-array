# -*- coding: utf-8 -*-
# Methylation per read: CpG sites on the same read are not independent.
#   A binomial interval over CpG counts takes the sample size to be the number of CpG sites rather than the number of reads,
#   and so understates the uncertainty badly. Measured: only 17 reads cover the chr18 array in the Potok egg sample,
#   yet 82 CpG sites gave an interval of +-0.107.
#   So one CpG methylation ratio is computed per read, and the mean and interval are taken over reads as the sampling unit.
#
#   Arguments: bam, label, regions.tsv, out.tsv, [minimum CpG sites per read (default 3)].
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import math
import os
import sys
import pysam

TAB = chr(9)
NL = chr(10)
BAM, LAB, REG, OUT = sys.argv[1:5]
MINC = int(sys.argv[5]) if len(sys.argv) > 5 else 3

b = pysam.AlignmentFile(BAM, "rb")
REFS = set(b.references)
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline=NL)
if new:
    o.write(TAB.join(["sample", "region", "n_read", "read_mean",
                      "read_sd", "read_se", "n_cpg_total"]) + NL)

for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    nm, ch, s, e = f[0], f[1], int(f[2]), int(f[3])
    if ch not in REFS:
        o.write(TAB.join([LAB, nm] + ["-"] * 5) + NL)
        continue
    vals = []
    tot = 0
    for r in b.fetch(ch, max(0, s - 1), e):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        try:
            xm = r.get_tag("XM")
        except KeyError:
            continue
        cg = me = 0
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp is None or qp is None or not (s <= rp + 1 <= e):
                continue
            c = xm[qp] if qp < len(xm) else "."
            if c in ("z", "Z"):
                cg += 1
                if c == "Z":
                    me += 1
        if cg >= MINC:
            vals.append(me / float(cg))
            tot += cg
    n = len(vals)
    if n == 0:
        o.write(TAB.join([LAB, nm, "0"] + ["-"] * 4) + NL)
        continue
    mu = sum(vals) / n
    if n > 1:
        sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / (n - 1))
        se = sd / math.sqrt(n)
    else:
        sd = se = float("nan")
    o.write(TAB.join([LAB, nm, str(n), "%.4f" % mu,
                      "%.4f" % sd, "%.4f" % se, str(tot)]) + NL)
o.close()
b.close()
print("%s done" % LAB)
