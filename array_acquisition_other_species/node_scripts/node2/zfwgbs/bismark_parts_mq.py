# -*- coding: utf-8 -*-
# Report CpG methylation with the array split into coding, spacer and outside (flank).
#
#   The intervals are highly fragmented (65 coding blocks on chr18), so instead of fetching per interval
#     the array plus flank is traversed once and each CpG site is assigned to a class.
#
#   The interval comes from a read-level bootstrap: CpG sites on the same read are not independent.
#     Reads are resampled with replacement, the ratio is recomputed, and the 2.5 / 97.5 percentiles are used.
#     A coding block can hold only a handful of CpG sites, so a binomial interval would be far too narrow.
#
#   Arguments: bam, label, regions.tsv, out.tsv, [bootstrap replicates (default 2000)].
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import random
import sys
import collections
import pysam

TAB = chr(9)
NL = chr(10)
BAM, LAB, REG, OUT = sys.argv[1:5]
NB = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
random.seed(0)

# region blocks -> the list of intervals per class
parts = collections.defaultdict(list)          # (base, kind) -> [(chrom,s,e)]
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    base, _, kind = f[0].rpartition("__")
    parts[(base, kind)].append((f[1], int(f[2]), int(f[3])))

b = pysam.AlignmentFile(BAM, "rb")
REFS = set(b.references)
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline=NL)
if new:
    o.write(TAB.join(["sample", "locus", "part", "n_read", "n_cpg", "n_meth",
                      "frac", "lo95", "hi95", "bp"]) + NL)


MQ=int(__import__("os").environ.get("MQ","0"))
def boot(per_read):
    """resample with replacement by read: the CpGs of one read move together"""
    n = len(per_read)
    if n < 2:
        return ("-", "-")
    out = []
    for _ in range(NB):
        c = m = 0
        for _ in range(n):
            a, d = per_read[random.randrange(n)]
            c += a
            m += d
        if c:
            out.append(m / float(c))
    if not out:
        return ("-", "-")
    out.sort()
    return ("%.4f" % out[int(0.025 * len(out))], "%.4f" % out[int(0.975 * len(out))])


for (base, kind) in sorted(parts):
    iv = parts[(base, kind)]
    ch = iv[0][0]
    if ch not in REFS:
        o.write(TAB.join([LAB, base, kind] + ["-"] * 7) + NL)
        continue
    iv_s = sorted((s, e) for _, s, e in iv)
    bp = sum(e - s + 1 for s, e in iv_s)
    # For speed, a binary search over sorted intervals is used instead of a set of positions.
    import bisect
    starts = [s for s, e in iv_s]
    ends = [e for s, e in iv_s]

    def inside(p):
        i = bisect.bisect_right(starts, p) - 1
        return i >= 0 and p <= ends[i]

    lo, hi = iv_s[0][0], iv_s[-1][1]
    per_read = []
    tot_c = tot_m = 0
    for r in b.fetch(ch, max(0, lo - 1), hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.mapping_quality < MQ:
            continue
        try:
            xm = r.get_tag("XM")
        except KeyError:
            continue
        c = m = 0
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp is None or qp is None:
                continue
            if not inside(rp + 1):
                continue
            ch2 = xm[qp] if qp < len(xm) else "."
            if ch2 in ("z", "Z"):
                c += 1
                if ch2 == "Z":
                    m += 1
        if c:
            per_read.append((c, m))
            tot_c += c
            tot_m += m
    fr = "%.4f" % (tot_m / float(tot_c)) if tot_c else "-"
    l, h = boot(per_read)
    o.write(TAB.join([LAB, base, kind, str(len(per_read)), str(tot_c), str(tot_m),
                      fr, l, h, str(bp)]) + NL)
o.close()
b.close()
print("%s done" % LAB)
