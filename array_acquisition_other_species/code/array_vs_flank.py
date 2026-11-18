# -*- coding: utf-8 -*-
# Measure the methylation of an array and of the 100 kb of flank on each side, from the same BAM.
#
#   The question: is the difference in methylation between arrays due to the surrounding local chromatin, or is it a property of the array itself?
#     - If it is the surroundings, the array and its flank move together, and a hypomethylated array has a hypomethylated flank.
#     - If it is the array, the flanks are similar and only the array falls away.
#   Measuring from the same BAM, the same tissue and the same basecaller removes those effects.
#   The flank is taken at GAP bp clear of the array boundary and spans WIN bp, which avoids the boundary effect.
#   No mapq filter; only supplementary and secondary alignments are removed.
#
#   Arguments: species key, sample, bam, array coordinates .tsv, out.tsv, [WIN = 100000], [GAP = 5000].
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import csv
import pysam

TAB = chr(9)
NL = chr(10)
SP, SAMPLE, BAM, COORD, OUT = sys.argv[1:6]
WIN = int(sys.argv[6]) if len(sys.argv) > 6 else 100000
GAP = int(sys.argv[7]) if len(sys.argv) > 7 else 5000
b = pysam.AlignmentFile(BAM, "rb")
REF = set(b.references)


def meth(chrom, lo, hi):
    if chrom not in REF or hi <= lo:
        return (0, 0, 0, 0)
    n = m = h = 0
    nr = 0
    for r in b.fetch(chrom, max(0, lo - 1), hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        mb = r.modified_bases
        if not mb:
            continue
        nr += 1
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
            if rp is None or not (lo <= rp + 1 <= hi):
                continue
            n += 1
            if pm >= 0.5:
                m += 1
            if ph >= 0.5:
                h += 1
    return (nr, n, m, h)


HDR = ["species_key", "sample", "array_id", "chrom", "array_start", "array_end", "array_len",
       "region", "reg_start", "reg_end", "reads", "n_c", "n_5mC", "frac_5mC", "n_5hmC", "frac_5hmC"]
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline=NL)
if new:
    o.write(TAB.join(HDR) + NL)
for L in io.open(COORD, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 5 or f[0] != SP:
        continue
    aid, ch, s, e = f[1], f[2], int(f[3]), int(f[4])
    regs = [("array", s, e),
            ("flank_left_100kb", max(1, s - GAP - WIN), max(1, s - GAP)),
            ("flank_right_100kb", e + GAP, e + GAP + WIN)]
    for nm, lo, hi in regs:
        nr, n, m, h = meth(ch, lo, hi)
        o.write(TAB.join([SP, SAMPLE, aid, ch, str(s), str(e), str(e - s + 1), nm, str(lo), str(hi),
                          str(nr), str(n), str(m),
                          "%.3f" % (m / float(n)) if n else "-", str(h),
                          "%.4f" % (h / float(n)) if n else "-"]) + NL)
o.close()
b.close()
print("%s %s done" % (SP, SAMPLE))
