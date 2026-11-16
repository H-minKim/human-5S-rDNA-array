# -*- coding: utf-8 -*-
# Report CpG methylation per region from a bismark BAM, for comparing an array with its flank.
#
#   bismark carries the methylation calls in the XM tag, not in MM/ML.
#     Z = methylated CpG, z = unmethylated CpG, X/x = CHG, H/h = CHH, U/u = unknown, . = not a C.
#   It is therefore kept separate from region_meth.py, which is for ONT.
#
#   No mapq filter: mapq 0 is normal in a tandem array.
#     The value using only mapq >= 10 is reported alongside, so any bias is visible.
#   Reads flagged as duplicates are excluded.
#
#   Arguments: bam, label, regions.tsv, out.tsv.
#     The region list is name, chrom, start, end, tab separated.
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import pysam

TAB = chr(9)
NL = chr(10)
BAM, LAB, REG, OUT = sys.argv[1:5]
b = pysam.AlignmentFile(BAM, "rb")
REFS = set(b.references)
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline=NL)
if new:
    o.write(TAB.join(["sample", "region", "chrom", "start", "end", "reads",
                      "n_cpg", "n_meth", "frac_CpG", "n_cpg_mq10", "frac_CpG_mq10",
                      "n_chg", "frac_CHG", "n_chh", "frac_CHH"]) + NL)
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    nm, ch, s, e = f[0], f[1], int(f[2]), int(f[3])
    if ch not in REFS:
        o.write(TAB.join([LAB, nm, ch, str(s), str(e)] + ["-"] * 10) + NL)
        continue
    nr = 0
    cg = me = 0
    cgq = meq = 0
    chg = chgm = 0
    chh = chhm = 0
    for r in b.fetch(ch, max(0, s - 1), e):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        try:
            xm = r.get_tag("XM")
        except KeyError:
            continue
        nr += 1
        q = r.mapping_quality
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp is None or qp is None or not (s <= rp + 1 <= e):
                continue
            c = xm[qp] if qp < len(xm) else "."
            if c in ("z", "Z"):
                cg += 1
                if c == "Z":
                    me += 1
                if q >= 10:
                    cgq += 1
                    if c == "Z":
                        meq += 1
            elif c in ("x", "X"):
                chg += 1
                if c == "X":
                    chgm += 1
            elif c in ("h", "H"):
                chh += 1
                if c == "H":
                    chhm += 1
    o.write(TAB.join([LAB, nm, ch, str(s), str(e), str(nr), str(cg), str(me),
                      "%.4f" % (me / float(cg)) if cg else "-",
                      str(cgq), "%.4f" % (meq / float(cgq)) if cgq else "-",
                      str(chg), "%.4f" % (chgm / float(chg)) if chg else "-",
                      str(chh), "%.4f" % (chhm / float(chh)) if chh else "-"]) + NL)
o.close()
b.close()
print("%s 완료" % LAB)
