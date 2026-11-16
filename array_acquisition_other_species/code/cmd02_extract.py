# -*- coding: utf-8 -*-
# ONT batch BAM (MM/ML) -> a table of CpG calls by read x array x class (coding / spacer).
#   CpG-only models alone (zebra finch dorado 5mCG_5hmCG, zebrafish guppy 5mC_CG), so a 5mC ('m') call is taken as a CpG without deciding the context.
#   A site is called methylated at p(5mC) = ML/255 >= 0.5, the same threshold modkit uses in the batch path.
#   Class: a CpG whose C lies inside a cmsearch gene interval is coding; anything else inside the array is spacer.
#   Supplementary and secondary alignments are excluded; there is no mapq filter (mapq 0 is normal inside a tandem array) and mapq is kept as a column.
#   Arguments: bam, sample, REGIONS.tsv, out.tsv.
# MIT License · Copyright (c) 2026 Hungmin Kim
import bisect
import io
import os
import sys
import pysam

T = "\t"
BAM, LAB, REG, OUT = sys.argv[1:5]
regs = []
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip("\n").split(T)
    if f[0] == "array_id":
        continue
    gi = []
    for x in f[7].split(","):
        if x:
            a, b = x.split("-")
            gi.append((int(a), int(b)))
    regs.append((f[0], f[1], int(f[2]), int(f[3]), gi))
b = pysam.AlignmentFile(BAM, "rb")
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline="\n")
if new:
    o.write(T.join(["sample", "array_id", "read_id", "mapq", "part", "n_calls", "n_meth"]) + "\n")
tot = 0
for aid, ch, a0, a1, gi in regs:
    gs = [x for x, _ in gi]
    ge = [y for _, y in gi]

    def coding(p):
        i = bisect.bisect_right(gs, p) - 1
        return i >= 0 and p <= ge[i]

    for r in b.fetch(ch, a0 - 1, a1):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        mb = r.modified_bases
        if not mb:
            continue
        q2r = dict(r.get_aligned_pairs(matches_only=True))
        cnt = {"코딩부": [0, 0], "스페이서": [0, 0]}
        for key, calls in mb.items():
            if key[2] != "m":
                continue
            rev = key[1] == 1
            for qp, ql in calls:
                rp = q2r.get(qp)
                if rp is None:
                    continue
                cpos = rp if not rev else rp - 1      # the C of the CpG (add 1 for 1-based)
                p1 = cpos + 1
                if p1 < a0 or p1 > a1:
                    continue
                part = "코딩부" if coding(p1) else "스페이서"
                cnt[part][0] += 1
                if ql / 255.0 >= 0.5:
                    cnt[part][1] += 1
        for part, (c, m) in cnt.items():
            if c:
                o.write(T.join([LAB, aid, r.query_name, str(r.mapping_quality), part, str(c), str(m)]) + "\n")
                tot += 1
o.close()
b.close()
print("%s 행 %d" % (LAB, tot))
