# -*- coding: utf-8 -*-
# For each group BAM, extract per read: locus class, region class, CpG count, methylated count. The bootstrap and the between-group test are done locally.
#   Locus classes: chr4 (a cluster of 6 arrays) and chr18.
#   Region classes: coding, spacer, outside.
#   Arguments: regions.tsv out.tsv bam:group [bam:group ...]
# MIT License. Copyright (c) 2026 Hyungmin Kim
import bisect
import collections
import io
import sys
import pysam

TAB = chr(9)
NL = chr(10)
REG, OUT = sys.argv[1], sys.argv[2]
parts = collections.defaultdict(list)   # (locus class, region class) -> [(chrom, start, end)]
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    base, _, kind = f[0].rpartition("__")
    br = "chr18" if base.startswith("chr18") else "chr4"
    k = kind
    parts[(br, k)].append((f[1], int(f[2]), int(f[3])))

o = io.open(OUT, "w", encoding="utf8", newline=NL)
o.write(TAB.join(["group", "branch", "part", "read", "n_cpg", "n_meth"]) + NL)
for arg in sys.argv[3:]:
    bam, grp = arg.split(":")
    b = pysam.AlignmentFile(bam, "rb")
    for (br, k), iv in sorted(parts.items()):
        bych = collections.defaultdict(list)
        for ch, s, e in iv:
            bych[ch].append((s, e))
        for ch, ivs in bych.items():
            ivs.sort()
            st = [s for s, e in ivs]
            en = [e for s, e in ivs]
            # The intervals are scattered, so fetch per interval but count each read only once.
            seen = {}
            for s, e in ivs:
                for r in b.fetch(ch, s - 1, e):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                        continue
                    key = (r.query_name, r.reference_start)
                    if key in seen:
                        continue
                    try:
                        xm = r.get_tag("XM")
                    except KeyError:
                        continue
                    c = m = 0
                    for qp, rp in r.get_aligned_pairs(matches_only=True):
                        if rp is None or qp is None:
                            continue
                        p = rp + 1
                        i = bisect.bisect_right(st, p) - 1
                        if i < 0 or p > en[i]:
                            continue
                        z = xm[qp] if qp < len(xm) else "."
                        if z in ("z", "Z"):
                            c += 1
                            m += (z == "Z")
                    seen[key] = (c, m)
            for (qn, rs), (c, m) in seen.items():
                if c:
                    o.write(TAB.join([grp, br, k, "%s:%d" % (qn, rs), str(c), str(m)]) + NL)
    b.close()
o.close()
print("done")
