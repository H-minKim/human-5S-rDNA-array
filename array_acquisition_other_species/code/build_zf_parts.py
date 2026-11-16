# -*- coding: utf-8 -*-
# Build a region table dividing the chr18 somatic-type array and the chr4 maternal-type arrays into coding, spacer and outside the array.
#   Coding: the 119 bp 5S gene hit intervals found by cmsearch.
#   Spacer: what remains inside the array once the coding part is removed.
#   Outside: the existing 100 kb flankL and flankR, held 5 kb clear of the boundary.
# The output is in the four-column form bismark_region_meth.py accepts, so an interval made of several pieces is written as one row per piece.
import io
import sys
import collections

TAB = chr(9)
NL = chr(10)
BASE = "$RDNA_OTHERSP_DIR/data_v2/"
HITS = BASE + "00_catalog/zebrafish_gene_types/zfgt.gene_hits.tsv"
REG = BASE + "00_catalog/zf_regions.tsv"
OUT = sys.argv[1]

# 1) read the array and the flanks from the existing region table
arrays = {}
flanks = []
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    nm, ch, s, e = f[0], f[1], int(f[2]), int(f[3])
    base, _, kind = nm.rpartition("__")
    if kind == "array":
        arrays[base] = (ch, s, e)
    else:
        flanks.append((base + "__바깥_" + kind, ch, s, e))

# 2) collect the gene hits per array
genes = collections.defaultdict(list)
for L in io.open(HITS, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if f[0] == "chrom":
        continue
    ch, gs, ge = f[0], int(f[1]), int(f[2])
    for base, (ach, as_, ae) in arrays.items():
        if ch == ach and gs >= as_ and ge <= ae:
            genes[base].append((gs, ge))

o = io.open(OUT, "w", encoding="utf8", newline=NL)
rep = []
for base in sorted(arrays):
    ch, as_, ae = arrays[base]
    g = sorted(genes.get(base, []))
    # merge overlaps
    mg = []
    for gs, ge in g:
        if mg and gs <= mg[-1][1]:
            mg[-1] = (mg[-1][0], max(mg[-1][1], ge))
        else:
            mg.append((gs, ge))
    cod = sum(e - s + 1 for s, e in mg)
    for gs, ge in mg:
        o.write(TAB.join([base + "__코딩부", ch, str(gs), str(ge)]) + NL)
    # the spacer is the array with the coding part removed
    p = as_
    sp = 0
    for gs, ge in mg:
        if gs > p:
            o.write(TAB.join([base + "__스페이서", ch, str(p), str(gs - 1)]) + NL)
            sp += gs - p
        p = max(p, ge + 1)
    if p <= ae:
        o.write(TAB.join([base + "__스페이서", ch, str(p), str(ae)]) + NL)
        sp += ae - p + 1
    rep.append((base, len(mg), cod, sp))
for nm, ch, s, e in flanks:
    o.write(TAB.join([nm, ch, str(s), str(e)]) + NL)
o.close()
print("array %d · flank 조각 %d" % (len(rep), len(flanks)))
for base, n, cod, sp in rep:
    print("  %-34s 유전자 %3d 개 · 코딩 %6d bp · 스페이서 %7d bp" % (base, n, cod, sp))
