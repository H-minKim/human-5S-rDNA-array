# -*- coding: utf-8 -*-
# Build the 8 k-mer lines per species for read selection: four 22-mers of the 5S gene plus their reverse complements.
#   The input is the units in <species>.units.fa from array_struct.py; the first 119 bp of a unit is the gene (the hit interval).
#   A majority vote over several units avoids the variants of any one unit.
#   Arguments: units.fa, gene length (default 119), output file.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import sys
import collections

NL = chr(10)
FA, GL, OUT = sys.argv[1], int(sys.argv[2] if len(sys.argv) > 2 else 119), sys.argv[3]
seqs, cur = [], []
for L in io.open(FA, encoding="utf8"):
    if L.startswith(">"):
        if cur:
            seqs.append("".join(cur))
        cur = []
    else:
        cur.append(L.strip())
if cur:
    seqs.append("".join(cur))
genes = [s[:GL].upper() for s in seqs if len(s) >= GL]
cons = []
for i in range(GL):
    c = collections.Counter(g[i] for g in genes if g[i] in "ACGT")
    cons.append(c.most_common(1)[0][0] if c else "N")
gene = "".join(cons)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
ks = [gene[i:i + 22] for i in (12, 39, 62, 84)]
o = io.open(OUT, "w", encoding="utf8", newline=NL)
for k in ks:
    if "N" in k:
        continue
    o.write(k + NL)
    o.write("".join(COMP[x] for x in reversed(k)) + NL)
o.close()
print("from %d units: consensus gene %d bp, k-mer %d lines" % (len(genes), len(gene), 2 * len([k for k in ks if "N" not in k])))
print(gene)
