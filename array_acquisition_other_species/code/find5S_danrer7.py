# -*- coding: utf-8 -*-
# Find the maternal-type (chr4) and somatic-type (chr18) 5S genes in the danRer7 chromosomes by sequence search, rather than by lifting coordinates.
#   The 119 bp of each type is slid along both strands, positions with at most 12 mismatches are collected, and each is assigned to the nearer type.
#   For speed, only candidates are examined, found by a 16-mer seed from the last 60 bp the two types share, positions 54-119, where they do not differ.
#   Arguments: chr.fa.gz out.tsv
import gzip
import io
import sys

TYPES = {}
name = None
for L in io.open(sys.argv[3] if len(sys.argv) > 3 else "../gsc/zf5S_types.fa"):
    L = L.strip()
    if L.startswith(">"):
        name = L[1:]
    elif L:
        TYPES[name] = L
RC = str.maketrans("ACGTN", "TGCAN")
seq = []
chrom = None
with gzip.open(sys.argv[1], "rt") as f:
    for L in f:
        if L.startswith(">"):
            chrom = L[1:].strip().split()[0]
        else:
            seq.append(L.strip().upper())
g = "".join(seq)
N = len(g)
mat = TYPES["maternal_chr4"]
som = TYPES["somatic_chr18"]
seeds = set()
for t in (mat, som):
    for i in range(60, 119 - 16 + 1):
        seeds.add(t[i:i + 16])
SEEDPOS = mat.find("GGGAGACCGCCTGGGA")
assert SEEDPOS == som.find("GGGAGACCGCCTGGGA") and SEEDPOS >= 0
out = io.open(sys.argv[2], "a", encoding="utf8", newline="\n")
seen = set()
hits = 0


def mm(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


for strand in "+-":
    s = g if strand == "+" else g.translate(RC)[::-1]
    i = 0
    while True:
        # the search uses a single seed, from the most conserved stretch, for speed
        j = s.find("GGGAGACCGCCTGGGA", i)
        if j < 0:
            break
        i = j + 1
        # the seed position is computed from the sequence, not assumed
        st = j - SEEDPOS
        if st < 0 or st + 119 > len(s):
            continue
        w = s[st:st + 119]
        dm, ds = mm(w, mat), mm(w, som)
        best = min(dm, ds)
        if best > 12:
            continue
        typ = "maternal" if dm < ds else ("somatic" if ds < dm else "tie")
        if strand == "+":
            a, b = st + 1, st + 119
        else:
            a, b = N - (st + 119) + 1, N - st
        if (a, b) in seen:
            continue
        seen.add((a, b))
        hits += 1
        out.write("\t".join([chrom, str(a), str(b), strand, typ, str(dm), str(ds)]) + "\n")
out.close()
print(chrom, "hits", hits)
