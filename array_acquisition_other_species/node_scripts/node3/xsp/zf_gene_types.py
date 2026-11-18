# -*- coding: utf-8 -*-
# Divide the zebrafish 5S gene sequences into types, exhaustively.
#   Every cmsearch (RF00001) hit is taken, its sequence extracted, and the number of each gene type counted per array.
#   The actual gene hit sequence is used, not the last 119 bp of the array's representative unit, so no boundary is assumed.
#   Identical sequences are grouped, and the number of positions each differs from the most common type, the candidate maternal type, is recorded.
#   Arguments: cmsearch.tbl, assembly.fa, array catalogue.tsv, species key, output prefix.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import csv
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
TBL, FA, CAT, SP, OUT = sys.argv[1:6]
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))
hits = []
for L in io.open(TBL, encoding="utf8"):
    if L.startswith("#"):
        continue
    f = L.split()
    if len(f) < 15:
        continue
    ch, s, e, st, sc = f[0], int(f[7]), int(f[8]), f[9], float(f[14])
    if st == "-":
        s, e = e, s
    hits.append((ch, s, e, st, sc))
print("cmsearch hit %d" % len(hits))
reg = OUT + ".reg.txt"
io.open(reg, "w", encoding="utf8", newline=NL).write(
    NL.join("%s:%d-%d" % (h[0], h[1], h[2]) for h in hits) + NL)
r = subprocess.run(["samtools", "faidx", "-r", reg, FA], capture_output=True, text=True)
seqs, cur = [], []
for L in r.stdout.splitlines():
    if L.startswith(">"):
        if cur:
            seqs.append("".join(cur))
        cur = []
    else:
        cur.append(L.strip())
if cur:
    seqs.append("".join(cur))
os.remove(reg)
G = []
for h, s in zip(hits, seqs):
    G.append((h, rc(s.upper()) if h[3] == "-" else s.upper()))
# assignment to an array
AR = []
for r2 in csv.DictReader(io.open(CAT, encoding="utf8"), delimiter=TAB):
    if r2["species_key"] != SP:
        continue
    AR.append((r2["chrom"], int(r2["start"]), int(r2["end"]), r2["array_id"], int(r2["n_copies"])))
byc = collections.defaultdict(list)
for a in AR:
    byc[a[0]].append(a)


def which(ch, s, e):
    for c, a0, a1, aid, n in byc.get(ch, []):
        if s >= a0 - 500 and e <= a1 + 500:
            return aid
    return "-"


cnt = collections.Counter(g[1] for g in G)
top = cnt.most_common(1)[0][0]
grp = {}
for i, (seq, c) in enumerate(cnt.most_common()):
    grp[seq] = "T%02d" % (i + 1)
o = io.open(OUT + ".gene_hits.tsv", "w", encoding="utf8", newline=NL)
o.write(TAB.join(["chrom", "start", "end", "strand", "score", "array_id", "gene_type",
                  "gene_len", "diff_vs_T01"]) + NL)
for (h, seq) in G:
    aid = which(h[0], h[1], h[2])
    d = sum(1 for x, y in zip(seq, top) if x != y) if len(seq) == len(top) else -1
    o.write(TAB.join([h[0], str(h[1]), str(h[2]), h[3], "%.1f" % h[4], aid, grp[seq],
                      str(len(seq)), str(d) if d >= 0 else "len%d" % len(seq)]) + NL)
o.close()
w = io.open(OUT + ".gene_types.tsv", "w", encoding="utf8", newline=NL)
w.write(TAB.join(["gene_type", "n_hits", "gene_len", "diff_vs_T01", "chroms", "arrays", "seq"]) + NL)
print("%-6s %6s %6s %8s  %s" % ("type", "hit", "길이", "T01차이", "염색체(array 수)"))
for seq, c in cnt.most_common():
    t = grp[seq]
    hh = [(h, s) for h, s in G if s == seq]
    ch = collections.Counter(h[0] for h, _s in hh)
    ar = set(which(h[0], h[1], h[2]) for h, _s in hh)
    d = sum(1 for x, y in zip(seq, top) if x != y) if len(seq) == len(top) else -1
    w.write(TAB.join([t, str(c), str(len(seq)), str(d) if d >= 0 else "len", ";".join(ch),
                      str(len(ar)), seq]) + NL)
    if c >= 3 or d in (-1,) or d > 3:
        print("%-6s %6d %6d %8s  %s" % (t, c, len(seq), d if d >= 0 else "길이다름",
                                        dict(ch)))
w.close()
