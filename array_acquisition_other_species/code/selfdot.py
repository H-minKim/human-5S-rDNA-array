# -*- coding: utf-8 -*-
# Read the period directly from a self dot plot, using no ruler, no annotation and none of our own units.
#   The array fragment is blastn'd against itself, and the diagonal offsets of the HSPs, sstart minus qstart, are collected weighted by aligned length.
#   A genuine tandem array gives an offset histogram standing as a ladder at every multiple of the period.
#   Human and zebra finch are run alongside as controls, to calibrate the scale.
#   Output: SELFDOT.tsv, holding the top 8 offsets per array with their total aligned length.
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import glob
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~/lit5s")
PICK = sys.argv[1] if len(sys.argv) > 1 else "rhesus,human_CHM13v2,orang_bornean,zfinch_v03_dip__CM109753,zebrafish_GRCz12tu__CM111185,zebrafish_GRCz12tu__CM111196,mouse"
TAGS = PICK.split(",")
BIN = 20          # the offset bucket
MINALN = 200      # only HSPs at least this long, which excludes short-fragment repeats
w = io.open(os.path.join(H, "SELFDOT.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "slice_len", "n_hsp", "aln_total_bp", "top_offsets(bp:정렬총량)"]) + NL)
for sfa in sorted(glob.glob(os.path.join(H, "slices", "*.arr.fa"))):
    aid = os.path.basename(sfa).replace(".arr.fa", "")
    if not any(t in aid for t in TAGS):
        continue
    hd = io.open(sfa, encoding="utf8").readline()
    kv = dict(x.split("=", 1) for x in hd.split() if "=" in x)
    slen = int(kv["slice_len"])
    r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "16", "-dust", "no", "-evalue", "1e-20",
                        "-query", sfa, "-subject", sfa, "-outfmt", "6 qstart qend sstart send pident length"],
                       capture_output=True, text=True)
    off = collections.Counter()
    n, tot = 0, 0
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        qs, qe, ss, se, pid, aln = int(f[0]), int(f[1]), int(f[2]), int(f[3]), float(f[4]), int(f[5])
        if se < ss:
            continue                      # reverse-orientation hits are excluded, not being a period
        d = ss - qs
        if d <= 0 or aln < MINALN:
            continue
        n += 1
        tot += aln
        off[int(round(d / float(BIN))) * BIN] += aln
    top = off.most_common(8)
    print("%-50s 조각 %6d · HSP %6d · 정렬 %9d" % (aid, slen, n, tot))
    print("      " + " | ".join("%d:%d" % (k, v) for k, v in top))
    w.write(TAB.join([aid, str(slen), str(n), str(tot),
                      ",".join("%d:%d" % (k, v) for k, v in top)]) + NL)
    sys.stdout.flush()
w.close()
print("표 : ~/lit5s/SELFDOT.tsv")
