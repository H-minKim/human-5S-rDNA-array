# -*- coding: utf-8 -*-
# Store the offset histogram of the array fragment's self-alignment as it is; the call is made separately.
#   The base period and its multiple, the higher-order repeat, must be reported separately: taking whichever of the two carries more aligned length
#                 records 1,520, that is 3 units, as the period of an array with a 500 bp unit.
#   Output: OFFHIST.tsv, holding array_id, the fragment length, the HSP count and every "offset:aligned bp" pair, for offsets up to 30000 and at least 300 bp aligned.
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
BIN = 20
MINALN = 200
w = io.open(os.path.join(H, "OFFHIST.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "slice_len", "n_hsp", "hist"]) + NL)
for sfa in sorted(glob.glob(os.path.join(H, "slices", "*.arr.fa"))):
    aid = os.path.basename(sfa).replace(".arr.fa", "")
    kv = dict(x.split("=", 1) for x in io.open(sfa, encoding="utf8").readline().split() if "=" in x)
    slen = int(kv["slice_len"])
    r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "16", "-dust", "no", "-evalue", "1e-20",
                        "-query", sfa, "-subject", sfa, "-outfmt", "6 qstart qend sstart send length"],
                       capture_output=True, text=True)
    off = collections.Counter()
    n = 0
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        qs, qe, ss, se, aln = int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4])
        if se < ss or aln < MINALN:
            continue
        d = ss - qs
        if d <= 0 or d > 30000:
            continue
        n += 1
        off[int(round(d / float(BIN))) * BIN] += aln
    hist = ",".join("%d:%d" % (k, off[k]) for k in sorted(off) if off[k] >= 300)
    w.write(TAB.join([aid, str(slen), str(n), hist]) + NL)
    w.flush()
    print("%-52s HSP %6d · 옵셋칸 %4d" % (aid, n, len(hist.split(",")) if hist else 0))
    sys.stdout.flush()
w.close()
print("표 : ~/lit5s/OFFHIST.tsv")
