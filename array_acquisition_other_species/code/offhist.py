# -*- coding: utf-8 -*-
# array 조각 자기-정렬 옵셋 히스토그램을 그대로 저장한다 (판정은 따로)
#   앞 판의 결함 : 기본 주기와 그 배수(고차반복) 가운데 정렬량이 큰 쪽을 골라, 500 bp 유닛 array 에
#                 1,520(=3 유닛) 이 주기로 찍혔다. 둘은 따로 내야 한다
#   산출 : OFFHIST.tsv — array_id · 조각길이 · HSP수 · "옵셋:정렬bp" 전부 (옵셋 <= 30000 · 정렬 >= 300)
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
