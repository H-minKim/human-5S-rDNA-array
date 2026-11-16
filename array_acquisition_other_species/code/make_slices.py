# -*- coding: utf-8 -*-
# Cut a fragment of up to 150 kb from the middle of the array, as the input for confirming the tandem repeat directly on the node by gapped alignment.
#   The header records the original array coordinates and the fragment's start offset, so that it can be traced back.
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import glob

NL = chr(10)
CAN = "$RDNA_OTHERSP_DIR/data_v2/01_canonical"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "slices")
WIN = 150000
os.makedirs(OUT, exist_ok=True)
n = 0
tot = 0
for p0 in sorted(glob.glob(os.path.join(CAN, "*", "*.fa"))):
    if p0.endswith(".unit.fa"):
        continue
    aid = os.path.basename(p0).replace(".fa", "")
    L = io.open(p0, encoding="utf8")
    hd = L.readline().strip()
    seq = "".join(x.strip() for x in L).upper()
    kv = dict(x.split("=", 1) for x in hd.split() if "=" in x)
    a, b = kv["array"].split("-")
    a0, a1 = int(a), int(b)
    span = a1 - a0 + 1
    fl = (len(seq) - span) // 2
    arr = seq[fl:fl + span]
    st = 0 if len(arr) <= WIN else (len(arr) - WIN) // 2
    sub = arr[st:st + WIN]
    with io.open(os.path.join(OUT, aid + ".arr.fa"), "w", encoding="utf8", newline=NL) as g:
        g.write(">%s array=%d-%d span=%d slice_off=%d slice_len=%d%s" % (aid, a0, a1, span, st, len(sub), NL))
        for i in range(0, len(sub), 80):
            g.write(sub[i:i + 80] + NL)
    n += 1
    tot += len(sub)
print("조각 %d 개 · 합 %.1f Mb" % (n, tot / 1e6))
