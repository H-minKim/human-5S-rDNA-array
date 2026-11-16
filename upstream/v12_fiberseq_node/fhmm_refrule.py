# -*- coding: utf-8 -*-
"""The rule of the reference Figure 5 applied to the FiberHMM calls.

Tullius 2024 Data S4 cell 47 (5S), as it is
  find_fp_starts(track, 20, 180, -55+250, -45+250)
  = a boundary (start or end) of a footprint of size 20-180 at diff index 196..204

track has the shape of the reference nc_df :
  one row = read x unit, columns = rel -250..+249
  value = size of the footprint covering that position; 2 = the read does not cover it; 0 = covered but no footprint
footprint = tf + nuc of FiberHMM (MSP is an open interval and is excluded)
"""
import io
import os
import collections
import numpy as np

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
W = H("~/fhmm")
TLO, THI = -250, 249
NT = THI - TLO + 1
NOCOV = 2
R_SMIN, R_SMAX = 20, 180
R_MINP, R_MAXP = -55 + 250, -45 + 250

SEG = collections.defaultdict(list)
f = io.open(os.path.join(W, "FHMM_UNITS.tsv"), encoding="utf8")
f.readline()
for L in f:
    x = L.rstrip(NL).split(TB)
    if x[4] in ("tf", "nuc"):
        SEG[(x[0], x[3], x[1], x[2])].append((int(x[5]), int(x[6]), int(x[7]), x[4]))

COV = {}
f = io.open(os.path.join(W, "FHMM_READCOV.tsv"), encoding="utf8")
f.readline()
for L in f:
    x = L.rstrip(NL).split(TB)
    COV[(x[0], x[3], x[1], x[2])] = (int(x[4]), int(x[5]))

HTF = np.zeros(NT)          # the distribution of TF footprint centres
HFP = np.zeros(NT)          # a cover by any footprint at all
NCO = np.zeros(NT)
nref = 0
ntot = 0
szs = []
for k, (lo, hi) in COV.items():
    t = np.full(NT, NOCOV, np.int32)
    a0, a1 = max(TLO, lo), min(THI, hi)
    if a1 < a0:
        continue
    t[a0 - TLO:a1 - TLO + 1] = 0
    for (s0, s1, sz, cls) in SEG.get(k, []):
        b0, b1 = max(TLO, s0), min(THI, s1)
        if b1 >= b0:
            t[b0 - TLO:b1 - TLO + 1] = sz
        if cls == "tf" and TLO <= (s0 + s1) // 2 <= THI:
            HTF[(s0 + s1) // 2 - TLO] += 1
            szs.append(sz)
    ntot += 1
    NCO[t != NOCOV] += 1
    HFP[(t != NOCOV) & (t != 0)] += 1
    c = t.copy()
    c[(c < R_SMIN) | (c > R_SMAX)] = 0
    b = c[1:] - c[:-1]
    idx = np.nonzero(b)[0]
    if np.any((idx > R_MINP) & (idx < R_MAXP)):
        nref += 1

print("read x 유닛 %d · 레퍼런스 규칙 양성 %d (%.2f%%)"
      % (ntot, nref, 100.0 * nref / max(ntot, 1)))
szs = np.array(szs)
if len(szs):
    print("TF footprint %d · 크기 중앙 %.0f (25~75%% %.0f~%.0f)"
          % (len(szs), np.median(szs), *np.percentile(szs, [25, 75])))
print()
print("== TF footprint 중심 분포 (20 bp 칸 · TSS 기준) ==")
for a in range(-200, 200, 20):
    i0, i1 = a - TLO, a + 20 - TLO
    n = HTF[i0:i1].sum()
    c = NCO[i0:i1].mean()
    print("  %+5d..%+5d  TF %6.0f   덮개 %7.0f   섬유당 %.4f"
          % (a, a + 19, n, c, n / max(c, 1)))
