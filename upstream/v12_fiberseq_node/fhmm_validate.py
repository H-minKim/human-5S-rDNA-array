# -*- coding: utf-8 -*-
"""Check of locus assignment: there should be no m6A inside a FiberHMM footprint.

The m6A positions of the M6A_READUNIT table (output of the placement step) are taken, and
the m6A density inside the intervals called by FiberHMM (tf, nuc, msp) is compared.
If the coordinate assignment were wrong, the contrast would vanish.

Coordinates : the table is in unit_pos (1-based), the FHMM output in rel.  rel = unit_pos - TSSU,
       TSSU = position (1-based) of the last occurrence of GTCTACGGC in the unit sequence, the rule of a_lib.
"""
import io
import os
import gzip
import glob
import collections

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
W = H("~/fhmm")
MLTHR = 128
ANCH = "GTCTACGGC"
SMPS = {"HG02486", "NA19036"}
MAPF = H("~/MAP_V11.tsv")


def rd_fa(p):
    d, nm, buf = {}, None, []
    for l in io.open(p, encoding="utf8"):
        if l.startswith(">"):
            if nm:
                d[nm] = "".join(buf).upper()
            nm, buf = l[1:].strip().split()[0], []
        else:
            buf.append(l.strip())
    if nm:
        d[nm] = "".join(buf).upper()
    return d


U = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TB)))
for L in f:
    r = L.rstrip(NL).split(TB)
    if r[ix["seg_type"]] == "UNIT":
        U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]),
                                r[ix["array_pos"]]))

TSSU = {}
for refp in glob.glob(H("~/pbmeth_v12d/bam/ref/*.fa")):
    smp = os.path.basename(refp).split(".")[0]
    if smp not in SMPS:
        continue
    for rn, s in rd_fa(refp).items():
        sm_, sl = rn.rsplit("_h", 1)
        for (us, ue, ap) in U.get((sm_, sl), []):
            u = s[us - 1:ue]
            p = u.rfind(ANCH)
            if p < 0:
                p = len(u) - 119
            TSSU[(sm_, sl, ap)] = p + 1

SEG = collections.defaultdict(lambda: collections.defaultdict(list))
f = io.open(os.path.join(W, "FHMM_UNITS.tsv"), encoding="utf8")
f.readline()
for L in f:
    x = L.rstrip(NL).split(TB)
    SEG[(x[0], x[3], x[1], x[2])][x[4]].append((int(x[5]), int(x[6])))

rows = 0
nA = collections.Counter()
nM = collections.Counter()
with gzip.open(H("~/pbmeth_v12d/M6A_READUNIT_ip-172-31-58-58.tsv.gz"),
               "rt", encoding="utf8") as g:
    g.readline()
    for L in g:
        x = L.rstrip(NL).split(TB)
        if len(x) < 14 or x[0] not in SMPS or x[12] != "Y" or not x[13]:
            continue
        k = (x[0], x[6], x[4], x[5])
        s = SEG.get(k)
        tu = TSSU.get((x[0], x[4], x[5]))
        if not s or tu is None:
            continue
        rows += 1
        for t in x[13].split(","):
            if not t:
                continue
            a_, b_ = t.split(":")
            rel = int(a_) - tu
            h = int(b_) >= MLTHR
            for cls in ("tf", "nuc", "msp"):
                for (p0, p1) in s.get(cls, []):
                    if p0 <= rel <= p1:
                        nA[cls] += 1
                        nM[cls] += h
                        break

print("matched read x unit (strand rows) %d" % rows)
print("  %-5s %10s %10s %8s" % ("region", "A 자리", "m6A", "비율"))
for cls in ("tf", "nuc", "msp"):
    a, m = nA[cls], nM[cls]
    print("  %-5s %10d %10d %8.4f" % (cls, a, m, (m / a) if a else 0))
if nA["msp"] and nA["nuc"]:
    print()
    print("  MSP / nuc ratio = %.2f-fold (should be large if the assignment is right)"
          % ((nM["msp"] / nA["msp"]) / max(nM["nuc"] / nA["nuc"], 1e-9)))
    print("  MSP / tf  ratio = %.2f-fold"
          % ((nM["msp"] / nA["msp"]) / max(nM["tf"] / max(nA["tf"], 1), 1e-9)))
