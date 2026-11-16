# -*- coding: utf-8 -*-
"""FiberHMM BED12 -> 5S unit coordinates

The FiberHMM output run on the placement BAMs (pbmeth_v12d/bam/*.G8v3.bam) is
cut at the unit boundaries of MAP_V11 and moved to rel coordinates from the TSS.

No realignment: the coordinates are those given by the placement step.
The TSS follows the rule of a_lib : the last occurrence of GTCTACGGC in the unit sequence.
"""
import io
import os
import sys
import glob
import collections

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
W = H("~/fhmm")
MAPF = H("~/MAP_V11.tsv")
ANCH = "GTCTACGGC"


def rd_fa(p):
    d = {}
    nm = None
    buf = []
    for l in io.open(p, encoding="utf8"):
        if l.startswith(">"):
            if nm:
                d[nm] = "".join(buf).upper()
            nm = l[1:].strip().split()[0]
            buf = []
        else:
            buf.append(l.strip())
    if nm:
        d[nm] = "".join(buf).upper()
    return d


# -- the unit boundaries
U = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TB)))
for L in f:
    r = L.rstrip(NL).split(TB)
    if r[ix["seg_type"]] != "UNIT":
        continue
    U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]),
                            r[ix["array_pos"]]))
for k in U:
    U[k].sort()

out = io.open(os.path.join(W, "FHMM_UNITS.tsv"), "w", encoding="utf8", newline=NL)
out.write(TB.join(["sample", "slot", "array_pos", "read_id", "cls",
                   "rel_start", "rel_end", "size"]) + NL)
cov = io.open(os.path.join(W, "FHMM_READCOV.tsv"), "w", encoding="utf8", newline=NL)
cov.write(TB.join(["sample", "slot", "array_pos", "read_id",
                   "rel_lo", "rel_hi"]) + NL)

TOT = collections.Counter()
for refp in sorted(glob.glob(os.path.join(H("~/pbmeth_v12d/bam/ref"), "*.fa"))
                   + glob.glob(os.path.join(H("~/fsmeth_v12c/bam/ref"), "*.fa"))):
    smp = os.path.basename(refp).split(".")[0]
    exd = os.path.join(W, "out", "ex_" + smp)
    if not os.path.isdir(exd):
        continue
    SEQ = rd_fa(refp)
    # the TSS per unit, 1-based in reference coordinates
    UT = {}
    for rn, s in SEQ.items():
        sm_, sl = rn.rsplit("_h", 1)
        for (us, ue, ap) in U.get((sm_, sl), []):
            u = s[us - 1:ue]
            p = u.rfind(ANCH)
            if p < 0:
                p = len(u) - 119
            UT.setdefault(rn, []).append((us, ue, ap, us + p))
    for rn in UT:
        UT[rn].sort()

    for cls, suf in (("tf", "_tf.bed"), ("nuc", "_nucleosome.bed"), ("msp", "_msp.bed")):
        bp = os.path.join(exd, smp + ".calls" + suf)
        if not os.path.exists(bp):
            continue
        for L in io.open(bp, encoding="utf8"):
            x = L.rstrip(NL).split(TB)
            if len(x) < 12:
                continue
            rn, cs, rid = x[0], int(x[1]), x[3]
            if rn not in UT:
                continue
            sm_, sl = rn.rsplit("_h", 1)
            bl = [int(v) for v in x[10].rstrip(",").split(",") if v != ""]
            bs = [int(v) for v in x[11].rstrip(",").split(",") if v != ""]
            if cls == "tf":
                # the range a read covers is also recorded once only
                r0, r1 = int(x[1]) + 1, int(x[2])
                for (us, ue, ap, tss) in UT[rn]:
                    a, b = max(r0, us), min(r1, ue)
                    if b - a >= 200:
                        cov.write(TB.join([sm_, sl, ap, rid,
                                           "%d" % (a - tss), "%d" % (b - tss)]) + NL)
            for sz, st in zip(bl, bs):
                g0 = cs + st + 1            # the 1-based reference start
                g1 = g0 + sz - 1
                for (us, ue, ap, tss) in UT[rn]:
                    if g1 < us or g0 > ue:
                        continue
                    out.write(TB.join([sm_, sl, ap, rid, cls,
                                       "%d" % (g0 - tss), "%d" % (g1 - tss),
                                       "%d" % sz]) + NL)
                    TOT[(smp, cls)] += 1
out.close()
cov.close()
for k in sorted(TOT):
    print("%-12s %-5s %8d" % (k[0], k[1], TOT[k]))
print("->", os.path.join(W, "FHMM_UNITS.tsv"))
