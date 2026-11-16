# -*- coding: utf-8 -*-
"""FiberHMM BED12 -> 5S unit coordinates, full set

Reads ex_<sample>__<pacbio|nanopore>/<tag>.calls_{tf,nucleosome,msp}.bed,
cuts at the MAP_V11 unit boundaries and moves to rel from the TSS.
No realignment: the coordinates are those given by the placement step (ao_g8cut).
TSS = the last occurrence of GTCTACGGC in the unit sequence (the rule of a_lib)
"""
import io
import os
import glob
import collections

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
W = H("~/fhmm")
MAPF = H("~/MAP_V11.tsv")
ANCH = "GTCTACGGC"
HOST = os.uname()[1]
REFD = [H("~/pbmeth_v12d/bam/ref"), H("~/fsmeth_v12c/bam/ref")]


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

REFP = {}
for d in REFD:
    for p in glob.glob(os.path.join(d, "*.fa")):
        REFP.setdefault(os.path.basename(p).split(".")[0], []).append(p)

out = io.open(os.path.join(W, "FHMM_UNITS_%s.tsv" % HOST), "w", encoding="utf8", newline=NL)
out.write(TB.join(["sample", "plat", "slot", "array_pos", "read_id", "cls",
                   "rel_start", "rel_end", "size"]) + NL)
cov = io.open(os.path.join(W, "FHMM_READCOV_%s.tsv" % HOST), "w", encoding="utf8", newline=NL)
cov.write(TB.join(["sample", "plat", "slot", "array_pos", "read_id",
                   "rel_lo", "rel_hi"]) + NL)

TOT = collections.Counter()
SEEN = set()
for exd in sorted(glob.glob(os.path.join(W, "out", "ex_*"))):
    tag = os.path.basename(exd)[3:]
    if "__" not in tag:
        continue
    smp, plat = tag.split("__", 1)
    paths = REFP.get(smp)
    if not paths:
        continue
    want = [p for p in paths if (("pbmeth" in p) == (plat == "pacbio"))]
    if not want:
        continue
    SEQ = rd_fa(want[0])
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
        bp = os.path.join(exd, tag + ".calls" + suf)
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
            # The span a read covers is written once, whichever class it is in
            #   (writing it only in the tf BED would drop reads with no TF from the denominator)
            r0, r1 = cs + 1, int(x[2])
            for (us, ue, ap, tss) in UT[rn]:
                if min(r1, ue) - max(r0, us) < 200:        # Only reads that cover this unit well enough
                    continue
                if (sm_, plat, sl, ap, rid) in SEEN:
                    continue
                SEEN.add((sm_, plat, sl, ap, rid))
                # The cover is not trimmed to the unit - the array is tandem, so the flank is the real sequence of the neighbouring unit.
                #   Trimming it would cut the denominator at the unit end (rel ~ +118) and raise a false peak
                a, b = max(r0, tss - 1200), min(r1, tss + 1200)
                cov.write(TB.join([sm_, plat, sl, ap, rid,
                                   "%d" % (a - tss), "%d" % (b - tss)]) + NL)
            for sz, st in zip(bl, bs):
                g0 = cs + st + 1
                g1 = g0 + sz - 1
                for (us, ue, ap, tss) in UT[rn]:
                    if g1 < tss - 1200 or g0 > tss + 1200:
                        continue
                    out.write(TB.join([sm_, plat, sl, ap, rid, cls,
                                       "%d" % (g0 - tss), "%d" % (g1 - tss),
                                       "%d" % sz]) + NL)
                    TOT[(plat, cls)] += 1
out.close()
cov.close()
for k in sorted(TOT):
    print("  %-9s %-4s %9d" % (k[0], k[1], TOT[k]))
print("[%s] -> FHMM_UNITS_%s.tsv" % (HOST, HOST))
