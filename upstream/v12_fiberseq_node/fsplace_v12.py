# -*- coding: utf-8 -*-
# Fiber-seq read placement: the reads are placed on the hifiasm consensus map (CONSENSUS_MAP_G6) by a chain score.
#
#  Why chain placement:
#    among the Fiber-seq reads only 3.6% span the array with at least 300 bp of flank on both sides.
#    Chain placement places a read with no flank as well, so the 31.6% with an array span of at least 12 units can be used.
#
#  The gates are the ones already validated on ordinary ONT data:
#    K=12 · Δ>=10 · MINARR=26,880
#    on a pairing test over 1,067 identical reads, K = 12 gave 97.52% against 97.63% for K = 16,
#    and of the 957 passing both at a margin of at least 10, none was called differently; the null was 0.00% for both.
#
#  MAXFL is left unlimited by default:
#    a Fiber-seq read is absent from HAP_READ_FINAL and so has no junction placement to collide with.
#    With nothing to overlap, there is no reason to discard a read carrying flank.
#
import io
import os
import re
import sys
import gzip
import glob
import random
import subprocess
import collections

W = os.path.expanduser(os.environ["W"])
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CMPD = os.path.expanduser(os.environ["CMPD"])
NL = "\n"
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
FL = 5000
K = int(os.environ.get("K", "16"))
DMIN = int(os.environ.get("DMIN", "10"))
MINARR = int(os.environ.get("MINARR", "36000"))
MAXFL = int(os.environ.get("MAXFL", "999999999"))
CANON = os.path.join(R, "unit_canon.fa")
C = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
random.seed(20260909)


def rc(s):
    return s.translate(C)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


FP = rd(os.path.join(R, "FLANK_proximal_50000.fa"))
FD = rd(os.path.join(R, "FLANK_distal_50000.fa"))
SEQF = {}
HSNF = os.path.expanduser(os.environ["HSNF"])
for x in io.open(HSNF, encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
SUMF = os.path.expanduser(os.environ["SUMF"])
for i, x in enumerate(io.open(SUMF, encoding="utf8")):
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
UB = collections.defaultdict(list)
MAPV = collections.defaultdict(dict)
MAPF = os.path.expanduser(os.environ["MAPF"])
f = io.open(MAPF, encoding="utf8")
h = f.readline().rstrip(NL).split("\t")
ix = dict((n, i) for i, n in enumerate(h))
for L in f:
    r = L.rstrip(NL).split("\t")
    if len(r) < 15 or r[ix["seg_type"]].startswith("PARTIAL"):
        continue
    k = (r[0], r[1])
    u = int(r[ix["array_pos"]])
    UB[k].append((u, int(r[ix["unit_start"]]) - 1, int(r[ix["unit_end"]]), int(r[ix["N"]])))
    MAPV[k][u] = set() if r[ix["varlist"]] == "-" else set(r[ix["varlist"]].split(","))
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r"([:*+-])([A-Za-z0-9]+)")


def tuc_targets(smp):
    """Choose the reads that hold enough of the array from the Fiber-seq call table.

    Difference from the ordinary ONT version: only the location of the call table
      ordinary ONT : ~/{s5s_ont_tuc,new44_tuc,kgp_tuc,giab_tuc,s5s_tuc}/out/<sample>.bound.tsv
      Fiber-seq    : ~/fs41/out/<sample>.bound.tsv   (same columns)

    No flank limit is applied (MAXFL unlimited by default)
      Fiber-seq reads have no junction placement: HAP_READ_FINAL is built from ordinary ONT reads
      only, so Fiber-seq reads are not in that table. With nothing to overlap,
      there is no reason to drop reads with flank. Measured: spanning reads (flank >=300 bp on both sides) are only 3.6%,
      so a flank limit would lose most of them.
    """
    out = set()
    p = os.path.expanduser("~/fs41/out/%s.bound.tsv" % smp)
    if not os.path.exists(p):
        return out
    f = io.open(p, encoding="utf8")
    hh = f.readline().rstrip(NL).split("\t")
    jx = dict((n, i) for i, n in enumerate(hh))
    if "arr_span" not in jx or "read" not in jx:
        sys.stderr.write("  열 없음 %s : %s%s" % (smp, ",".join(hh), NL))
        return out
    for L in f:
        r = L.rstrip(NL).split("\t")
        if len(r) < len(hh):
            continue
        try:
            if (int(r[jx["arr_span"]]) >= MINARR
                    and int(r[jx["flankL"]]) < MAXFL
                    and int(r[jx["flankR"]]) < MAXFL):
                out.add(r[jx["read"]])
        except ValueError:
            pass
    return out


def varof(seqs, tag):
    q = os.path.join(W, "u_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm5", "-c", "--cs", "-t", TH,
                        "--secondary=no", CANON, q], capture_output=True, text=True)
    V = {}
    for L in r.stdout.splitlines():
        f = L.split("\t")
        if len(f) < 12:
            continue
        i = int(f[0][1:])
        if i in V:
            continue
        cs = ""
        for x in f[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        tp = int(f[7])
        out = set()
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == ":":
                tp += int(v)
            elif op == "*":
                out.add("%d%s" % (tp + 1, v[1].upper()))
                tp += 1
            elif op == "-":
                out.add("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                out.add("%dI%d" % (tp + 1, len(v)))
        V[i] = out
    os.remove(q)
    return V


def sc(S, V):
    return len(S & V) - len(S - V)


def place(chain, M):
    keys = sorted(M)
    best = []
    for i0 in keys:
        if not all((i0 + t) in M for t in range(len(chain))):
            continue
        best.append((sum(sc(M[i0 + t], chain[t]) for t in range(len(chain))), i0))
    if len(best) < 2:
        return None
    best.sort(reverse=True)
    return best[0][1], best[0][0] - best[1][0]


DONE = set()
_dp = os.path.join(W, "done.txt")
if os.path.exists(_dp):
    DONE = set(x.strip() for x in io.open(_dp, encoding="utf8") if x.strip())
OUT = io.open(os.path.join(W, "fsjx_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
# the per-unit variant list, the material for the read consensus map; positions alone cannot be compared as sequence
VOU = io.open(os.path.join(W, "fsjv_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
# "others" is rebuilt for each individual inside the loop below; a module-level cache would be wrong here.
for smp in [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]:
    if smp in DONE:
        continue
    want = tuc_targets(smp)
    slots = sorted(s for (a, s) in UB if a == smp)
    if not want or len(slots) < 2:
        continue
    # the Fiber-seq reads are uncompressed .hit.fq files, with the MM/ML methylation tags on the header line.
    #   The read id is the first whitespace-separated token, so it reads correctly even with the tags present.
    p = os.path.expanduser("~/fs41/out/%s.hit.fq" % smp)
    if not (os.path.exists(p) and os.path.getsize(p) > 1000):
        continue
    ref = os.path.join(W, "%s.ref.fa" % smp)
    ARR = {}
    with io.open(ref, "w", newline=NL) as g:
        for sl in slots:
            q = os.path.join(CMPD, SEQF.get((smp, sl), "_"))
            if not os.path.exists(q):
                continue
            a = rd(q)
            if RCF.get((smp, sl)):
                a = rc(a)
            ARR[sl] = a
            g.write(">%s_h%s%s%s%s" % (smp, sl, NL, rc(FD[:FL]) + a + rc(FP[-FL:]), NL))
    if len(ARR) < 2:
        continue
    sub = os.path.join(W, "%s.fa" % smp)
    nw = 0
    with io.open(sub, "w", newline=NL) as g:
        with io.open(p, encoding="utf8") as fh:
            j = 0
            rid = None
            for L in fh:
                m = j & 3
                j += 1
                if m == 0:
                    rid = L[1:].split()[0]
                elif m == 1 and rid in want:
                    g.write(">%s%s%s%s" % (rid, NL, L.strip(), NL))
                    nw += 1
    if nw == 0:
        for q in (ref, sub):
            try:
                os.remove(q)
            except OSError:
                pass
        continue
    pr = subprocess.Popen(["minimap2", "-a", "-x", "map-ont", "-t", TH, "-K", "100M",
                           "--secondary=no", ref, sub],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    REC = []
    for L in pr.stdout:
        if L.startswith("@"):
            continue
        f = L.split("\t")
        if len(f) < 10 or f[2] == "*":
            continue
        if int(f[1]) & 0x900 or int(f[4]) < 20:
            continue
        ctg = f[2].rsplit("_h", 1)[1]
        a = ARR.get(ctg)
        if a is None:
            continue
        seq = f[9]
        rpos = int(f[3]) - 1
        q2 = 0
        m2 = {}
        for m in CIG.finditer(f[5]):
            k2, op = int(m.group(1)), m.group(2)
            if op == "S":
                q2 += k2
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k2):
                    m2[rpos + t] = q2 + t
                rpos += k2
                q2 += k2
            elif op == "I":
                q2 += k2
            elif op in "DN":
                rpos += k2
        if not m2:
            continue
        seg = []
        for u, s0, e0, N in UB[(smp, ctg)]:
            A, B = s0 + FL, e0 + FL - 1
            if A in m2 and B in m2:
                seg.append((u, seq[m2[A]:m2[B] + 1]))
        if len(seg) >= K:
            REC.append((f[0], ctg, seg, int(f[4])))
    pr.stdout.close()
    pr.wait()
    for q in (ref, sub):
        try:
            os.remove(q)
        except OSError:
            pass
    if not REC:
        io.open(_dp, "a", encoding="utf8", newline=NL).write(smp + NL)
        continue
    # recomputed per individual
    others = sorted(set(a for (a, s) in UB if a != smp))
    flat = [s for _r, _c, seg, _q in REC for _u, s in seg]
    V = varof(flat, "%s_%s" % (TAG, smp))
    j = 0
    npass = nfail = 0
    for rid, ctg, seg, mq in REC:
        ch = []
        for t in range(len(seg)):
            ch.append((seg[t][0], V.get(j)))
            j += 1
        if any(x[1] is None for x in ch):
            continue
        vs = [x[1] for x in ch[:K]]
        r0 = place(vs, MAPV[(smp, ctg)])
        if not r0:
            continue
        # the null control is applied at placement time rather than afterwards
        o = random.Random(rid).choice(others)
        os_ = sorted(s for (a, s) in UB if a == o)
        nm = 0
        if os_:
            M2 = MAPV.get((o, os_[0]))
            if M2:
                r1 = place(vs, M2)
                if r1:
                    nm = r1[1]
        dlt = r0[1] - nm
        if dlt < DMIN:
            nfail += 1
            continue
        npass += 1
        # placement succeeded: the units this read covered are recorded
        N = UB[(smp, ctg)][0][3]
        # the unit record is not truncated at K.
        #   The position was already fixed by the first K units, and the remaining units being consecutive,
        #   their positions follow automatically as r0[0] + t.
        #   A read spans a median of 23.6 units; recording only 16 would lose 35.9%.
        #   The units beyond the Kth did not contribute to the placement score, however,
        #     so they are marked in the scored column, to keep them apart.
        for t, (u, _s) in enumerate(seg):
            # scored = 1 marks the first K units, which contributed to the placement score, and 0 those beyond
            sc_ = 1 if t < K else 0
            OUT.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d%s"
                      % (smp, rid[:16], ctg, r0[0] + t, N + 1 - (r0[0] + t),
                         r0[1], nm, dlt, mq, sc_, NL))
            # building a read consensus map requires the variant list:
            #   the position (array_pos) alone cannot be used to compare sequence.
            vl = ",".join(sorted(ch[t][1])) if ch[t][1] else "-"
            VOU.write("%s\t%s\t%s\t%d\t%s\t%d%s"
                      % (smp, rid[:16], ctg, r0[0] + t, vl, sc_, NL))
    OUT.flush()
    VOU.flush()
    io.open(_dp, "a", encoding="utf8", newline=NL).write(smp + NL)
    sys.stderr.write("  %-10s TUC대상 %-5d 정렬 %-5d · 채택 %-5d 탈락 %d%s"
                     % (smp, len(want), len(REC), npass, nfail, NL))
    sys.stderr.flush()
OUT.close()
print("done")
