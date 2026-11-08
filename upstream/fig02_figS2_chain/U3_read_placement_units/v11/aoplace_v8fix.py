# -*- coding: utf-8 -*-
# ARRAY_ONLY actual placement, using only the gates set by validation (2026-09-09)
#
#  Gates (fixed in validation, no post hoc adjustment)
#      K = 16  ·  Δ >= 10
#      delta = (margin on own individual's map) - (margin on another individual's map)
#    validation measured: retained 95.3%, correct 98.08%, other's map better 0.00%
#
#  Why delta (not the absolute margin)
#    With the absolute margin the null did not go down (it went up: 8.2% -> 16.2%).
#    Unit numbers run 1~100, so the share of the same number appearing by chance was large.
#    The median margin on another's map is 1 (own map 25),
#    so delta, which asks "does it attach confidently only to the own individual", is the real discriminating axis.
#
#  Units are cut by alignment coordinates (no blastn)
#    The cause that broke chain placement was blastn cutting (unit count agreement 70.8%).
#
#  Target = reads that TUC judged to have caught the array and that have no flank
#    arr_span >= 36000 (16 units = K met), flankL,flankR < 500
import io
import os
import re
import sys
import gzip
import glob
import random
import subprocess
import collections

W = os.path.expanduser(os.environ.get("WDIR", "~/aojx8"))
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
NL = "\n"
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
FL = 5000
K = int(os.environ.get("K", "16"))
DMIN = int(os.environ.get("DMIN", "10"))
MINARR = int(os.environ.get("MINARR", "36000"))
MAXFL = int(os.environ.get("MAXFL", "500"))
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
for x in io.open(os.path.expanduser("~/hap_seq_node.tsv"), encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
for i, x in enumerate(io.open(os.path.expanduser("~/CONSENSUS_SUM_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
SEQO = {}
if os.environ.get("SEQOVR"):
    for x in io.open(os.path.expanduser(os.environ["SEQOVR"]), encoding="utf8"):
        f = x.split()
        if len(f) >= 3:
            SEQO[(f[0], f[1])] = f[2]
UB = collections.defaultdict(list)
MAPV = collections.defaultdict(dict)
f = io.open(os.path.expanduser(os.environ.get("MAPF", "~/CONSENSUS_MAP_G6.tsv")), encoding="utf8")
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
    """Selects reads that TUC judged to have no flank and to have caught enough of the array"""
    out = set()
    # giab_tuc was missing, so the 4 GIAB individuals were silently skipped (fixed).
    #   If even one TUC call table directory is left out, tuc_targets returns an empty set
    #   and they disappear without error via "if not want: continue". Always check the directory list.
    for d in ("s5s_ont_tuc", "new44_tuc", "kgp_tuc", "giab_tuc", "s5s_tuc"):
        p = os.path.expanduser("~/%s/out/%s.bound.tsv" % (d, smp))
        if not os.path.exists(p):
            continue
        f = io.open(p, encoding="utf8")
        hh = f.readline().rstrip(NL).split("\t")
        jx = dict((n, i) for i, n in enumerate(hh))
        for L in f:
            r = L.rstrip(NL).split("\t")
            if len(r) < 13:
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
OUT = io.open(os.path.join(W, "aojx_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
# per-unit variant list, material for the read-consensus map (positions alone do not allow sequence comparison)
VOU = io.open(os.path.join(W, "aojv_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
others = None
for smp in [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]:
    if smp in DONE:
        continue
    want = tuc_targets(smp)
    slots = sorted(s for (a, s) in UB if a == smp)
    if not want or len(slots) < 2:
        continue
    p = None
    for _p in ("~/s5s_ont/hits/%s.hits.fastq.gz" % smp,
               "~/s5s_new/hits/%s.hits.fastq.gz" % smp,
               "~/kgp_ont/out/%s.array.fastq.gz" % smp,
               "~/cancer5s/hits/%s.hits.fastq.gz" % smp,
               "~/hg008/out/%s.hits.fastq.gz" % smp):
        _q = os.path.expanduser(_p)
        if os.path.exists(_q) and os.path.getsize(_q) > 1000:
            p = _q
            break
    if not p:
        # GIAB and HG008 have fastq split into pieces (hg001_1..hg001_4).
        #   The handling of rjfreq2.py is copied as is, concatenated and used as one.
        ch = sorted(glob.glob(os.path.expanduser("~/giab5s_seq/out/%s_*.hits.fastq.gz" % smp))
                    + glob.glob(os.path.expanduser("~/hg008/out/%s_*.hits.fastq.gz" % smp)))
        if ch:
            p = os.path.join(W, "%s.cat.fastq.gz" % smp)
            if not (os.path.exists(p) and os.path.getsize(p) > 1000):
                with io.open(p, "wb") as g:
                    for x in ch:
                        g.write(io.open(x, "rb").read())
    if not p:
        continue
    ref = os.path.join(W, "%s.ref.fa" % smp)
    ARR = {}
    with io.open(ref, "w", newline=NL) as g:
        for sl in slots:
            # slots whose representative sequence was changed by sequence consensus 2 (SEQOVR). Already in the RCF-flipped frame, so not flipped again.
            if (smp, sl) in SEQO:
                a = rd(SEQO[(smp, sl)])
            else:
                q = os.path.expanduser("~/cmp4d/" + SEQF.get((smp, sl), "_"))
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
        with gzip.open(p, "rt") as fh:
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
    others = sorted(set(a for (a, s) in UB if a != smp))   # anew for each individual
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
        # the null control is applied at every placement (in real time, not post hoc)
        o = random.Random(rid).choice(others)   # fixed by read name
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
        # placement succeeded, record the units this read covered
        N = UB[(smp, ctg)][0][3]
        # unit records are not cut at K (fixed)
        #   the position is already fixed by the first K units, and the remaining units are contiguous, so
        #   their positions follow automatically as r0[0]+t. But originally they were discarded.
        #   measured: median 23.6 units actually spanned by a read, but only 16 were written, a 35.9% loss.
        #   However, units after the K-th are an interval not used for the placement score, so
        #     they are marked separately in the scored column (so they are not mixed in before validation).
        for t, (u, _s) in enumerate(seg):
            # scored = 1 means the first K units used for the placement score, 0 means those after
            sc_ = 1 if t < K else 0
            OUT.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d%s"
                      % (smp, rid[:16], ctg, r0[0] + t, N + 1 - (r0[0] + t),
                         r0[1], nm, dlt, mq, sc_, NL))
            # making the read-consensus map needs the variant list.
            #   positions (array_pos) alone do not allow sequence comparison.
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
