# -*- coding: utf-8 -*-
# PacBio Fiber-seq read placement: the ONT version, pbontplace.py, carried across unchanged.
#
#  What differs from pbontplace.py is only what the platform requires:
#    1) the preset is minimap2 -x map-hifi, against map-ont for ONT;
#    2) the gates are K = 4, a margin of at least 5, and MINARR = 8,960, that is 4 units,
#               the values validated on HiFi data (hifitest and hifigate over 211 individuals and 3,824 reads:
#                 99.08% correct, 0.00% where another individual's map scored higher, and no call differing between K = 4 and K = 6, 8 or 10).
#               Applying the ONT gates, K = 12 and MINARR 26,880, to PacBio would drop most of the reads
#                 (a PacBio Fiber-seq read has a median length of about 17 kb, roughly 7 units).
#    3) the whole read name is recorded, PacBio names being movie/zmw/ccs so that the first 16 characters coincide;
#    4) the outputs are pbfsplace/pfjx_<TAG>.tsv and pfjv_<TAG>.tsv,
#               pbu_<TAG>.tsv, the 17-column AO table for the recall, including the read coordinates f0 and f1 of each unit,
#               and <SAMD>/<sample>.<PBJOB>.sam, holding those reads with flag 0 and the SEQ in the reference orientation.
#               They are then recalled by recall_pad2.py AO pbu_*.tsv out.tsv at PAD = 100 and CTX = 200.
#                 The variant tokens of pfjv are only a fingerprint for the placement, taken against the bare canonical unit with asm5 and no context.
#                 The variants entering the merge must always come from the recall.
#    5) the map paths come from environment variables, defaulting to the canonical ones, so that omitting them leaves the behaviour unchanged:
#         MAPF  CONSENSUS_MAP (for hg002, slot 2 must be given the map corrected to the T2T paternal assembly)
#         SUMF  CONSENSUS_SUM (revcomp is located by column name: the 6th in G6, and
#                                the 5th in G8v4, which has no seq_file column)
#         HSNF  hap_seq_node, and CMPD the folder of array .fa files
#  Inputs: pbfs/out/<label>.hit.fq and <label>.bound.tsv, where the label is <map individual>__<data source>.
#  Everything else is the same: the units are cut on the alignment coordinates, the variants called with asm5, the score is |S intersect V| - |S - V|, the null runs at placement time, and every unit is recorded with its scored flag.
import io
import os
import re
import sys
import random
import subprocess
import collections

W = os.path.expanduser(os.environ["W"])
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
D = os.path.expanduser(os.environ["READD"])
NL = "\n"
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
FL = 5000
K = int(os.environ.get("K", "4"))
DMIN = int(os.environ.get("DMIN", "5"))
MINARR = int(os.environ.get("MINARR", "8960"))
MAXFL = int(os.environ.get("MAXFL", "999999999"))
PRESET = "map-hifi"
MAPF = os.path.expanduser(os.environ["MAPF"])
SUMF = os.path.expanduser(os.environ["SUMF"])
HSNF = os.path.expanduser(os.environ["HSNF"])
CMPD = os.path.expanduser(os.environ["CMPD"])
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
for x in io.open(HSNF, encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
RCI = None
for i, x in enumerate(io.open(SUMF, encoding="utf8")):
    f = x.rstrip(NL).split("\t")
    if i == 0:
        # revcomp is located by column name: it is the 6th column of the G6 SUM but in the G8v4 SUM
        #   it is the 5th, having no seq_file column. Reading it by index would read off_prox as revcomp
        #   and silently treat every array as forward.
        if "revcomp" not in f:
            raise SystemExit("SUMF has no revcomp column: %s" % SUMF)
        RCI = f.index("revcomp")
        continue
    if len(f) > RCI:
        RCF[(f[0], f[1])] = (f[RCI] == "Y")
UB = collections.defaultdict(list)
MAPV = collections.defaultdict(dict)
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


def tuc_targets(lab):
    out = set()
    p = os.path.join(D, "%s.bound.tsv" % lab)
    if not os.path.exists(p):
        return out
    f = io.open(p, encoding="utf8")
    hh = f.readline().rstrip(NL).split("\t")
    jx = dict((n, i) for i, n in enumerate(hh))
    if "arr_span" not in jx or "read" not in jx:
        sys.stderr.write("  열 없음 %s : %s%s" % (lab, ",".join(hh), NL))
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
OUT = io.open(os.path.join(W, "pfjx_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
VOU = io.open(os.path.join(W, "pfjv_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
# the output for the recall (recall_pad2.py): the table of read coordinates per unit, plus the SAM of those reads
UOU = io.open(os.path.join(W, "pbu_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
PBJOB = os.environ.get("PBJOB", "PBFS")
DSRC = os.environ.get("DSRC", "PacBio_FiberSeq")
SAMD = os.path.expanduser(os.environ["SAMD"])
os.makedirs(SAMD, exist_ok=True)
SAMR = collections.defaultdict(list)
# "others" is rebuilt for each individual inside the loop below.
for lab in [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]:
    if lab in DONE:
        continue
    smp = lab.split("__")[0]
    want = tuc_targets(lab)
    slots = sorted(s for (a, s) in UB if a == smp)
    p = os.path.join(D, "%s.hit.fq" % lab)
    if not want or len(slots) < 2 or not (os.path.exists(p) and os.path.getsize(p) > 1000):
        sys.stderr.write("  %-18s 건너뜀 (TUC대상 %d · 지도 %s 슬롯 %d · read %s)%s"
                         % (lab, len(want), smp, len(slots), os.path.exists(p), NL))
        continue
    ref = os.path.join(W, "%s.ref.fa" % lab)
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
        sys.stderr.write("  %-18s 건너뜀 (array 파일 %d 개)%s" % (lab, len(ARR), NL))
        continue
    sub = os.path.join(W, "%s.fa" % lab)
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
    pr = subprocess.Popen(["minimap2", "-a", "-x", PRESET, "-t", TH, "-K", "100M",
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
                # f0 and f1 are in SAM SEQ coordinates, that is the reference orientation, and f1 is exclusive, equal to the f0 of the next unit.
                #   These coordinates are passed to recall_pad2.py unchanged and the SAM is written with flag 0,
                #   so its reverse branch is never taken and an orientation flip cannot arise in the first place.
                seg.append((u, seq[m2[A]:m2[B] + 1], m2[A], m2[B] + 1))
        if len(seg) >= K:
            REC.append((f[0], ctg, seg, int(f[4]), seq))
    pr.stdout.close()
    pr.wait()
    for q in (ref, sub):
        try:
            os.remove(q)
        except OSError:
            pass
    if not REC:
        io.open(_dp, "a", encoding="utf8", newline=NL).write(lab + NL)
        continue
    # The competitor candidates are rebuilt for every individual.
    others = sorted(set(a for (a, s) in UB if a != smp))
    flat = [s for _r, _c, seg, _q, _sq in REC for _u, s, _a, _b in seg]
    V = varof(flat, "%s_%s" % (TAG, lab))
    j = 0
    npass = nfail = 0
    for rid, ctg, seg, mq, rseq in REC:
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
        # competitors are selected deterministically by read name.
        #   The test is still a comparison against one competitor, so the K = 4 and margin >= 5 calibration holds.
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
        N = UB[(smp, ctg)][0][3]
        SAMR[smp].append((rid, rseq))
        for t, (u, _s, f0_, f1_) in enumerate(seg):
            sc_ = 1 if t < K else 0
            # the 17-column AO table for the recall, in the same shape as the input to recall_pad2.py
            #   1 individual, 2 full read name, 3 short read name, 4 slot, 5 job, 6 source, 7 unused, 8 array_pos
            #   9 array_pos 10 f0 11 f1 12 - 13 - 14 verify 15 varlist 16 - 17 Y
            #   verify and varlist are overwritten by the recall; only the positions are set here.
            UOU.write("\t".join([
                smp, rid, rid[:16], ctg, PBJOB, DSRC, "-",
                str(r0[0] + t), str(r0[0] + t), str(f0_), str(f1_),
                "RUN", "RUN", "RUN", "-", "-", "Y"]) + NL)
            OUT.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s%s"
                      % (smp, rid, ctg, r0[0] + t, N + 1 - (r0[0] + t),
                         r0[1], nm, dlt, mq, sc_, lab, NL))
            vl = ",".join(sorted(ch[t][1])) if ch[t][1] else "-"
            VOU.write("%s\t%s\t%s\t%d\t%s\t%d\t%s%s"
                      % (smp, rid, ctg, r0[0] + t, vl, sc_, lab, NL))
    OUT.flush()
    VOU.flush()
    UOU.flush()
    # the SAM for recall_pad2.py to read: flag 0, with the SEQ exactly as used for the placement, in the reference orientation.
    #   f0 and f1 are in the same coordinate frame, so nothing needs reversing at the other end.
    if SAMR[smp]:
        with io.open(os.path.join(SAMD, "%s.%s.sam" % (smp, PBJOB)),
                     "w", encoding="utf8", newline=NL) as sg:
            for rid_, sq_ in SAMR[smp]:
                sg.write("\t".join([rid_, "0", "*", "0", "0", "*", "*", "0", "0",
                                    sq_, "*"]) + NL)
        del SAMR[smp]
    io.open(_dp, "a", encoding="utf8", newline=NL).write(lab + NL)
    sys.stderr.write("  %-18s TUC대상 %-5d 정렬 %-5d · 채택 %-5d 탈락 %d%s"
                     % (lab, len(want), len(REC), npass, nfail, NL))
    sys.stderr.flush()
OUT.close()
VOU.close()
print("done")
