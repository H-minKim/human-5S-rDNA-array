# -*- coding: utf-8 -*-
# Reproduction version (2026-09-11): only additions to rjfreq2_g8.py (original rjfreq2 + input environment variables). Decision rules, alignment command and filters unchanged
#   (1) the input fastq is first reduced to the reads of the assignment table (HAP_READ_FINAL). The original also discards unassigned reads after alignment, so the result is the same
#       (minimap2 aligns each read independently)
#   (2) added outputs: readunit_<TAG>.tsv (read x unit, coordinates within the read, both G6 and G8v2)
#                     <individual>.placed.sam (original alignments of accepted reads + MM/ML tags from the fastq)
#                     rjf8_<TAG>.tsv (position table of the same reads and same alignments cut at G8v2 boundaries, the half split uses a separate random stream)
#   (3) G8v2 cutting only for slots whose reference equals the G8v2 canonical (slots in C2SKIP have a different reference and are skipped)
# G8v2 version (2026-09-11): differs from the original only in its inputs (placement rules, thresholds and scores unchanged)
#   Changes: map MAPF (G8v2), orientation table RCFF/RCFCOL (G8v2 SUM revcomp = strand majority), canonical list HSN, canonical folder CMPD, work folder
#   all are required environment variables (stops with KeyError if missing, so it never runs silently on the old G6 inputs)
#   G8v2 cuts units at the end of the gene, so read units are also cut as [spacer][gene] and the gene part (2122-2240) is visible
# junction anchor placement v3: slot label correspondence decided per individual + token frequencies kept
#   
#
#  The defect being fixed
#    There is no guarantee that the slot numbers of the assignment table (HAP_READ_FINAL) and those of the reference assembly (hap_seq_node)
#    point to the same molecule. Measured: in 83 of 303 individuals
#    the two numberings were swapped wholesale and all reads were discarded (e.g. HG00099 0/75).
#    side_hap itself was confirmed to use the existing HPRC split (hap_reads_v5) 100% unchanged.
#    The only problem is that the numbers assigned at the pairing step were never matched to the assembly numbers.
#
#  How it is decided: not left to minimap2 read by read
#    For each individual both orientations are counted and the correspondence is decided on the aggregate.
#      agreement >= 0.90 -> keep, <= 0.10 -> swap, in between -> marked AMBIG
#    The decision rests on unanimous splits such as 0/75 and 88/88. It is not a per-read decision.
#    AMBIG individuals use the majority orientation but a record is kept so they can be filtered later.
#
#  Token frequencies are kept (the threshold is applied later)
#    per position "token:total reads:half A:half B", not fixed by majority vote.
#    At the same threshold the self-reproducibility (A vs B) ceiling can also be measured.
#
#  Position comes from the junction anchor (fixed by flank sequence). The map is not consulted.
import io
import os
import re
import sys
import glob
import gzip
import random
import subprocess
import collections

W = os.path.expanduser(os.environ["RJW"])
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
NL = "\n"
TH = os.environ.get("TH", "3")
TAG = os.environ.get("TAG", "x")
FL = 5000
HI = float(os.environ.get("HI", "0.90"))
LO = float(os.environ.get("LO", "0.10"))
CANON = os.path.join(R, "unit_canon.fa")
C = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
random.seed(20260907)


def rc(s):
    return s.translate(C)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


FP = rd(os.path.join(R, "FLANK_proximal_50000.fa"))
FD = rd(os.path.join(R, "FLANK_distal_50000.fa"))
SEQF = {}
for x in io.open(os.path.expanduser(os.environ["HSN"]), encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
for i, x in enumerate(io.open(os.path.expanduser(os.environ["RCFF"]), encoding="utf8")):
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[int(os.environ["RCFCOL"])] == "Y")
UB = collections.defaultdict(list)
VAR = {}
f = io.open(os.path.expanduser(os.environ["MAPF"]), encoding="utf8")
h = f.readline().rstrip(NL).split("\t")
ix = dict((n, i) for i, n in enumerate(h))
for L in f:
    r = L.rstrip(NL).split("\t")
    if len(r) < 15 or r[ix["seg_type"]].startswith("PARTIAL"):
        continue
    k = (r[0], r[1])
    u = int(r[ix["array_pos"]])
    UB[k].append((u, int(r[ix["unit_start"]]) - 1, int(r[ix["unit_end"]]), int(r[ix["N"]])))
    VAR[(r[0], r[1], u)] = set() if r[ix["varlist"]] == "-" else set(r[ix["varlist"]].split(","))
ASGH = {}
SIDE = {}
for i, L in enumerate(io.open(os.path.expanduser("~/HAP_READ_FINAL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    r = L.rstrip(NL).split("\t")
    if len(r) >= 7 and r[4] != "-":
        ASGH[(r[0], r[1])] = r[4]
        SIDE[(r[0], r[1])] = r[2]

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r"([:*+-])([A-Za-z0-9]+)")
UB8 = collections.defaultdict(list)
N8 = {}
VAR8 = {}
f8 = io.open(os.path.expanduser(os.environ["G8MAP"]), encoding="utf8")
h8 = f8.readline().rstrip(NL).split(chr(9))
i8 = dict((n, i) for i, n in enumerate(h8))
for L8 in f8:
    r8 = L8.rstrip(NL).split(chr(9))
    if r8[i8["seg_type"]] != "UNIT":
        continue
    UB8[(r8[0], r8[1])].append((int(r8[i8["array_pos"]]), int(r8[i8["unit_start"]]) - 1, int(r8[i8["unit_end"]])))
    N8[(r8[0], r8[1])] = int(r8[i8["N"]])
    VAR8[(r8[0], r8[1], int(r8[i8["array_pos"]]))] = set() if r8[i8["varlist"]] == "-" else set(r8[i8["varlist"]].split(","))
C2SKIP = set(tuple(x.rsplit("_", 1)) for x in os.environ.get("C2SKIP", "").split(",") if x)
TABJ = lambda x: chr(9).join(x)


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


DONE = set()
_dp = os.path.join(W, "done.txt")
if os.path.exists(_dp):
    DONE = set(x.strip() for x in io.open(_dp, encoding="utf8") if x.strip())
OUT = io.open(os.path.join(W, "rjf_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
FLP = io.open(os.path.join(W, "flip_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
PRD = io.open(os.path.join(W, "readunit_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
OUT8 = io.open(os.path.join(W, "rjf8_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
R8 = random.Random(20260907)

for smp in [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]:
    if smp in DONE:
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
    if p is None:
        ch = sorted(x for x in
                    glob.glob(os.path.expanduser("~/giab5s_seq/out/%s_*.hits.fastq.gz" % smp))
                    + glob.glob(os.path.expanduser("~/hg008/out/%s_*.hits.fastq.gz" % smp))
                    if os.path.getsize(x) > 1000)
        if ch:
            p = os.path.join(W, "%s.cat.fastq.gz" % smp)
            if not (os.path.exists(p) and os.path.getsize(p) > 1000):
                with io.open(p, "wb") as g:
                    for x in ch:
                        g.write(io.open(x, "rb").read())
    slots = sorted(s for (a, s) in UB if a == smp)
    if not p or not slots:
        continue
    subq = os.path.join(W, "%s.lab.fastq" % smp)
    nlab = 0
    TG = {}
    with io.open(subq, "w", encoding="utf8", newline=NL) as g_, gzip.open(p, "rt", encoding="utf8") as fh_:
        rec_ = []
        for L_ in fh_:
            rec_.append(L_)
            if len(rec_) == 4:
                hd_ = rec_[0][1:].split()
                rid_ = hd_[0]
                if (smp, rid_) in ASGH or (smp, rid_[:12]) in ASGH:
                    g_.write("".join(rec_))
                    TG[rid_] = [x for x in hd_[1:] if x.startswith("MM:Z:") or x.startswith("ML:B:")]
                    nlab += 1
                rec_ = []
    p = subq
    ref = os.path.join(W, "%s.ref.fa" % smp)
    ARR = {}
    with io.open(ref, "w", newline=NL) as g:
        for sl in slots:
            q = os.path.expanduser(os.environ["CMPD"] + "/" + SEQF.get((smp, sl), "___"))
            if not os.path.exists(q):
                continue
            a = rd(q)
            if RCF.get((smp, sl)):
                a = rc(a)
            ARR[sl] = a
            g.write(">%s_h%s%s%s%s" % (smp, sl, NL, rc(FD[:FL]) + a + rc(FP[-FL:]), NL))
    pr = subprocess.Popen(["minimap2", "-a", "-x", "map-ont", "-t", TH, "-K", "100M",
                           "--secondary=no", ref, p],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    RXI = []
    SAMH = []
    REC = []            # (assigned slot, aligned contig slot, [(unit, sequence)...])
    AG = collections.Counter()
    for L in pr.stdout:
        if L.startswith("@"):
            SAMH.append(L)
            continue
        f = L.split("\t")
        if len(f) < 10 or f[2] == "*":
            continue
        if int(f[1]) & 0x900 or int(f[4]) < 20:
            continue
        asg = ASGH.get((smp, f[0])) or ASGH.get((smp, f[0][:12]))
        if asg is None:
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
            k, op = int(m.group(1)), m.group(2)
            if op == "S":
                q2 += k
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k):
                    m2[rpos + t] = q2 + t
                rpos += k
                q2 += k
            elif op == "I":
                q2 += k
            elif op in "DN":
                rpos += k
        if not m2:
            continue
        lo, hi = min(m2), max(m2)
        Ld, Lp = FL, FL + len(a)
        if not ((lo < Ld - 300 and hi > Ld + 300) or (lo < Lp - 300 and hi > Lp + 300)):
            continue
        sd = SIDE.get((smp, f[0])) or SIDE.get((smp, f[0][:12])) or "?"
        AG[(sd, asg, ctg)] += 1
        seg = []
        segc = []
        for u, s0, e0, N in UB[(smp, ctg)]:
            A, B = s0 + Ld, e0 + Ld - 1
            if A in m2 and B in m2:
                seg.append((u, seq[m2[A]:m2[B] + 1]))
                segc.append((u, m2[A], m2[B]))
        seg8 = []
        if (smp, ctg) not in C2SKIP:
            for u8, s8, e8 in UB8[(smp, ctg)]:
                A8, B8 = s8 + Ld, e8 + Ld - 1
                if A8 in m2 and B8 in m2:
                    seg8.append((u8, m2[A8], m2[B8], seq[m2[A8]:m2[B8] + 1]))
        if seg:
            REC.append((sd, asg, ctg, seg))
            RXI.append((f[0], int(f[1]), segc, seg8, L))
    pr.stdout.close()
    rcode = pr.wait()
    try:
        os.remove(ref)
    except OSError:
        pass
    if rcode != 0:
        sys.stderr.write("  MM_FAIL %s rc=%s%s" % (smp, rcode, NL))
        sys.stderr.flush()
        continue
    # orientation is decided separately for each side (one bit per individual cannot fix the HG00106 type)
    ORI = {}
    STS = {}
    for sd in ("prox", "dist", "?"):
        ns = sum(v for (x, a, c), v in AG.items() if x == sd)
        sm = sum(v for (x, a, c), v in AG.items() if x == sd and a == c)
        if not ns:
            continue
        fr2 = sm / float(ns)
        if fr2 >= HI:
            ORI[sd], STS[sd] = False, "OK"
        elif fr2 <= LO:
            ORI[sd], STS[sd] = True, "FLIP"
        else:
            ORI[sd], STS[sd] = (fr2 < 0.5), "AMBIG"
        FLP.write("%s\t%s\t%d\t%d\t%.4f\t%s%s" % (smp, sd, ns, sm, fr2, STS[sd], NL))
    n = sum(AG.values())
    same = sum(v for (x, a, c), v in AG.items() if a == c)
    fr = same / float(n) if n else -1.0
    # individual summary: if the two sides have different orientations, it signals that the pairing was wrong
    pr_, ds_ = ORI.get("prox"), ORI.get("dist")
    if pr_ is None or ds_ is None:
        st = "ONESIDE"
    elif pr_ != ds_:
        st = "PAIRERR"
    elif pr_:
        st = "FLIP"
    else:
        st = "OK"
    if n == 0:
        st = "NOREAD"
    FLP.write("%s\tALL\t%d\t%d\t%.4f\t%s%s" % (smp, n, same, fr, st, NL))
    FLP.flush()
    # with the correspondence decided, collect only the reads matching that orientation
    ANC = collections.defaultdict(list)
    ANC8 = collections.defaultdict(list)
    SAMO = io.open(os.path.join(W, "%s.placed.sam" % smp), "w", encoding="utf8", newline=NL)
    SAMO.write("".join(SAMH))
    for ri, (sd, asg, ctg, seg) in enumerate(REC):
        fl = ORI.get(sd)
        if fl is None:
            continue
        ok = (asg != ctg) if fl else (asg == ctg)
        if not ok:
            continue
        for u, sq in seg:
            ANC[(ctg, u)].append(sq)
        rid_, flg_, segc_, seg8_, sam_ = RXI[ri]
        strand_ = "1" if flg_ & 16 else "0"
        for u, qs_, qe_ in segc_:
            PRD.write(TABJ([smp, rid_, sd, asg, ctg, strand_, "G6", str(u), str(qs_), str(qe_)]) + NL)
        for u8, qs8, qe8, sq8 in seg8_:
            PRD.write(TABJ([smp, rid_, sd, asg, ctg, strand_, "G8", str(u8), str(qs8), str(qe8)]) + NL)
            ANC8[(ctg, u8)].append(sq8)
        tg_ = TG.get(rid_) or TG.get(rid_[:12]) or []
        SAMO.write(sam_.rstrip(NL) + (chr(9) + chr(9).join(tg_) if tg_ else "") + NL)
    SAMO.close()
    keys = [k for k, v in ANC.items() if len(v) >= 3]
    for (sl, u) in sorted(keys):
        lst = ANC[(sl, u)]
        V = varof(lst, "%s_%s_%s_%d" % (TAG, smp, sl, u))
        idx = list(range(len(lst)))
        random.shuffle(idx)
        hm = len(idx) // 2
        SA = set(idx[:hm])
        c = collections.Counter()
        ca = collections.Counter()
        cb = collections.Counter()
        for i in range(len(lst)):
            for x in V.get(i, set()):
                c[x] += 1
                if i in SA:
                    ca[x] += 1
                else:
                    cb[x] += 1
        M = VAR.get((smp, sl, u))
        N = UB[(smp, sl)][0][3]
        upos = N + 1 - u
        tok = ";".join("%s:%d:%d:%d" % (x, c[x], ca[x], cb[x]) for x in sorted(c))
        OUT.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s%s" %
                  (smp, sl, upos, u, len(lst), hm, len(idx) - hm, st,
                   tok or "-",
                   (",".join(sorted(M)) or "-") if M is not None else "NA", NL))
    OUT.flush()
    keys8 = [k for k, v in ANC8.items() if len(v) >= 3]
    for (sl, u) in sorted(keys8):
        lst = ANC8[(sl, u)]
        V = varof(lst, "%s_%s_%s_%d_g8" % (TAG, smp, sl, u))
        idx = list(range(len(lst)))
        R8.shuffle(idx)
        hm = len(idx) // 2
        SA = set(idx[:hm])
        c = collections.Counter()
        ca = collections.Counter()
        cb = collections.Counter()
        for i in range(len(lst)):
            for x in V.get(i, set()):
                c[x] += 1
                if i in SA:
                    ca[x] += 1
                else:
                    cb[x] += 1
        M = VAR8.get((smp, sl, u))
        N = N8[(smp, sl)]
        tok = ";".join("%s:%d:%d:%d" % (x, c[x], ca[x], cb[x]) for x in sorted(c))
        OUT8.write(TABJ([smp, sl, str(N + 1 - u), str(u), str(len(lst)), str(hm), str(len(idx) - hm), st, tok or "-",
                         (",".join(sorted(M)) or "-") if M is not None else "NA"]) + NL)
    OUT8.flush()
    PRD.flush()
    try:
        os.remove(subq)
    except OSError:
        pass
    io.open(_dp, "a", encoding="utf8", newline=NL).write(smp + NL)
    sys.stderr.write("  %-10s 정렬 %4d · 라벨 %-5s(%.2f) · 자리 %d%s"
                     % (smp, n, st, fr, len(keys), NL))
    sys.stderr.flush()
OUT.close()
FLP.close()
print("done")
