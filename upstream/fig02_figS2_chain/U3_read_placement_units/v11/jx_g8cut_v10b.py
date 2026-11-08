# -*- coding: utf-8 -*-
# cut junction-placed reads on the G8v3 frame (2026-09-11), run on the node, nothing is re-aligned
#   material: alignment records of accepted reads stored by the reproduction run (<sample>.placed.sam, original read sequence, quality and MM/ML unchanged, passed comparison with the original)
#   reference: rebuilt exactly as used for that alignment = rc(distal flank 5000) + canonical + rc(proximal flank 5000)
#          main edition: hap_seq_node, cmp4d, orientation CONSENSUS_SUM_G6 revcomp (exactly the original rjfreq2 input)
#          c2 edition  : hap_seq_node_v7, cmp4d_v7, orientation G8v2 SUM revcomp (consensus2 canonical 6 slots, latest canonical)
#   boundaries: G8v3 unit coordinates (canonical) + 5000 transferred onto the read by alignment. For a gene-end T followed by the terminator TTTT run,
#          the two non-T bases outside the T run are used as reference points, and the same number of Ts as in the canonical is given to the end of the preceding unit (same rule as ao_g8cut)
#   numbering: for junction reads the alignment coordinate is the position (the original rjfreq2 also used the aligned unit numbers as they were) -> G8v3 array_pos
#   variants: same method as the G8v3 map (asm20, --cs=long, end fill, restoration check)
#   the main edition skips the 6 consensus2 slots (the G8 coordinates of those slots are based on the consensus2 canonical and cannot be used with the old reference alignment)
#   Arguments: main|c2, output folder, number of parallel jobs
import io
import os
import re
import sys
import glob
import subprocess
import collections
import multiprocessing
import socket

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
PASS, OUTD = sys.argv[1], H(sys.argv[2])
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 6
FL = 5000
CMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CS = re.compile(r"([:*+\-=])([A-Za-z0-9]+)")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
REF = H("~/s5s_tuc/ref")
CANON = os.path.join(REF, "unit_canon.fa")
C2 = set([("HG02683", "1"), ("NA19090", "1"), ("HG03126", "1"), ("HG03366", "1"), ("NA19120", "1"), ("HG03198", "1")])


def rc(s):
    return s.translate(CMP)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


CANSEQ = rd(CANON).upper()
FD = rd(os.path.join(REF, "FLANK_distal_50000.fa"))
FP = rd(os.path.join(REF, "FLANK_proximal_50000.fa"))


def call_variants(seqs, tag, wd):
    q = os.path.join(wd, "q_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "2", "--secondary=no", CANON, q],
                       capture_output=True, text=True)
    OUT = {}
    for L in r.stdout.splitlines():
        f2 = L.split(TAB)
        if len(f2) < 12:
            continue
        i = int(f2[0][1:])
        if i in OUT:
            continue
        cs = ""
        for x in f2[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        qs, qe, ts, te = int(f2[2]), int(f2[3]), int(f2[7]), int(f2[8])
        toks, reb, tp = [], [], ts
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == "=":
                reb.append(v.upper())
                tp += len(v)
            elif op == ":":
                n = int(v)
                reb.append(CANSEQ[tp:tp + n])
                tp += n
            elif op == "*":
                toks.append("%d%s" % (tp + 1, v[1].upper()))
                reb.append(v[1].upper())
                tp += 1
            elif op == "-":
                toks.append("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                toks.append("%dI%d" % (tp + 1, len(v)))
                reb.append(v.upper())
        qseq = seqs[i]
        kt = 0
        while te + kt < 2240 and qe + kt < len(qseq):
            if qseq[qe + kt] != CANSEQ[te + kt]:
                toks.append("%d%s" % (te + kt + 1, qseq[qe + kt]))
            kt += 1
        if te + kt < 2240:
            toks.append("%dD%d" % (te + kt + 1, 2240 - (te + kt)))
        kh = 0
        while ts - 1 - kh >= 0 and qs - 1 - kh >= 0:
            if qseq[qs - 1 - kh] != CANSEQ[ts - 1 - kh]:
                toks.append("%d%s" % (ts - kh, qseq[qs - 1 - kh]))
            kh += 1
        OUT[i] = dict(toks=sorted(set(toks)), qs=qs, qe=qe, reb="".join(reb), fill=kh + kt, cov_e=te + kt, cov_s=ts - kh)
    os.remove(q)
    return OUT


def verify(seq, res):
    if res is None:
        return "NOALN"
    if res["reb"] != seq[res["qs"]:res["qe"]]:
        return "SEQ_MISMATCH"
    if res["cov_e"] < 2240 or res["cov_s"] > 0:
        return "NO_FULL"
    return "OK_FILLED" if res["fill"] else "OK"


def boundary(A, m2, b, side):
    """Boundary after gene end b (1-based) of reference A (0-based) -> read coordinate (0-based on SAM SEQ). Same rule as ao_g8cut.py"""
    i = b - 1
    if 0 <= i < len(A) and A[i] == "T":
        rs = i
        while rs - 1 >= 0 and A[rs - 1] == "T":
            rs -= 1
        re_ = i
        while re_ + 1 < len(A) and A[re_ + 1] == "T":
            re_ += 1
        a0, a1 = rs - 1, re_ + 1
        k1 = i - rs + 1
        if a0 in m2 and a1 in m2:
            q0, q1 = m2[a0], m2[a1]
            lr = q1 - q0 - 1
            qe = q0 + min(k1, max(lr, 0))
            return (qe if side == "end" else qe + 1), "RUN"
    j = i if side == "end" else i + 1
    if j in m2:
        return m2[j], "DIRECT"
    for d in range(1, 50):
        k = j - d if side == "end" else j + d
        if k in m2:
            return m2[k], "NEAREST"
    return None, "FAIL"


def init():
    global SEQF, RCF, CMPD, U8, SRC
    if PASS == "main":
        hsn, sumf, col, CMPD = H("~/hap_seq_node.tsv"), H("~/CONSENSUS_SUM_G6.tsv"), 5, H("~/cmp4d")
    else:
        hsn, sumf, col, CMPD = H("~/g8in/hap_seq_node_v7.tsv"), H("~/g8in/CONSENSUS_SUM_G8v2_full.tsv"), 4, H("~/cmp4d_v7")
    SEQF = {}
    for L in io.open(hsn, encoding="utf8"):
        f = L.rstrip(NL).split(TAB)
        if len(f) >= 3 and f[0] != "sample":
            SEQF[(f[0], f[1])] = f[2]
    RCF = {}
    for i, L in enumerate(io.open(sumf, encoding="utf8")):
        if i:
            f = L.rstrip(NL).split(TAB)
            RCF[(f[0], f[1])] = (f[col] == "Y")
    U8 = collections.defaultdict(list)
    f = io.open(H("~/g8jxmap_v10/CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
    ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    for L in f:
        r = L.rstrip(NL).split(TAB)
        if r[ix["seg_type"]] == "UNIT":
            U8[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]), int(r[ix["array_pos"]])))
    SRC = {}
    for L in io.open(H("~/srcmap_%s.tsv" % socket.gethostname()), encoding="utf8"):
        r = L.rstrip(NL).split(TAB)
        SRC[(r[1], r[2])] = r[3]


def refseq(smp, sl):
    a = rd(os.path.join(CMPD, SEQF[(smp, sl)]))
    if RCF.get((smp, sl)):
        a = rc(a)
    return (rc(FD[:FL]) + a + rc(FP[-FL:])).upper()


def work(sam):
    smp = os.path.basename(sam).split(".")[0]
    wd = os.path.join(OUTD, "tmp_" + smp)
    os.makedirs(wd, exist_ok=True)
    job = "JX" if PASS == "main" else "JX_C2"
    lab = SRC.get((job, smp), "UNKNOWN")
    REFS = {}
    seqs, meta, stat = [], [], collections.Counter()
    for L in io.open(sam, encoding="utf8"):
        if L.startswith("@"):
            continue
        f = L.rstrip(NL).split(TAB)
        smp_, sl = f[2].rsplit("_h", 1)
        if (PASS == "main" and (smp_, sl) in C2) or (PASS == "c2" and (smp_, sl) not in C2):
            stat["건너뜀(참조 틀 다름)"] += 1
            continue
        if (smp_, sl) not in REFS:
            REFS[(smp_, sl)] = refseq(smp_, sl)
        A = REFS[(smp_, sl)]
        rpos = int(f[3]) - 1
        q2 = 0
        m2 = {}
        for m in CIG.finditer(f[5]):
            k, op = int(m.group(1)), m.group(2)
            if op == "S":
                q2 += k
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
        seq = f[9]
        L0 = len(seq)
        rev = int(f[1]) & 16
        stat["read"] += 1
        for s8, e8, p8 in U8[(smp_, sl)]:
            sb, eb = s8 - 1 + FL, e8 + FL            # reference 1-based: gene end of the preceding unit, gene end of this unit
            if not (lo < sb - 1 and eb < hi):
                continue
            qs_, bs = boundary(A, m2, sb, "start")
            qe_, be = boundary(A, m2, eb, "end")
            if qs_ is None or qe_ is None or qe_ < qs_:
                stat["경계 실패"] += 1
                continue
            stat["경계 " + bs] += 1
            stat["경계 " + be] += 1
            u = seq[qs_:qe_ + 1]
            f0, f1 = (L0 - (qe_ + 1), L0 - qs_) if rev else (qs_, qe_ + 1)
            seqs.append(u)
            meta.append([smp_, f[0], sl, job, lab, "-" if rev else "+", str(p8), str(f0), str(f1), bs, be])
    V = call_variants(seqs, smp, wd) if seqs else {}
    rows = []
    for i, mt in enumerate(meta):
        res = V.get(i)
        vd = verify(seqs[i], res)
        toks = res["toks"] if (res and vd in ("OK", "OK_FILLED")) else None
        rows.append(mt + [vd, (",".join(toks) or "-") if toks is not None else "-"])
        stat["검산 " + vd] += 1
    for fn in glob.glob(os.path.join(wd, "*")):
        os.remove(fn)
    os.rmdir(wd)
    return smp, rows, stat


if __name__ == "__main__":
    init()
    os.makedirs(OUTD, exist_ok=True)
    src = H("~/rjrepro") if PASS == "main" else H("~/rjrepro_c2")
    sams = sorted(glob.glob(os.path.join(src, "*.placed.sam")))
    O = io.open(os.path.join(OUTD, "JX_G8_UNITS_%s_%s.tsv" % (PASS, socket.gethostname())), "w", encoding="utf8", newline=NL)
    TOT = collections.Counter()
    with multiprocessing.Pool(NP, initializer=init) as pool:
        for smp, rows, stat in pool.imap_unordered(work, sams):
            for r in rows:
                O.write(TAB.join(r) + NL)
            TOT.update(stat)
            O.flush()
    O.close()
    print("%s %s · 개체 %d · %s" % (socket.gethostname(), PASS, len(sams), " · ".join("%s %d" % kv for kv in sorted(TOT.items()))))
