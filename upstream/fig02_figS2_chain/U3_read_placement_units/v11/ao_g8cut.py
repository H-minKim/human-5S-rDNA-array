# -*- coding: utf-8 -*-
# Cut ARRAY_ONLY placed reads in the G8 frame (2026-09-11), run on local WSL
#   No re-placement: read, slot and unit numbers are from the saved v7 placement (AOJV)
#   ① raw reads: fastq pulled from the nodes (header MM/ML, sequence and quality as is)
#   ② alignment: to the canonical sequence of its own slot, in the saved unit range (shrunk into the canonical) +-250bp window (minimap2 -a -x map-ont --secondary=no -y)
#            (fixed) the old +-3 unit window let reads slip by 1~3 units in arrays of consecutive similar units (ARRAY_ONLY 11%)
#            saved placement check: per read, compare spacer variants with saved number g+d, and if clearly off (same number <0.6, neighbouring number >=0.8)
#            realign to the whole-read footprint window (+-300bp) moved to the saved position, and accept if it matches the saved number (recorded in READS column 11)
#            with -y the MM/ML of the fastq header follow into the SAM, and since the input is FASTQ, QUAL is kept as well
#   ③ boundaries: transfer the canonical (G8v3) unit boundaries to the read via the alignment. After the gene end (T) comes the terminator TTTT, so
#            the boundary lies in the middle of a T run -> use the two non-T bases outside the T run as anchors,
#            and in the T run of the read give the same number of T as the canonical (usually 1) to the end of the previous unit, the rest to the head of the next unit.
#            Boundaries not in a T run are transferred as is, and if no corresponding base exists, the nearest corresponding base is used (flag NEAREST)
#   ④ variants: same way as the G8v3 map (asm20, --cs=long, end sequence filling, restore check)
#   ⑤ numbering: canonical coordinate from the alignment -> G6 number (checked against the saved value) -> G8 number. Spacer variants confirm agreement with the saved number
#   args: work folder (aog8), individual list, parallel count
import io
import os
import re
import sys
import gzip
import glob
import bisect
import subprocess
import collections
import multiprocessing

NL = chr(10)
TAB = chr(9)
D = sys.argv[1]
LST = sys.argv[2]
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 6
CMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CS = re.compile(r"([:*+\-=])([A-Za-z0-9]+)")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
NUM = re.compile(r"\d+")
WIN = 3            # (old value, not used) 2026-09-11 with a +-3 unit window reads slipped by 1~3 units (ARRAY_ONLY 8%) -> replaced by PADBP
PADBP = 250        # window = saved G6 unit range (min~max) in canonical coordinates +- 250bp: removes the room to slip by one unit (2,240bp)
SHIFTMAX = 3       # range (+-units) for checking per-read number differences against the saved placement


def rc(s):
    return s.translate(CMP)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


REF = os.path.join(D, "s5s_tuc", "ref")
CANON = os.path.join(REF, "unit_canon.fa")
CANSEQ = rd(CANON).upper()
assert len(CANSEQ) == 2240


def samples(p):
    return set(L.split(TAB)[0] for L in io.open(p, encoding="utf8") if L.strip())


def loadmap(p, types=("UNIT",)):
    U = collections.defaultdict(list)
    f = io.open(p, encoding="utf8")
    ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    for L in f:
        r = L.rstrip(NL).split(TAB)
        if r[ix["seg_type"]] in types:
            U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]), int(r[ix["array_pos"]])))
    for k in U:
        U[k].sort()
    return U


def call_variants(seqs, tag, wd):
    q = os.path.join(wd, "q_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "1", "--secondary=no", CANON, q],
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
    """Transfer the boundary after gene end b (1-based) of canonical A (0-based string) to the read.
    side='end' gives the read end coordinate (inclusive) of the previous unit, 'start' the read start coordinate of the next unit. (read coordinates are 0-based, on SAM SEQ)"""
    i = b - 1                                   # 0-based gene end
    if 0 <= i < len(A) and A[i] == "T":
        rs = i
        while rs - 1 >= 0 and A[rs - 1] == "T":
            rs -= 1
        re_ = i
        while re_ + 1 < len(A) and A[re_ + 1] == "T":
            re_ += 1
        a0, a1 = rs - 1, re_ + 1                # anchor outside the T run (canonical 0-based)
        k1 = i - rs + 1                         # number of T belonging to the previous unit in the canonical
        if a0 in m2 and a1 in m2:
            q0, q1 = m2[a0], m2[a1]
            lr = q1 - q0 - 1                    # length of the T run of the read
            qe = q0 + min(k1, max(lr, 0))       # end of the previous unit (inclusive)
            return (qe if side == "end" else qe + 1), "RUN"
    # not a T run or no anchor found: transfer as is, otherwise the nearest corresponding base
    j = i if side == "end" else i + 1
    if j in m2:
        return m2[j], "DIRECT"
    for d in range(1, 50):
        k = j - d if side == "end" else j + d
        if k in m2:
            return m2[k] + (0 if side == "end" else 0), "NEAREST"
    return None, "FAIL"


def align_cut(smp, slot, rid, r16, cfg, lab, rec, A, u6, ws, we, wd):
    """Align the read to window [ws, we] (canonical 1-based) and cut it into G8 units.
    Returns: (one SAM line (canonical coordinates), [(unit sequence, meta)], (aligned canonical start0, aligned canonical end0 (exclusive), left clip, right clip)) or None"""
    wfa = os.path.join(wd, "w.fa")
    rfq = os.path.join(wd, "r.fq")
    io.open(wfa, "w", newline=NL).write(">%s_h%s%s%s%s" % (smp, slot, NL, A[ws - 1:we], NL))
    io.open(rfq, "w", newline=NL).write(NL.join(rec) + NL)
    r = subprocess.run(["minimap2", "-a", "-x", "map-ont", "--secondary=no", "-y", "-t", "1", wfa, rfq],
                       capture_output=True, text=True)
    prim = None
    for L in r.stdout.splitlines():
        if L.startswith("@"):
            continue
        f = L.split(TAB)
        if f[2] != "*" and not (int(f[1]) & 0x900):
            prim = f
            break
    if prim is None:
        return None
    # SAM moved to canonical coordinates (POS += window start - 1). The raw read sequence, quality and MM/ML are not touched
    prim2 = list(prim)
    prim2[3] = str(int(prim[3]) + ws - 1)
    rpos = int(prim[3]) - 1 + ws - 1      # canonical 0-based
    ref0 = rpos
    q2 = 0
    m2 = {}
    cig = [(int(m.group(1)), m.group(2)) for m in CIG.finditer(prim[5])]
    for k, op in cig:
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
    lclip = cig[0][0] if cig and cig[0][1] in "SH" else 0
    rclip = cig[-1][0] if cig and cig[-1][1] in "SH" else 0
    seq = prim[9]
    flag = int(prim[1])
    L0 = len(seq)
    lo_r, hi_r = min(m2), max(m2)
    units = []
    for s8, e8, p8 in U8[(smp, slot)]:
        # unit p8 = [s8, e8] : start boundary = after the gene end of the previous unit (s8-1), end boundary = after the gene end of this unit (e8)
        if not (lo_r < s8 - 2 and e8 < hi_r):
            continue
        qs_, bs = boundary(A, m2, s8 - 1, "start")
        qe_, be = boundary(A, m2, e8, "end")
        if qs_ is None or qe_ is None or qe_ < qs_:
            continue
        useq = seq[qs_:qe_ + 1]
        # coordinates on the raw read (fastq)
        if flag & 16:
            f0, f1 = L0 - (qe_ + 1), L0 - qs_
        else:
            f0, f1 = qs_, qe_ + 1
        mid = s8 + 1000                     # center of the spacer of this unit (canonical 1-based)
        i6 = bisect.bisect_right([x[0] for x in u6], mid) - 1
        g6 = u6[i6][2] if (i6 >= 0 and u6[i6][0] <= mid <= u6[i6][1]) else None
        units.append((useq, [smp, rid, r16, slot, cfg, lab, "-" if not (flag & 16) else "+rc", p8, g6, f0, f1, bs, be]))
    return TAB.join(prim2), units, (ref0, rpos, lclip, rclip)


def spacer(toks):
    return set(t for t in toks if int(NUM.match(t).group()) <= 2121)


def jac(a, b):
    un = a | b
    return len(a & b) / float(len(un)) if un else 1.0


def read_offset(items, P):
    """items = [(g6, spacer variant set)] (only units with called variants) -> {d: mean Jaccard with saved number g6+d}"""
    S = {}
    for d in range(-SHIFTMAX, SHIFTMAX + 1):
        v = [jac(t, spacer(P[g + d])) for g, t in items if g is not None and g + d in P]
        if len(v) >= max(2, len(items) // 2):
            S[d] = sum(v) / len(v)
    return S


def shift_pos(u6, pos1, d):
    """Canonical coordinate pos1 (1-based) moved by d G6 units (the relative position within the unit is kept)"""
    st = [x[0] for x in u6]
    i = bisect.bisect_right(st, pos1) - 1
    i = max(i, 0)
    j = i + d
    if j < 0 or j >= len(u6):
        return None
    return pos1 + (u6[j][0] - u6[i][0])


def relabel(res, bd):
    """Reads that do not attach to the saved position even after realignment: the cut read coordinates and variants are kept, only the unit numbers are shifted by the saved difference bd.
    The SAM record keeps the actual alignment and is flagged with a ZR tag (the POS of that read is off by bd units from the saved numbering)"""
    sam, units, foot = res
    units2 = [(us, mt[:7] + [mt[7] + bd, (mt[8] + bd) if mt[8] is not None else None] + mt[9:]) for us, mt in units]
    return (sam + TAB + "ZR:Z:relabeled_d=%+d" % bd, units2, foot)


def work(smp):
    wd = os.path.join(D, "tmp", smp)
    os.makedirs(wd, exist_ok=True)
    cfg = CFG[smp]
    sv = SAVED[(cfg, smp)]
    files = sorted(glob.glob(os.path.join(D, "aoreads", "%s.%s.*.fastq.gz" % (smp, cfg))))
    out_u, out_r, out_s = [], [], []
    FQ = {}
    for fp in files:
        lab = fp.split(".")[-3]
        with gzip.open(fp, "rt", encoding="utf8") as fh:
            rec = []
            for x in fh:
                rec.append(x.rstrip(NL))
                if len(rec) == 4:
                    rid = rec[0][1:].split()[0]
                    FQ.setdefault(rid[:16], (rid, rec, lab))
                    rec = []
    RES = {}                  # r16 -> (sam, units, foot)
    NOTE = {}                 # r16 -> realignment record
    # ── pass 1: window = saved G6 unit range (shrunk into the canonical) +- PADBP
    for r16, (slot, P) in sorted(sv.items()):
        if r16 not in FQ:
            out_r.append([smp, r16, "-", slot, cfg, "-", "NOT_IN_FASTQ"])
            continue
        rid, rec, lab = FQ[r16]
        A = ARR[(smp, slot)]
        u6 = U6[(smp, slot)]
        p6 = dict((p, (s, e)) for s, e, p in u6)
        inr = sorted(p for p in p6 if min(P) <= p <= max(P))      # units of the saved range present in the canonical (saved numbers beyond the array end are clipped)
        if not inr:
            out_r.append([smp, r16, rid, slot, cfg, lab, "NO_WINDOW"])
            continue
        if max(P) > inr[-1] or min(P) < inr[0]:
            NOTE[r16] = "saved_beyond_canon(%d-%d>%d-%d)" % (min(P), max(P), inr[0], inr[-1])
        ws, we = max(1, p6[inr[0]][0] - PADBP), min(len(A), p6[inr[-1]][1] + PADBP)
        x = align_cut(smp, slot, rid, r16, cfg, lab, rec, A, u6, ws, we, wd)
        if x is None:
            out_r.append([smp, r16, rid, slot, cfg, lab, "NOALN"])
            continue
        RES[r16] = x
    order = [(r16, i) for r16 in sorted(RES) for i in range(len(RES[r16][1]))]
    V = call_variants([RES[r][1][i][0] for r, i in order], smp, wd) if order else {}
    CALL = {}                 # (r16, i) -> (vd, toks)
    for k, (r16, i) in enumerate(order):
        res = V.get(k)
        vd = verify(RES[r16][1][i][0], res)
        CALL[(r16, i)] = (vd, set(res["toks"]) if (res and vd in ("OK", "OK_FILLED")) else None)
    # ── pass 2: if the saved placement and the per-read numbers are clearly off (same number < 0.6, neighbouring number >= 0.8)
    #          realign to the 'whole-read footprint' window moved to the saved position (+-300bp at both ends). Slipping by one unit leaves 2kb of the read unaligned, which is penalized
    REDO = {}
    RELAB = {}
    for r16, (sam, units, foot) in RES.items():
        slot, P = sv[r16]
        items = [(units[i][1][8], spacer(CALL[(r16, i)][1])) for i in range(len(units)) if CALL[(r16, i)][1] is not None]
        S = read_offset(items, P)
        if 0 not in S:
            continue
        bd = max(S, key=lambda d: (S[d], -abs(d)))
        if bd == 0 or not (S[bd] >= 0.8 and S[0] < 0.6):
            continue
        u6 = U6[(smp, slot)]
        A = ARR[(smp, slot)]
        ref0, ref1, lcl, rcl = foot
        fs, fe = shift_pos(u6, ref0 + 1, bd), shift_pos(u6, ref1, bd)
        if fs is None or fe is None:
            NOTE[r16] = (NOTE.get(r16, "") + ";" if r16 in NOTE else "") + "shift_d=%+d_no_room" % bd
            RELAB[r16] = bd
            continue
        ws2, we2 = max(1, fs - lcl - 300), min(len(A), fe + rcl + 300)
        rid, rec, lab = FQ[r16]
        y = align_cut(smp, slot, rid, r16, cfg, lab, rec, A, u6, ws2, we2, wd)
        if y is None or not y[1]:
            NOTE[r16] = (NOTE.get(r16, "") + ";" if r16 in NOTE else "") + "shift_d=%+d_realign_noaln" % bd
            RELAB[r16] = bd
            continue
        REDO[r16] = (y, bd, S[0], S[bd])
    if REDO:
        order2 = [(r16, i) for r16 in sorted(REDO) for i in range(len(REDO[r16][0][1]))]
        V2 = call_variants([REDO[r][0][1][i][0] for r, i in order2], smp + "_re", wd)
        CALL2 = {}
        for k, (r16, i) in enumerate(order2):
            res = V2.get(k)
            vd = verify(REDO[r16][0][1][i][0], res)
            CALL2[(r16, i)] = (vd, set(res["toks"]) if (res and vd in ("OK", "OK_FILLED")) else None)
        for r16, (y, bd, s0, sb) in REDO.items():
            slot, P = sv[r16]
            units = y[1]
            items = [(units[i][1][8], spacer(CALL2[(r16, i)][1])) for i in range(len(units)) if CALL2[(r16, i)][1] is not None]
            S2 = read_offset(items, P)
            if 0 in S2 and S2[0] >= 0.7 and S2[0] > s0 + 0.1 and max(S2, key=lambda d: (S2[d], -abs(d))) == 0:
                RES[r16] = y
                for i in range(len(units)):
                    CALL[(r16, i)] = CALL2[(r16, i)]
                NOTE[r16] = (NOTE.get(r16, "") + ";" if r16 in NOTE else "") + "realigned_d=%+d(%.3f->%.3f)" % (bd, s0, S2[0])
            else:
                NOTE[r16] = (NOTE.get(r16, "") + ";" if r16 in NOTE else "") + "shift_d=%+d_realign_failed(%.3f,%.3f->%s)" % (bd, s0, sb, "%.3f" % S2[0] if 0 in S2 else "-")
                RELAB[r16] = bd
    for r16, bd in RELAB.items():
        RES[r16] = relabel(RES[r16], bd)
        NOTE[r16] = NOTE[r16] + ";relabeled_d=%+d" % bd
    # ── output
    per = collections.defaultdict(list)
    for r16 in sorted(RES):
        sam, units, foot = RES[r16]
        out_s.append(sam)
        slot, P = sv[r16]
        for i, (useq, mt) in enumerate(units):
            vd, toks = CALL[(r16, i)]
            sc = "-"
            if toks is not None and mt[8] in P:
                sc = "%.3f" % jac(spacer(toks), spacer(P[mt[8]]))
                per[r16].append(float(sc))
            out_u.append(mt + [vd, (",".join(sorted(toks)) or "-") if toks is not None else "-", sc,
                               "Y" if mt[8] in P else "N"])
    for r16, (slot, P) in sorted(sv.items()):
        if r16 in FQ and not any(x[1] == r16 for x in out_r):
            s = per.get(r16, [])
            out_r.append([smp, r16, FQ[r16][0], slot, cfg, FQ[r16][2], "OK" if s else "NOUNIT",
                          str(len(s)), "%.3f" % (sum(s) / len(s)) if s else "-", str(len(P)), NOTE.get(r16, "-")])
    for fn in glob.glob(os.path.join(wd, "*")):
        os.remove(fn)
    os.rmdir(wd)
    return smp, out_u, out_r, out_s


def init():
    global CFG, SAVED, ARR, U6, U8
    AGS = samples(os.path.join(D, "AOJV_GIAB.tsv"))
    ACS = samples(os.path.join(D, "AOJV8.tsv"))
    AFS = samples(os.path.join(D, "AOJV5.tsv"))
    SAVED = collections.defaultdict(dict)
    for cfg, fn in (("AOJV3", "AOJV3.tsv"), ("AOJV5", "AOJV5.tsv"), ("AOJV8", "AOJV8.tsv"), ("AOJVG", "AOJV_GIAB.tsv")):
        for L in io.open(os.path.join(D, fn), encoding="utf8"):
            r = L.rstrip(NL).split(TAB)
            if len(r) < 5:
                continue
            e = SAVED[(cfg, r[0])].setdefault(r[1], [r[2], {}])
            e[1][int(r[3])] = set() if r[4] == "-" else set(r[4].split(","))
    CFG = {}
    for s in set(k[1] for k in SAVED):
        CFG[s] = "AOJVG" if s in AGS else "AOJV8" if s in ACS else "AOJV5" if s in AFS else "AOJV3"
    SEQF = {}
    for L in io.open(os.path.join(D, "g8in", "hap_seq_node_v7.tsv"), encoding="utf8"):
        f = L.rstrip(NL).split(TAB)
        if len(f) >= 3:
            SEQF[(f[0], f[1])] = f[2]
    FLIP = {}
    f = io.open(os.path.join(D, "SUM_G8.tsv"), encoding="utf8")
    f.readline()
    for L in f:
        r = L.rstrip(NL).split(TAB)
        FLIP[(r[0], r[1])] = (r[4] == "Y")
    ARR = {}
    for k, fn in SEQF.items():
        a = rd(os.path.join(D, "cmp4d_v7", fn)).upper()
        ARR[k] = rc(a) if FLIP.get(k) else a
    U6 = loadmap(os.path.join(D, "MAP_G6C2.tsv"), ("UNIT", "UNIT_LAST"))   # the last complete unit of the array (UNIT_LAST) is also used for the window and G6 numbers (without it, reads whose saved max number is the last unit become NO_WINDOW)
    U8 = loadmap(os.path.join(D, "MAP_G8.tsv"))


if __name__ == "__main__":
    init()
    smps = [x.strip() for x in io.open(LST, encoding="utf8") if x.strip()]
    smps = [s for s in smps if s in CFG]
    OU = io.open(os.path.join(D, "AO_G8_UNITS.tsv"), "a", encoding="utf8", newline=NL)
    OR = io.open(os.path.join(D, "AO_G8_READS.tsv"), "a", encoding="utf8", newline=NL)
    os.makedirs(os.path.join(D, "sam"), exist_ok=True)
    with multiprocessing.Pool(NP, initializer=init) as pool:
        for smp, uu, rr, ss in pool.imap_unordered(work, smps):
            for x in uu:
                OU.write(TAB.join(str(y) for y in x) + NL)
            for x in rr:
                OR.write(TAB.join(str(y) for y in x) + NL)
            io.open(os.path.join(D, "sam", "%s.%s.sam" % (smp, CFG[smp])), "w", encoding="utf8", newline=NL).write(NL.join(ss) + (NL if ss else ""))
            OU.flush()
            OR.flush()
            sys.stderr.write("  %s 유닛 %d · read %d%s" % (smp, len(uu), len(rr), NL))
    OU.close()
    OR.close()
    print("done")
