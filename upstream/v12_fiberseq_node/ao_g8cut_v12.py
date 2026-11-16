# -*- coding: utf-8 -*-
# Cut the ARRAY_ONLY placed reads into the G8 frame, run locally under WSL.
#   The placement is not redone: the read, the slot and the unit index are taken from the stored v7 placement (AOJV) unchanged.
# ① source reads : fastq taken from the compute node (header MM/ML, sequence and quality unchanged)
#   2) Alignment: against its own slot's canonical sequence, over the stored unit interval (clipped into the canonical sequence) plus or minus 250 bp (minimap2 -a -x map-ont --secondary=no -y).
#            A window of plus or minus 3 units lets a read slip by 1 to 3 units in an array of similar consecutive units, affecting 11% of ARRAY_ONLY reads.
#            A check against the stored placement: the spacer variants of each read are compared with the stored index g+d, and where they disagree clearly (the same index below 0.6 and a neighbouring index at 0.8 or above)
#            the read is realigned against the whole-footprint window at the stored position, plus or minus 300 bp, and accepted if it then matches the stored index; this is recorded in column 11 of READS.
#            -y carries the MM/ML from the fastq header through into the SAM, and feeding FASTQ in preserves the QUAL as well.
#   3) Boundaries: the canonical (G8v3) unit boundaries are carried onto the read through the alignment. The gene ends in T and is followed by the TTTT termination signal,
#            so the boundary falls in the middle of a run of Ts. The two non-T bases outside that run are taken as the anchors,
#            and of the read's run of Ts the same number as in the canonical sequence, usually one, is given to the end of the preceding unit and the rest to the start of the next.
#            A boundary not in a run of Ts is carried straight across, and where no corresponding base exists the nearest one is used, marked NEAREST.
#   4) Variants: the same way as the G8v3 map (asm20, --cs=long, the terminal sequence filled in, and a reconstruction check).
#   5) Indexing: the canonical coordinate from the alignment -> the G6 index, checked against the stored value -> the G8 index. The spacer variants confirm the agreement with the stored index.
#   Arguments: the working folder (aog8), the list of individuals, and the number of parallel jobs.
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
WIN = 3            # not used; at a window of plus or minus 3 units a read slips by 1 to 3 units, affecting 8% of ARRAY_ONLY reads, hence PADBP
PADBP = 250        # the window is the canonical-coordinate span of the stored G6 units, minimum to maximum, plus or minus 250 bp, which removes the room to slip by a whole unit of 2,240 bp
SHIFTMAX = 3       # the range, in units, over which the per-read index is compared with the stored placement


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
    """Carry the boundary after gene end b (1-based) of canonical A (0-based string) over to the read.
    side='end' gives the read end coordinate (inclusive) of the preceding unit, 'start' the read start coordinate of the next unit. (read coordinates 0-based, in SAM SEQ)"""
    i = b - 1                                   # the 0-based gene end
    if 0 <= i < len(A) and A[i] == "T":
        rs = i
        while rs - 1 >= 0 and A[rs - 1] == "T":
            rs -= 1
        re_ = i
        while re_ + 1 < len(A) and A[re_ + 1] == "T":
            re_ += 1
        a0, a1 = rs - 1, re_ + 1                # the anchors outside the run of Ts, 0-based in the canonical sequence
        k1 = i - rs + 1                         # the number of Ts belonging to the preceding unit in the canonical sequence
        if a0 in m2 and a1 in m2:
            q0, q1 = m2[a0], m2[a1]
            lr = q1 - q0 - 1                    # the length of the run of Ts in the read
            qe = q0 + min(k1, max(lr, 0))       # the end of the preceding unit, inclusive
            return (qe if side == "end" else qe + 1), "RUN"
    # where it is not a run of Ts, or the anchors cannot be found: carried straight across, and failing that the nearest corresponding base
    j = i if side == "end" else i + 1
    if j in m2:
        return m2[j], "DIRECT"
    for d in range(1, 50):
        k = j - d if side == "end" else j + d
        if k in m2:
            return m2[k] + (0 if side == "end" else 0), "NEAREST"
    return None, "FAIL"


def align_cut(smp, slot, rid, r16, cfg, lab, rec, A, u6, ws, we, wd):
    """Align the read to window [ws, we] (canonical, 1-based) and cut it into G8 units.
    Returns : (one SAM line (canonical coordinates), [(unit sequence, meta)], (aligned canonical start 0-based, aligned canonical end 0-based (exclusive), left clip, right clip)) or None"""
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
    # the SAM carried into canonical coordinates (POS += window start - 1); the original read sequence, quality and MM/ML are untouched
    prim2 = list(prim)
    prim2[3] = str(int(prim[3]) + ws - 1)
    rpos = int(prim[3]) - 1 + ws - 1      # 0-based in the canonical sequence
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
        # unit p8 = [s8, e8]: the start boundary lies after the gene end of the preceding unit (s8-1) and the end boundary after this unit's gene end (e8)
        if not (lo_r < s8 - 2 and e8 < hi_r):
            continue
        qs_, bs = boundary(A, m2, s8 - 1, "start")
        qe_, be = boundary(A, m2, e8, "end")
        if qs_ is None or qe_ is None or qe_ < qs_:
            continue
        useq = seq[qs_:qe_ + 1]
        # coordinates on the original read, as in the fastq
        if flag & 16:
            f0, f1 = L0 - (qe_ + 1), L0 - qs_
        else:
            f0, f1 = qs_, qe_ + 1
        mid = s8 + 1000                     # the middle of this unit's spacer, 1-based in the canonical sequence
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
    """items = [(g6, set of spacer variants)] (only units with called variants) -> {d: mean Jaccard against stored number g6+d}"""
    S = {}
    for d in range(-SHIFTMAX, SHIFTMAX + 1):
        v = [jac(t, spacer(P[g + d])) for g, t in items if g is not None and g + d in P]
        if len(v) >= max(2, len(items) // 2):
            S[d] = sum(v) / len(v)
    return S


def shift_pos(u6, pos1, d):
    """canonical coordinate pos1 (1-based) moved by d G6 units (the relative position inside the unit is kept)"""
    st = [x[0] for x in u6]
    i = bisect.bisect_right(st, pos1) - 1
    i = max(i, 0)
    j = i + d
    if j < 0 or j >= len(u6):
        return None
    return pos1 + (u6[j][0] - u6[i][0])


def relabel(res, bd):
    """A read that does not attach to the stored position even after realignment : the cut read coordinates and variants are kept, and only the unit number is moved by the stored difference bd.
    The SAM record keeps the actual alignment and is marked with the ZR tag (the POS of that read is bd units off the stored numbering)"""
    sam, units, foot = res
    units2 = [(us, mt[:7] + [mt[7] + bd, (mt[8] + bd) if mt[8] is not None else None] + mt[9:]) for us, mt in units]
    return (sam + TAB + "ZR:Z:relabeled_d=%+d" % bd, units2, foot)


PBCFG = ("AOJV_PB",)       # this cfg does not truncate the read name, the first 16 characters of a PacBio movie name being shared
def rkey(rid, cfg):
    return rid if cfg in PBCFG else rid[:16]


def work(arg):
    smp, cfg = arg                      # an individual can have several cfgs
    wd = os.path.join(D, "tmp", "%s.%s" % (smp, cfg))
    os.makedirs(wd, exist_ok=True)
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
                    FQ.setdefault(rkey(rid, cfg), (rid, rec, lab))
                    rec = []
    RES = {}                  # r16 -> (sam, units, foot)
    NOTE = {}                 # r16 -> the realignment record
    # -- first pass: the window is the canonical span of the stored G6 units, clipped into the canonical sequence, plus or minus PADBP
    for r16, (slot, P) in sorted(sv.items()):
        if r16 not in FQ:
            out_r.append([smp, r16, "-", slot, cfg, "-", "NOT_IN_FASTQ"])
            continue
        rid, rec, lab = FQ[r16]
        A = ARR[(smp, slot)]
        u6 = U6[(smp, slot)]
        p6 = dict((p, (s, e)) for s, e, p in u6)
        inr = sorted(p for p in p6 if min(P) <= p <= max(P))      # the units of the stored range that exist in the canonical sequence; a stored index beyond the end of the array is clipped
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
    # -- second pass: where the per-read index disagrees clearly with the stored placement (the same index below 0.6 and a neighbouring index at 0.8 or above),
    #          realign against the whole-footprint window at the stored position, plus or minus 300 bp; slipping by one unit leaves 2 kb of the read over, which scores against it
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
    # -- output
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
    return smp, cfg, out_u, out_r, out_s


def init():
    global CFG, SAVED, ARR, U6, U8
    AGS = samples(os.path.join(D, "AOJV_GIAB.tsv"))
    ACS = samples(os.path.join(D, "AOJV8.tsv"))
    AFS = samples(os.path.join(D, "AOJV5.tsv"))
    SAVED = collections.defaultdict(dict)
    for cfg, fn in (("AOJV3", "AOJV3.tsv"), ("AOJV5", "AOJV5.tsv"), ("AOJV8", "AOJV8.tsv"), ("AOJVG", "AOJV_GIAB.tsv"),
                    ("AOJV_FS", "AOJV_FS.tsv"), ("AOJV_PB", "AOJV_PB.tsv")):   # Fiber-seq
        if not os.path.exists(os.path.join(D, fn)):
            continue                   # a data source that is absent is skipped
        for L in io.open(os.path.join(D, fn), encoding="utf8"):
            r = L.rstrip(NL).split(TAB)
            if len(r) < 5:
                continue
            e = SAVED[(cfg, r[0])].setdefault(r[1], [r[2], {}])
            e[1][int(r[3])] = set() if r[4] == "-" else set(r[4].split(","))
    # all available data sources are kept side by side rather than one being chosen,
    #   except AOJVG, AOJV8, AOJV5 and AOJV3, of which one is chosen by priority so that the same read is never duplicated
    CFG = {}
    for s in set(k[1] for k in SAVED):
        base = "AOJVG" if s in AGS else "AOJV8" if s in ACS else "AOJV5" if s in AFS else "AOJV3"
        CFG[s] = [base] + [c for c in ("AOJV_FS", "AOJV_PB") if (c, s) in SAVED]
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
    U6 = loadmap(os.path.join(D, "MAP_G6C2.tsv"), ("UNIT", "UNIT_LAST"))   # the last complete unit of the array (UNIT_LAST) is used for the window and the G6 index as well; without it a read whose highest stored index is the last unit comes out NO_WINDOW
    U8 = loadmap(os.path.join(D, "MAP_G8.tsv"))


if __name__ == "__main__":
    init()
    smps = [x.strip() for x in io.open(LST, encoding="utf8") if x.strip()]
    smps = [s for s in smps if s in CFG]
    OU = io.open(os.path.join(D, "AO_G8_UNITS.tsv"), "a", encoding="utf8", newline=NL)
    OR = io.open(os.path.join(D, "AO_G8_READS.tsv"), "a", encoding="utf8", newline=NL)
    os.makedirs(os.path.join(D, "sam"), exist_ok=True)
    tasks = [(s, c) for s in smps for c in CFG[s]]
    with multiprocessing.Pool(NP, initializer=init) as pool:
        for smp, cfg, uu, rr, ss in pool.imap_unordered(work, tasks):
            for x in uu:
                OU.write(TAB.join(str(y) for y in x) + NL)
            for x in rr:
                OR.write(TAB.join(str(y) for y in x) + NL)
            io.open(os.path.join(D, "sam", "%s.%s.sam" % (smp, cfg)), "w", encoding="utf8", newline=NL).write(NL.join(ss) + (NL if ss else ""))
            OU.flush()
            OR.flush()
            sys.stderr.write("  %s %s 유닛 %d · read %d%s" % (smp, cfg, len(uu), len(rr), NL))
    OU.close()
    OR.close()
    print("done")
