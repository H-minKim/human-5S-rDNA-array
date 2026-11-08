# -*- coding: utf-8 -*-
# cut the flank from the newly assembled contig and attach it to the canonical (2026-09-15)
#
#   Invariant: the canonical array sequence is not changed by a single character. Things are only attached at both ends.
#
#   Why anchor alignment
#     the new assembly is a different run, so the canonical may not be contained character for character.
#     So the 20kb before and the 20kb after the canonical are aligned to the contigs as queries to locate the position,
#     and the 50kb outside that is used as the flank. The identity is recorded so it can be judged whether it is the same haplotype.
#
#   Output: ~/flank_v12/<canonical file name>   header line with flank length, contig and identity
#          ~/FLANK_REASM_<host>.tsv
#   Arguments: none (environment variables ASM default ~/reasm_v12/out, OUT default ~/flank_v12)
import io
import os
import glob
import socket
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~")
ASM = os.path.expanduser(os.environ.get("ASM", "~/reasm_v12/out"))
OUT = os.path.expanduser(os.environ.get("OUT", "~/flank_v12"))
CMPD = os.path.join(H, "cmp4d_v11")
TMP = os.path.join(H, "frtmp")
os.makedirs(OUT, exist_ok=True)
os.makedirs(TMP, exist_ok=True)
FL = 50000
ANC = 20000
CMPX = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(CMPX.get(c, "N") for c in reversed(s))


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8") if not x.startswith(">")).upper()


def readfa(p):
    d = {}
    nm = None
    buf = []
    for L in io.open(p, encoding="utf8"):
        if L.startswith(">"):
            if nm:
                d[nm] = "".join(buf).upper()
            nm = L[1:].strip().split()[0]
            buf = []
        else:
            buf.append(L.strip())
    if nm:
        d[nm] = "".join(buf).upper()
    return d


SEQF = {}
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 3:
        SEQF.setdefault(r[0], []).append((r[1], r[2]))
# only slots from our hifiasm (asmunits3, hg002asm, tier1ref and cons2seq excluded)
OWN = {}
for smp, lst in SEQF.items():
    for slot, fn in lst:
        if not fn.startswith(("asmunits3__", "hg002asm__", "tier1ref__", "cons2seq__")):
            OWN.setdefault(smp, []).append((slot, fn))

LOG = io.open(os.path.join(H, "FLANK_REASM_%s.tsv" % socket.gethostname()), "w",
              encoding="utf8", newline=NL)
LOG.write(TAB.join(["sample", "slot", "file", "contig", "strand", "id5", "id3",
                    "len5", "len3", "canon_bp", "verdict"]) + NL)


def aln(q, tfa, tag):
    """Aligns the query to the contig collection -> (contig, direction, start, end, identity)"""
    qf = os.path.join(TMP, "q_%s.fa" % tag)
    io.open(qf, "w", newline=NL).write(">q" + NL + q + NL)
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--secondary=no", "-t", "4", tfa, qf],
                       capture_output=True, text=True)
    best = None
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        if len(f) < 12:
            continue
        qs, qe = int(f[2]), int(f[3])
        if (qe - qs) < 0.8 * len(q):
            continue
        nm, bl = int(f[9]), int(f[10])
        v = 100.0 * nm / bl if bl else 0
        if best is None or v > best[4]:
            best = (f[5], f[4], int(f[7]), int(f[8]), v)
    try:
        os.remove(qf)
    except OSError:
        pass
    return best


st = collections.Counter()
for smp in sorted(OWN):
    # hifiasm outputs the same sequence several times, as bp.p_ctg, bp.hap1, bp.hap2 and the copy p_ctg.
    #   With all of them in, the front and back anchors attach to different copies and are misjudged as 'different contigs' (measured 2026-09-15).
    #   Only hap1 and hap2 are used. If absent, only bp.p_ctg.
    cands = sorted(glob.glob(os.path.join(ASM, "%s.bp.hap?.p_ctg.fa" % smp)))
    if not cands:
        cands = [x for x in [os.path.join(ASM, "%s.bp.p_ctg.fa" % smp)] if os.path.exists(x)]
    if not cands:
        st["조립본 없음"] += 1
        continue
    # all contigs into one file
    allfa = os.path.join(TMP, "%s_all.fa" % smp)
    g = io.open(allfa, "w", newline=NL)
    CT = {}
    for c in cands:
        tagn = os.path.basename(c).replace(".fa", "")
        for k, v in readfa(c).items():
            key = tagn + "|" + k
            CT[key] = v
            g.write(">%s%s%s%s" % (key, NL, v, NL))
    g.close()
    for slot, fn in OWN[smp]:
        outp = os.path.join(OUT, fn)
        if os.path.exists(outp):
            continue
        canon = rd(os.path.join(CMPD, fn))
        a5 = aln(canon[:ANC], allfa, "%s%s5" % (smp, slot))
        a3 = aln(canon[-ANC:], allfa, "%s%s3" % (smp, slot))
        if not a5 or not a3 or a5[0] != a3[0] or a5[1] != a3[1]:
            st["앵커 불일치"] += 1
            LOG.write(TAB.join([smp, slot, fn, "-", "-", "-", "-", "0", "0",
                                str(len(canon)), "ANCHOR_FAIL"]) + NL)
            continue
        ctg, strand = a5[0], a5[1]
        s = CT[ctg]
        if strand == "-":
            s = rc(s)
            L_ = len(CT[ctg])
            p5 = L_ - a5[3]
            p3 = L_ - a3[2]
        else:
            p5 = a5[2]
            p3 = a3[3]
        if p3 <= p5:
            st["좌표 역전"] += 1
            continue
        lf = s[max(0, p5 - FL):p5]
        rf = s[p3:p3 + FL]
        new = lf + canon + rf
        assert new[len(lf):len(lf) + len(canon)] == canon, "정본 손상"
        io.open(outp, "w", newline=NL).write(
            ">%s|flank5=%d|flank3=%d|src=REASM:%s:%s|id5=%.2f|id3=%.2f%s%s%s"
            % (fn.replace(".fa", ""), len(lf), len(rf), ctg, strand, a5[4], a3[4], NL, new, NL))
        st["OK"] += 1
        LOG.write(TAB.join([smp, slot, fn, ctg, strand, "%.2f" % a5[4], "%.2f" % a3[4],
                            str(len(lf)), str(len(rf)), str(len(canon)), "OK"]) + NL)
        LOG.flush()
    try:
        os.remove(allfa)
    except OSError:
        pass
LOG.close()
print("%s  %s" % (socket.gethostname(), dict(st)))
