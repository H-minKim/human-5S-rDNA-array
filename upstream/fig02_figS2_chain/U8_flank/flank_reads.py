# -*- coding: utf-8 -*-
# Extract flank sequences from junction reads + test whether the two haplotypes separate by flank alone (2026-09-14)
#
#   Reference frame : rc(distal flank 5000) + [canonical] + rc(proximal flank 5000)
#   Reads were not trimmed, the part outside the reference is still in the SAM as soft-clip (measured median 50kb, max 345kb)
#
#   What is extracted (anchored at the canonical boundary, going outward)
#     left (distal)    : for reads touching reference position 1, [left soft-clip + part aligned to reference 1..5000]
#     right (proximal) : for reads touching the reference end, [part aligned to reference (L-4999)..L + right soft-clip]
#   Anchor = just outside the first/last canonical base. Coordinates are counted outward from here (0 = right next to the canonical).
#
#   Output : <OUT>/<sample>_h<slot>.<L|R>.fa   (one line per read, length in the header)
#          <OUT>/FLANKREAD_<host>.tsv       (sample, slot, side, read, length)
#   Arguments : SAM folder [SAM folder ...]   (environment variable OUT, default ~/flankread)
import io
import os
import re
import sys
import glob
import socket

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~")
OUT = os.path.expanduser(os.environ.get("OUT", "~/flankread"))
os.makedirs(OUT, exist_ok=True)
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
FL = 5000                      # flank length contained in the reference
EDGE = 200                     # slack for counting a read as having 'touched' the reference end
CMPX = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(CMPX.get(c, "N") for c in reversed(s))


LOG = io.open(os.path.join(OUT, "FLANKREAD_%s.tsv" % socket.gethostname()), "w",
              encoding="utf8", newline=NL)
LOG.write(TAB.join(["sample", "slot", "side", "read", "len", "clip", "aln_in_flank"]) + NL)
BUF = {}
for d in sys.argv[1:]:
    for p in sorted(glob.glob(os.path.join(os.path.expanduser(d), "*.placed.sam"))):
        smp = os.path.basename(p).split(".")[0]
        LN = {}
        for L in io.open(p, encoding="utf8"):
            if L.startswith("@"):
                if L.startswith("@SQ"):
                    f = dict(x.split(":", 1) for x in L.rstrip(NL).split(TAB)[1:])
                    LN[f["SN"]] = int(f["LN"])
                continue
            f = L.rstrip(NL).split(TAB)
            if len(f) < 11 or f[2] == "*":
                continue
            flag = int(f[1])
            if flag & 0x900:                      # skip supplementary and secondary alignments (avoids flank contamination)
                continue
            ref = LN.get(f[2])
            if not ref:
                continue
            slot = f[2].rsplit("_h", 1)[1]
            cig = [(int(a), b) for a, b in CIG.findall(f[5])]
            if not cig:
                continue
            pos = int(f[3])
            seq = f[9]
            sl = cig[0][0] if cig[0][1] == "S" else 0
            sr = cig[-1][0] if cig[-1][1] == "S" else 0
            alen = sum(a for a, b in cig if b in "MDN=X")
            end = pos + alen - 1
            # -- left (distal) : touches reference position 1
            if pos <= EDGE:
                # count the read length consumed up to reference FL(=5000)
                q = sl
                r = pos - 1
                take = None
                for a, b in cig:
                    if b == "S":
                        continue
                    if b in "M=X":
                        if r + a >= FL:
                            q += (FL - r)
                            take = q
                            break
                        q += a
                        r += a
                    elif b == "I":
                        q += a
                    elif b in "DN":
                        if r + a >= FL:
                            take = q
                            break
                        r += a
                if take is None:
                    take = len(seq)
                s = seq[:take].upper()
                if len(s) >= 100:
                    BUF.setdefault((smp, slot, "L"), []).append((f[0], s))
                    LOG.write(TAB.join([smp, slot, "L", f[0], str(len(s)), str(sl),
                                        str(max(0, take - sl))]) + NL)
            # -- right (proximal) : touches the reference end
            if end >= ref - EDGE:
                start_ref = ref - FL + 1
                q = sl
                r = pos - 1
                begin = None
                for a, b in cig:
                    if b == "S":
                        continue
                    if b in "M=X":
                        if r + a >= start_ref - 1 and begin is None:
                            begin = q + max(0, (start_ref - 1) - r)
                        q += a
                        r += a
                    elif b == "I":
                        q += a
                    elif b in "DN":
                        if r + a >= start_ref - 1 and begin is None:
                            begin = q
                        r += a
                if begin is None:
                    begin = sl
                s = seq[begin:].upper()
                if len(s) >= 100:
                    BUF.setdefault((smp, slot, "R"), []).append((f[0], s))
                    LOG.write(TAB.join([smp, slot, "R", f[0], str(len(s)), str(sr),
                                        str(max(0, len(s) - sr))]) + NL)
LOG.close()
n = 0
for (smp, slot, side), lst in sorted(BUF.items()):
    lst.sort(key=lambda x: -len(x[1]))
    g = io.open(os.path.join(OUT, "%s_h%s.%s.fa" % (smp, slot, side)), "w",
                encoding="utf8", newline=NL)
    for rid, s in lst:
        g.write(">%s|len=%d%s%s%s" % (rid, len(s), NL, s, NL))
    g.close()
    n += 1
print("%s  flank read 묶음 %d 개 · read %d" % (socket.gethostname(), n,
                                              sum(len(v) for v in BUF.values())))
