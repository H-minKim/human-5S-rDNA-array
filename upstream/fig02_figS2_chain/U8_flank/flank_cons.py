# -*- coding: utf-8 -*-
# build flanks from read consensus (2026-09-14)
#
#   material: <sample>_h<slot>.<L|R>.fa extracted by flank_reads.py
#          L = up to just before the canonical start (anchor at the end side), R = from just after the canonical end (anchor at the front side)
#   method: backbone polishing, the longest read is the backbone and the rest are attached with minimap2,
#          then majority vote per position (substitution, deletion, insertion), 2 rounds. mafft is not used.
#   slots are not mixed. Many samples have nearly identical flanks on the two haplotypes (separation test median +0.2%p),
#     so mixing would not show. Only reads placed on the slot are used.
#   length: at most CAP (default 50,000) bp outward from the anchor
#   output: <OUT>/<sample>_h<slot>.<L|R>.cons.fa   depth and length in the header
#   Arguments: flankread folder   (environment variables OUT, CAP, MINDEP)
import io
import os
import re
import sys
import glob
import socket
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
D = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/flankread")
OUT = os.path.expanduser(os.environ.get("OUT", "~/flankcons"))
CAP = int(os.environ.get("CAP", "50000"))
MINDEP = int(os.environ.get("MINDEP", "2"))
os.makedirs(OUT, exist_ok=True)
W = os.path.expanduser("~/flankconstmp")
os.makedirs(W, exist_ok=True)
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def rdfa(p):
    out = []
    nm = None
    buf = []
    for L in io.open(p, encoding="utf8"):
        if L.startswith(">"):
            if nm:
                out.append((nm, "".join(buf)))
            nm = L[1:].strip().split("|")[0]
            buf = []
        else:
            buf.append(L.strip())
    if nm:
        out.append((nm, "".join(buf)))
    return out


def polish(bb, reads, tag):
    """Attach reads to backbone bb and take a majority vote per position. Substitutions, deletions and insertions are applied within the backbone length."""
    bf = os.path.join(W, "bb_%s.fa" % tag)
    io.open(bf, "w", newline=NL).write(">bb" + NL + bb + NL)
    rf = os.path.join(W, "rd_%s.fa" % tag)
    g = io.open(rf, "w", newline=NL)
    for i, s in enumerate(reads):
        g.write(">r%d%s%s%s" % (i, NL, s, NL))
    g.close()
    r = subprocess.run(["minimap2", "-a", "-x", "map-ont", "--secondary=no", "-t", "2", bf, rf],
                       capture_output=True, text=True)
    col = [collections.Counter() for _ in range(len(bb))]
    ins = collections.defaultdict(collections.Counter)
    nrd = 0
    for L in r.stdout.splitlines():
        if L.startswith("@"):
            continue
        f = L.split(TAB)
        if len(f) < 11 or f[2] == "*":
            continue
        if int(f[1]) & 0x900:
            continue
        nrd += 1
        pos = int(f[3]) - 1
        q = 0
        seq = f[9]
        for a, b in CIG.findall(f[5]):
            a = int(a)
            if b == "S" or b == "H":
                if b == "S":
                    q += a
                continue
            if b in "M=X":
                for k in range(a):
                    if pos + k < len(col):
                        col[pos + k][seq[q + k]] += 1
                pos += a
                q += a
            elif b == "I":
                ins[pos][seq[q:q + a]] += 1
                q += a
            elif b in "DN":
                for k in range(a):
                    if pos + k < len(col):
                        col[pos + k]["-"] += 1
                pos += a
    for x in (bf, rf):
        try:
            os.remove(x)
        except OSError:
            pass
    out = []
    for i, c in enumerate(col):
        if not c:
            out.append(bb[i])
        else:
            ch = c.most_common(1)[0][0]
            if ch != "-":
                out.append(ch)
        if i in ins:
            s, n = ins[i].most_common(1)[0]
            if n > (sum(col[i].values()) / 2.0):
                out.append(s)
    return "".join(out), nrd


n_ok = 0
n_skip = 0
LOG = io.open(os.path.join(OUT, "FLANKCONS_%s.tsv" % socket.gethostname()), "w",
              encoding="utf8", newline=NL)
LOG.write(TAB.join(["sample", "slot", "side", "n_reads", "backbone_len", "cons_len", "used"]) + NL)
for p in sorted(glob.glob(os.path.join(D, "*.?.fa"))):
    b = os.path.basename(p)[:-3]
    smp_slot, side = b.rsplit(".", 1)
    smp, slot = smp_slot.rsplit("_h", 1)
    outp = os.path.join(OUT, "%s_h%s.%s.cons.fa" % (smp, slot, side))
    if os.path.exists(outp):
        continue
    lst = rdfa(p)
    if not lst:
        continue
    # only CAP from the anchor: for L the anchor is at the end -> CAP from the back, for R the anchor is at the front -> CAP from the front
    cut = [(n, (s[-CAP:] if side == "L" else s[:CAP])) for n, s in lst]
    cut.sort(key=lambda x: -len(x[1]))
    if len(cut) < MINDEP:
        io.open(outp, "w", newline=NL).write(
            ">%s_h%s_%s|dep=1|len=%d|single%s%s%s" % (smp, slot, side, len(cut[0][1]), NL, cut[0][1], NL))
        LOG.write(TAB.join([smp, slot, side, "1", str(len(cut[0][1])), str(len(cut[0][1])), "single"]) + NL)
        n_skip += 1
        continue
    bb = cut[0][1]
    others = [s for _, s in cut[1:]]
    cons = bb
    used = 0
    for _ in range(2):
        cons, used = polish(cons, others, "%s%s%s" % (smp, slot, side))
    # CAP again so the anchor-side end is not cut
    cons = cons[-CAP:] if side == "L" else cons[:CAP]
    io.open(outp, "w", newline=NL).write(
        ">%s_h%s_%s|dep=%d|len=%d%s%s%s" % (smp, slot, side, len(cut), len(cons), NL, cons, NL))
    LOG.write(TAB.join([smp, slot, side, str(len(cut)), str(len(bb)), str(len(cons)), str(used)]) + NL)
    n_ok += 1
LOG.close()
print("%s  합의 %d · 단일read %d" % (socket.gethostname(), n_ok, n_skip))
