# -*- coding: utf-8 -*-
# HG00133 h2 (N 153→149, 347,971→337,083 bp), HG03270 h1 (N 78→77, 178,074→175,843 bp)
#   Step 3b: canonical-unit repair (revision 1, 2 sites).
# NODE_HOME has no default.
#       to run : export NODE_HOME=/home/<user>

# -*- coding: utf-8 -*-
# Revision 1: source cmp4d_v11, output cmp4d_v11_fix1 as a real copy and never a symbolic link, map g8canonv10, anchors given explicitly.
# Canonical unit v10: a surgical replacement.
#   Starting from v8 (the original assembly), only the 7 broken sites are spliced as [assembly left anchor][read consensus middle][assembly right anchor].
#   Stored orientation: the .fa from cmp4d holds the sequence before flipping, so a window built in map orientation is reverse-complemented before it is written.
import io
import os
import glob
import shutil
import collections

H = os.environ["NODE_HOME"]   # the compute node's home directory
NL = chr(10)
TAB = chr(9)
CMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
SITES = [("HG00133", "2", 119, 129), ("HG03270", "1", 40, 47)]


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


def rdfa(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()


SEQF = {}
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) >= 3:
        SEQF[(f[0], f[1])] = f[2]
FLIP = {}
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_SUM_G8_full.tsv"), encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    FLIP[(r[0], r[1])] = (r[4] == "Y")
U = collections.defaultdict(list)
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] == "UNIT":
        U[(r[0], r[1])].append((int(r[ix["array_pos"]]), int(r[ix["unit_start"]]), int(r[ix["unit_end"]])))
for k in U:
    U[k].sort()

V10 = os.path.join(H, "cmp4d_v11_fix1")
assert not os.path.exists(V10), "already exists - not overwritten"
os.makedirs(V10)
n = 0
for p in glob.glob(os.path.join(H, "cmp4d_v11", "*.fa")):
    b = os.path.basename(p)
    t = os.path.join(V10, b)
    shutil.copyfile(os.path.realpath(p), t)          # a real copy
    n += 1
print("files actually copied from v11 %d" % n)

LOG = io.open(os.path.join(H, "fix1", "CANON_V11_FIX1_PATCH.tsv"), "w", encoding="utf8", newline=NL)
LOG.write(TAB.join(["sample", "slot", "after_unit", "file", "anchor_start", "anchor_end",
                    "asm_win_bp", "splice_win_bp", "hap_bp_old", "hap_bp_new", "flip"]) + NL)
for (smp, sl, ap, apR) in SITES:
    fn = SEQF[(smp, sl)]
    raw = rdfa(os.path.join(H, "cmp4d_v11", fn))
    A = rc(raw) if FLIP.get((smp, sl)) else raw
    v = U[(smp, sl)]
    i = [j for j, x in enumerate(v) if x[0] == ap][0]
    iR = [j for j, x in enumerate(v) if x[0] == apR][0]
    a0, a1 = v[i][1], v[iR][2]
    spl = rdfa(os.path.join(H, "fix1", "splice", "spl_%s_%s_%d.fa" % (smp, sl, ap)))
    # Confirm the anchor sequence matches the assembly character for character; the surgical replacement depends on it.
    assert spl.startswith(A[v[i][1] - 1:v[i][2]]), "front anchor mismatch %s %s" % (smp, sl)
    assert spl.endswith(A[v[iR][1] - 1:v[iR][2]]), "rear anchor mismatch %s %s" % (smp, sl)
    B = A[:a0 - 1] + spl + A[a1:]
    out = rc(B) if FLIP.get((smp, sl)) else B
    t = os.path.join(V10, fn)
    assert not os.path.islink(t)
    io.open(t, "w", newline=NL).write(">%s_splice_v11fix1%s%s%s" % (fn.replace(".fa", ""), NL, out, NL))
    print("  %s h%s unit%d/%d: window %d -> %d bp, total %d -> %d bp"
          % (smp, sl, v[i][0], v[iR][0], a1 - a0 + 1, len(spl), len(A), len(B)))
    LOG.write(TAB.join([smp, sl, str(ap), fn, str(a0), str(a1), str(a1 - a0 + 1),
                        str(len(spl)), str(len(A)), str(len(B)),
                        "Y" if FLIP.get((smp, sl)) else "N"]) + NL)
LOG.close()
real = sum(1 for x in os.listdir(V10) if not os.path.islink(os.path.join(V10, x)))
print("cmp4d_v10: total %d, actual files %d" % (len(os.listdir(V10)), real))
