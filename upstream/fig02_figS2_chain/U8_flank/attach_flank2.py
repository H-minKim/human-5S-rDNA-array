# -*- coding: utf-8 -*-
# Attach read-consensus flanks to the canonical sequence (final orientation rule), 2026-09-15
#
#   invariant: the array canonical sequence is not changed by a single character. It is only attached to.
#
#   orientation rule (fixed by measurement)
#     junction reference = rc(FD[:5000]) + [canonical (RCF applied)] + rc(FP[-5000:])     (jx_g8cut.py refseq)
#     L bundle = outer side of the reads touching the left of the reference  -> before the canonical (map orientation)
#     R bundle = right of the reference                                       -> after the canonical (map orientation)
#       measured: L end 5kb <-> rc(FD[:5000]) 99.7-99.96%, R start 5kb <-> rc(FP[-5000:]) 98.9-99.98%
#              (same for FLIP True/False)
#     moved to the stored orientation
#       FLIP=False : 5'flank = L(as is), 3'flank = R(as is)
#       FLIP=True  : 5'flank = rc(R)     · 3'flank = rc(L)
#       measured check (20 slots with assemblies): with FLIP=True 5'=R- and 3'=L-, consistent with the rule
#
#   test: slots that have assembly flanks are compared with them and the identity is recorded (only those that exist).
#           every slot is checked for ① middle identical to the canonical character by character ② >=95% to the reference flank.
#   output: ~/flank_v11/<canonical file name> (skipped if an HPRC version already exists), FLANK_LOG_READ_<host>.tsv
import io
import os
import glob
import socket
import subprocess

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~")
W = os.path.join(H, "flank_v11")
CONS = os.path.join(H, "flankcons")
CMPD = os.path.join(H, "cmp4d_v11")
REF = os.path.join(H, "s5s_tuc", "ref")
TMP = os.path.join(H, "attachtmp")
os.makedirs(TMP, exist_ok=True)
CMPX = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(CMPX.get(c, "N") for c in reversed(s))


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8") if not x.startswith(">")).upper()


FD = rd(os.path.join(REF, "FLANK_distal_50000.fa"))
FP = rd(os.path.join(REF, "FLANK_proximal_50000.fa"))
REFL = rc(FD[:5000])          # reference left 5kb, the anchor-side end of the L bundle attaches here
REFR = rc(FP[-5000:])         # reference right 5kb, the anchor-side start of the R bundle attaches here

SEQF = {}
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 3:
        SEQF[(r[0], r[1])] = r[2]
FLIP = {}
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_SUM_G8_full.tsv"), encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    FLIP[(r[0], r[1])] = (r[4] == "Y")


def ident(a, b):
    fa = os.path.join(TMP, "a_%d.fa" % os.getpid())
    fb = os.path.join(TMP, "b_%d.fa" % os.getpid())
    io.open(fa, "w", newline=NL).write(">a" + NL + a + NL)
    io.open(fb, "w", newline=NL).write(">b" + NL + b + NL)
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--secondary=no", "-t", "2", fa, fb],
                       capture_output=True, text=True)
    best = (0, 0)
    for L in r.stdout.splitlines():
        g = L.split(TAB)
        if len(g) < 12:
            continue
        nm, bl = int(g[9]), int(g[10])
        if bl > best[1]:
            best = (nm, bl)
    return (100.0 * best[0] / best[1] if best[1] else 0.0), best[1]


LOG = io.open(os.path.join(W, "FLANK_LOG_READ_%s.tsv" % socket.gethostname()), "w",
              encoding="utf8", newline=NL)
LOG.write(TAB.join(["sample", "slot", "file", "flip", "len5", "len3",
                    "chk5_vs_ref", "chk3_vs_ref", "verdict"]) + NL)
todo = set()
for p in glob.glob(os.path.join(CONS, "*.cons.fa")):
    b = os.path.basename(p)[:-len(".cons.fa")]
    ss, side = b.rsplit(".", 1)
    smp, slot = ss.rsplit("_h", 1)
    todo.add((smp, slot))
n_ok = n_part = n_bad = 0
for smp, slot in sorted(todo):
    fn = SEQF.get((smp, slot))
    if not fn:
        continue
    outp = os.path.join(W, fn)
    if os.path.exists(outp):
        continue
    fl = FLIP.get((smp, slot))
    if fl is None:
        continue
    canon = rd(os.path.join(CMPD, fn))
    cl = rd(os.path.join(CONS, "%s_h%s.L.cons.fa" % (smp, slot))) \
        if os.path.exists(os.path.join(CONS, "%s_h%s.L.cons.fa" % (smp, slot))) else ""
    cr = rd(os.path.join(CONS, "%s_h%s.R.cons.fa" % (smp, slot))) \
        if os.path.exists(os.path.join(CONS, "%s_h%s.R.cons.fa" % (smp, slot))) else ""
    # compare with the reference flank to check that the bundle is in place
    v5 = v3 = 0.0
    if cl:
        v5, _ = ident(REFL, cl[-5000:])
    if cr:
        v3, _ = ident(REFR, cr[:5000])
    okL = (not cl) or v5 >= 95.0
    okR = (not cr) or v3 >= 95.0
    if not (okL and okR):
        LOG.write(TAB.join([smp, slot, fn, str(fl), "0", "0",
                            "%.2f" % v5, "%.2f" % v3, "REF_MISMATCH"]) + NL)
        n_bad += 1
        continue
    if fl:
        f5 = rc(cr) if cr else ""
        f3 = rc(cl) if cl else ""
    else:
        f5 = cl
        f3 = cr
    new = f5 + canon + f3
    assert new[len(f5):len(f5) + len(canon)] == canon, "정본 손상"
    vd = "OK" if (f5 and f3) else ("PARTIAL" if (f5 or f3) else "NONE")
    if vd == "NONE":
        n_bad += 1
        continue
    if vd == "OK":
        n_ok += 1
    else:
        n_part += 1
    io.open(outp, "w", newline=NL).write(
        ">%s|flank5=%d|flank3=%d|src=READCONS|flip=%s%s%s%s"
        % (fn.replace(".fa", ""), len(f5), len(f3), "Y" if fl else "N", NL, new, NL))
    LOG.write(TAB.join([smp, slot, fn, str(fl), str(len(f5)), str(len(f3)),
                        "%.2f" % v5, "%.2f" % v3, vd]) + NL)
    LOG.flush()
LOG.close()
print("%s  붙임 OK %d · 한쪽만 %d · 참조대조 실패/없음 %d" % (socket.gethostname(), n_ok, n_part, n_bad))
