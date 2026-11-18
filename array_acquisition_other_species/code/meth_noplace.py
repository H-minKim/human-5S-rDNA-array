# -*- coding: utf-8 -*-
# Measure methylation without placement: the 5S methylation state is obtained without placing the reads onto the array.
#
#   Why: the array is tandem, so reads are placed at mapq 0 (69-72% of the placed zebra finch reads).
#        If which unit a read was placed on cannot be trusted, neither can the conclusion.
#        But the methylation calls are carried on the read itself, in MM/ML, and need no reference.
#        All that is needed is whether this stretch of the read is 5S and which kind of 5S,
#        and that follows directly from a sequence comparison. Neither array coordinates nor mapq are used.
#
#   method
#     1) The read is cut into windows of WIN bp, using no alignment and no coordinates.
#     2) Each window is aligned only against the reference unit types and labelled with the highest-scoring type, using the context version.
#     3) The methylation fraction over the C modifications in that window is computed from MM/ML, counting mod_qual >= 0.5 as methylated.
#        The context is decided from the read sequence directly: the next base G gives CG, the one after that G gives CHG, otherwise CHH.
#        (A plant calls all three contexts, so pooling them dilutes a CG value of 0.84 to 0.5.)
#        Where the calling model is CpG-only, as for the great apes, the zebra finch and the zebrafish, every C in MM is a CpG.
#          Splitting by context from the read sequence then splits only the read errors (measured in gorilla: 18% leaked into CHH).
#          So the kind of model is decided from the data, by the fraction of MM positions whose next base is G,
#          and above 0.75 everything is taken as CG.
#     4) The fractions are collected per type and per read.
#   The BAM is used only as a container for MM/ML; the alignment coordinates, the mapq and the unit index are never read.
#
#   Outputs: <out>/<species>__<sample>.noplace_window.tsv, one row per window,
#          and <out>/<species>__<sample>.noplace_read.tsv, one row per read and type.
#   Arguments: species key, sample, bam, reference unit types .fa, output folder, [WIN = 500], [MINID = 0.70], [CTX = 120].
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import csv
import subprocess
import collections
import pysam

NL = chr(10)
TAB = chr(9)
SP, SAMPLE, BAM, REFT, OUT = sys.argv[1:6]
WIN = int(sys.argv[6]) if len(sys.argv) > 6 else 500
MINID = float(sys.argv[7]) if len(sys.argv) > 7 else 0.70
CTX = int(sys.argv[8]) if len(sys.argv) > 8 else 120
# The calling threshold: modkit discards the least confident 10% by default and applies 0.5 to what remains.
#   Applying 0.5 here directly would put every borderline call into the methylated class and inflate a weakly methylated context such as CHH.
#   To match the same ruler the threshold is raised, and several values are reported together as a sensitivity.
MODTHR = float(sys.argv[9]) if len(sys.argv) > 9 else 0.5
os.makedirs(OUT, exist_ok=True)
TMP = os.path.join(OUT, "_t_%s_%s" % (SP, SAMPLE))
os.makedirs(TMP, exist_ok=True)

# reference unit type -> its context version, the array being tandem so that a window straddles a unit boundary
refs = collections.OrderedDict()
name = None
for L in io.open(REFT, encoding="utf8"):
    L = L.rstrip(NL)
    if L.startswith(">"):
        name = L[1:].split()[0]
        refs[name] = []
    elif name:
        refs[name].append(L.strip())
refs = collections.OrderedDict((k, "".join(v).upper()) for k, v in refs.items())
cp = os.path.join(TMP, "ctx.fa")
with io.open(cp, "w", encoding="utf8", newline=NL) as fh:
    for k, s in refs.items():
        r = (s * 4)[:max(len(s) * 2, WIN + 2 * CTX)]   # the context version, concatenated tandem
        fh.write(">%s%s%s%s" % (k, NL, r, NL))

bam = pysam.AlignmentFile(BAM, "rb", check_sq=False)
wf = io.open(os.path.join(TMP, "win.fa"), "w", encoding="utf8", newline=NL)
META = []
MOD = {}
nread = 0
SEQ = {}
for rd in bam.fetch(until_eof=True):
    if rd.is_secondary or rd.is_supplementary:
        continue
    q = rd.query_sequence
    if not q:
        continue
    try:
        mb = rd.modified_bases
    except (ValueError, KeyError):
        mb = None
    if not mb:
        continue
    nread += 1
    # CpG 5mC: C+m and C+h are collected from the (canonical, strand, code) key
    call = {}
    for k, v in mb.items():
        code = k[2]
        if code not in ("m", "h", 27551):
            continue
        st = k[1]          # 0 = the read strand; 1 = the opposite strand, where the read base at that position is G
        for pos, ql in v:
            p = call.get(pos, [0.0, 0.0, st])
            p[2] = st
            if code == "m":
                p[0] = max(p[0], ql / 255.0)
            else:
                p[1] = max(p[1], ql / 255.0)
            call[pos] = p
    if not call:
        continue
    MOD[rd.query_name] = call
    SEQ[rd.query_name] = q
    for off in range(0, len(q) - WIN + 1, WIN):
        META.append((rd.query_name, off))
        wf.write(">w%d%s%s%s" % (len(META) - 1, NL, q[off:off + WIN], NL))
wf.close()
bam.close()
best = {}
rm = subprocess.run("minimap2 -cx map-ont -N 5 -p 0.6 -t 8 %s %s 2>/dev/null"
                    % (cp, os.path.join(TMP, "win.fa")), shell=True, capture_output=True, text=True)
for L in rm.stdout.splitlines():
    f = L.split(TAB)
    if len(f) < 12:
        continue
    nm, tgt = f[0], f[5]
    nmatch, blk = int(f[9]), int(f[10])
    sc = 0
    for x in f[12:]:
        if x.startswith("AS:i:"):
            sc = int(x[5:])
    if nm not in best or sc > best[nm][0]:
        best[nm] = (sc, tgt, nmatch / float(blk) if blk else 0, int(f[2]), int(f[3]))

NG = [0, 0]
for rid, call in MOD.items():
    sq = SEQ.get(rid, "")
    for p, v in call.items():
        NG[0] += 1
        if v[2] == 1:      # on the opposite strand, a G in the read preceded by a C is a CpG
            if p - 1 >= 0 and sq[p - 1] == "C":
                NG[1] += 1
        elif p + 1 < len(sq) and sq[p + 1] == "G":
            NG[1] += 1
CPGONLY = (NG[1] / float(NG[0])) >= 0.75 if NG[0] else False


COMP2 = {"A": "T", "C": "G", "G": "C", "T": "A"}


def ctx_of(seq, p, st):
    if CPGONLY:
        return "CG"
    # deciding the context without placement: only the read sequence is used, read backwards on the opposite strand
    if st == 1:
        a = COMP2.get(seq[p - 1], "N") if p - 1 >= 0 else "N"
        b = COMP2.get(seq[p - 2], "N") if p - 2 >= 0 else "N"
    else:
        a = seq[p + 1] if p + 1 < len(seq) else "N"
        b = seq[p + 2] if p + 2 < len(seq) else "N"
    if a == "G":
        return "CG"
    if b == "G":
        return "CHG"
    return "CHH"


HW = ["species_key", "sample", "read_id", "win_offset", "ref_type", "context", "identity", "aln_frac",
      "n_c", "n_5mC", "frac_5mC", "n_5hmC"]
ow = io.open(os.path.join(OUT, "%s__%s.noplace_window.tsv" % (SP, SAMPLE)), "w",
             encoding="utf8", newline=NL)
ow.write(TAB.join(HW) + NL)
byread = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0]))  # read -> (type, context)
tot = collections.Counter()
for i, (rid, off) in enumerate(META):
    b = best.get("w%d" % i)
    if not b or b[2] < MINID:
        tot["unassigned"] += 1
        continue
    sc, tgt, idt, qs, qe = b
    call = MOD.get(rid, {})
    sq = SEQ.get(rid, "")
    acc = collections.defaultdict(lambda: [0, 0, 0])
    for p, (pm, ph, st) in call.items():
        if off <= p < off + WIN:
            c = acc[ctx_of(sq, p, st)]
            if max(pm, ph) < MODTHR and max(pm, ph) > 1.0 - MODTHR:
                continue          # borderline calls are discarded, which corresponds to modkit's removal of the least confident
            c[0] += 1
            if pm >= MODTHR:
                c[1] += 1
            if ph >= MODTHR:
                c[2] += 1
    if not acc:
        continue
    for c, (n, m, h) in sorted(acc.items()):
        if n == 0:          # a context left with only borderline calls is skipped, which happens as the threshold is raised
            continue
        ow.write(TAB.join([SP, SAMPLE, rid, str(off), tgt, c, "%.3f" % idt,
                           "%.3f" % ((qe - qs) / float(WIN)), str(n), str(m),
                           "%.3f" % (m / float(n)), str(h)]) + NL)
        a = byread[rid][(tgt, c)]
        a[0] += n
        a[1] += m
        a[2] += h
    tot[tgt] += 1
ow.close()
HR = ["species_key", "sample", "read_id", "ref_type", "context", "n_c", "n_5mC", "frac_5mC", "n_5hmC"]
orr = io.open(os.path.join(OUT, "%s__%s.noplace_read.tsv" % (SP, SAMPLE)), "w",
              encoding="utf8", newline=NL)
orr.write(TAB.join(HR) + NL)
cnt = collections.Counter()
agg = collections.defaultdict(lambda: [0, 0, 0, 0])
for rid, d in byread.items():
    for (tgt, c), (n, m, h) in d.items():
        orr.write(TAB.join([SP, SAMPLE, rid, tgt, c, str(n), str(m),
                            "%.3f" % (m / float(n)) if n else "-", str(h)]) + NL)
        a = agg[(tgt, c)]
        a[0] += 1
        a[1] += n
        a[2] += m
        a[3] += h
orr.close()
for f in os.listdir(TMP):
    os.remove(os.path.join(TMP, f))
os.rmdir(TMP)
print("%-22s %-14s threshold %.2f, model %s(next base G among MM positions %.3f), read %5d, windows %7d, unassigned %6d"
      % (SP, SAMPLE, MODTHR, "CpG only" if CPGONLY else "전문맥",
         NG[1] / float(NG[0]) if NG[0] else 0, nread, len(META), tot["unassigned"]))
for t in sorted(agg):
    a = agg[t]
    print("   %-46s %-4s read %5d · C %9d · 5mC %.3f · 5hmC %.4f"
          % (t[0][:46], t[1], a[0], a[1], a[2] / float(a[1]) if a[1] else 0,
             a[3] / float(a[1]) if a[1] else 0))
