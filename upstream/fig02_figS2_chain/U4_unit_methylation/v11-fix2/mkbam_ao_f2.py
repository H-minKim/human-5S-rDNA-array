# -*- coding: utf-8 -*-
# ARRAY_ONLY G8 cut SAM (~/aog8/sam/<individual>.<cfg>.sam) -> sorted and indexed BAM (run on the node)
#   SAM records exactly as written by ao_g8cut.py (original read sequence, QUAL and MM/ML untouched, POS already in canonical coordinates)
#   What is added: header line @SQ (canonical name = <individual>_h<slot>, length = canonical length built the same way as ao_g8cut)
#   canonical direction = G8v3 direction (reverse complement if the flip in SUM_G8 is Y), the reference sequence is also saved as ref/<individual>.fa
#   data source in the file name: <individual>.<cfg>.<source>.G8v3.bam
#   one output line = node, cfg, individual, source, number of records, records with QUAL, records with MM, number of OK reads in the READS table
import io
import os
import sys
import glob
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
D = H("~/fix2/run1/aog8")
OUT = H("~/fix2/run1/placed_bam_g8")
os.makedirs(os.path.join(OUT, "ref"), exist_ok=True)
HOST = os.uname()[1]
CMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


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
LAB = {}
for L in io.open(H("~/srcmap_%s.tsv" % HOST), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 4:
        LAB[(r[1], r[2])] = r[3]
NOK = collections.Counter()
for L in io.open(os.path.join(D, "AO_G8_READS.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 7 and r[6] == "OK":
        NOK[(r[0], r[4])] += 1
for sp in sorted(glob.glob(os.path.join(D, "sam", "*.sam"))):
    smp, cfg = os.path.basename(sp)[:-4].rsplit(".", 1)
    lab = LAB.get((cfg, smp), "UNKNOWN")
    recs = [L.rstrip(NL) for L in io.open(sp, encoding="utf8") if L.strip() and not L.startswith("@")]
    # the third renumbering (ZR tag) was reverted (number = alignment position) -> only the marker tag is removed. Sequence, quality, MM/ML and alignment unchanged
    nzr = sum(1 for x in recs if (TAB + "ZR:Z:") in x)
    recs = [TAB.join(f for f in x.split(TAB) if not f.startswith("ZR:Z:")) for x in recs]
    used = sorted(set(x.split(TAB, 3)[2] for x in recs))
    hdr = ["@HD" + TAB + "VN:1.6" + TAB + "SO:unsorted"]
    fa = io.open(os.path.join(OUT, "ref", "%s.%s.G8v3ref.fa" % (smp, cfg)), "w", encoding="utf8", newline=NL)
    bad = False
    for nm in used:
        s_, sl = nm.rsplit("_h", 1)
        if (s_, sl) not in SEQF:
            bad = True
            continue
        a = rd(os.path.join(D, "cmp4d_v7", SEQF[(s_, sl)])).upper()
        a = rc(a) if FLIP.get((s_, sl)) else a
        hdr.append("@SQ" + TAB + "SN:" + nm + TAB + "LN:%d" % len(a))
        fa.write(">" + nm + NL + a + NL)
    fa.close()
    if bad:
        print(TAB.join([HOST, cfg, smp, lab, "NO_REFSEQ"]))
        continue
    hdr.append("@CO" + TAB + "ARRAY_ONLY placed reads re-aligned to haplotype canonical (v7, G8v3 orientation) window; minimap2 -a -x map-ont --secondary=no -y; POS shifted to canonical coordinates; SEQ/QUAL/MM/ML untouched")
    tmp = "/tmp/mkbam_ao_%s.%s.sam" % (smp, cfg)
    io.open(tmp, "w", encoding="utf8", newline=NL).write(NL.join(hdr + recs) + NL)
    out = os.path.join(OUT, "%s.%s.%s.G8v3.bam" % (smp, cfg, lab))
    r1 = subprocess.run(["samtools", "sort", "-@2", "-o", out, tmp], capture_output=True, text=True)
    r2 = subprocess.run(["samtools", "index", out], capture_output=True, text=True) if r1.returncode == 0 else r1
    os.remove(tmp)
    if r1.returncode or r2.returncode:
        print(TAB.join([HOST, cfg, smp, lab, "SORT_FAIL", (r1.stderr or r2.stderr)[:200].replace(NL, " ")]))
        continue
    v = subprocess.run(["samtools", "view", out], capture_output=True, text=True).stdout.splitlines()
    q = sum(1 for x in v if x.split(TAB, 11)[10] != "*")
    m = sum(1 for x in v if (TAB + "MM:Z:") in x or (TAB + "Mm:Z:") in x)
    print(TAB.join([HOST, cfg, smp, lab, str(len(v)), str(q), str(m), str(NOK[(smp, cfg)]), str(len(recs)), "ZR지움%d" % nzr]))
