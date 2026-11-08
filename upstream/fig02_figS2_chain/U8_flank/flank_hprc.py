# -*- coding: utf-8 -*-
# flank restoration (HPRC assemblies): cut 50kb on each side of the canonical from the contig and attach it (2026-09-14)
#
#   Invariant: the canonical array sequence is not changed by a single character. Things are only attached before and after.
#   The slot number and the assembly haplotype number are not the same, both haplotype contigs are searched and the one that
#     contains the canonical character for character is chosen (reverse direction also checked). If not found, it is only recorded and not touched.
#
#   Outputs: ~/flank_v11/<canonical file name>      = [5' flank][canonical][3' flank]
#          ~/flank_v11/FLANK_LOG.tsv     = per slot source, direction, lengths before and after, canonical md5 (for checking)
#   Arguments: individual [individual ...]
import io
import os
import sys
import glob
import hashlib
import subprocess

H = os.path.expanduser("~")
NL = chr(10)
TAB = chr(9)
FL = 50000
W = os.path.join(H, "flank_v11")
os.makedirs(W, exist_ok=True)
TMP = os.path.join(H, "flanktmp")
os.makedirs(TMP, exist_ok=True)
IDX = os.path.join(H, "telomere_pilot/data/assemblies_release2_v1.0.index.csv")
CMPD = os.path.join(H, "cmp4d_v11")
CMPX = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(CMPX.get(c, "N") for c in reversed(s))


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8") if not x.startswith(">")).upper()


def md5(s):
    return hashlib.md5(s.encode()).hexdigest()


URL = {}
f = io.open(IDX, encoding="utf8")
col = dict((c, k) for k, c in enumerate(f.readline().rstrip(NL).split(",")))
for L in f:
    r = L.rstrip(NL).split(",")
    if len(r) <= col["assembly"]:
        continue
    u = r[col["assembly"]].replace("s3://human-pangenomics/",
                                   "https://s3-us-west-2.amazonaws.com/human-pangenomics/")
    URL.setdefault(r[col["sample_id"]], []).append((r[col["haplotype"]], r[col["assembly_name"]], u))

CTG = {}
for p in glob.glob(os.path.join(H, "asmunits3", "asm_units_*.tsv")):
    for L in io.open(p, encoding="utf8"):
        r = L.rstrip(NL).split(TAB)
        if len(r) >= 9 and r[8] == "OK":
            CTG[(r[0], r[1])] = r[3]

SEQF = {}
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 3:
        SEQF.setdefault(r[0], []).append((r[1], r[2]))

LOGP = os.path.join(W, "FLANK_LOG_%s.tsv" % os.environ.get("TAG", "x"))
LOG = io.open(LOGP, "a", encoding="utf8", newline=NL)


def fetch(u, ctg, tag):
    fai = os.path.join(TMP, tag + ".fai")
    gzi = os.path.join(TMP, tag + ".gzi")
    try:
        subprocess.run(["curl", "-sfL", "-o", fai, u + ".fai"], check=True, timeout=900)
        subprocess.run(["curl", "-sfL", "-o", gzi, u + ".gzi"], check=True, timeout=900)
    except Exception:
        return None
    SZ = {}
    for L in io.open(fai, encoding="utf8"):
        r = L.rstrip(NL).split(TAB)
        if len(r) >= 2:
            SZ[r[0]] = int(r[1])
    if ctg not in SZ:
        return None
    rr = subprocess.run(["samtools", "faidx", "--fai-idx", fai, "--gzi-idx", gzi,
                         u, "%s:1-%d" % (ctg, SZ[ctg])], capture_output=True, text=True)
    for x in (fai, gzi):
        try:
            os.remove(x)
        except OSError:
            pass
    if rr.returncode != 0:
        return None
    return "".join(rr.stdout.split(NL)[1:]).upper()


for smp in sys.argv[1:]:
    if smp not in URL or smp not in SEQF:
        print("  %s : 인덱스/정본 없음" % smp)
        continue
    # fetch this individual's haplotype contig sequences in advance (so they are not fetched again for each slot)
    CACHE = {}
    for hap, anm, u in URL[smp]:
        c = CTG.get((smp, hap))
        if not c:
            continue
        s = fetch(u, c, "%s_%s" % (smp, hap))
        if s:
            CACHE[hap] = (c, s)
    if not CACHE:
        print("  %s : contig 못 받음" % smp)
        continue
    for slot, fn in SEQF[smp]:
        outp = os.path.join(W, fn)
        if os.path.exists(outp):
            continue
        canon = rd(os.path.join(CMPD, fn))
        hit = None
        for hap, (c, s) in sorted(CACHE.items()):
            i = s.find(canon)
            if i >= 0:
                hit = (hap, c, s, i, "+")
                break
            j = s.find(rc(canon))
            if j >= 0:
                hit = (hap, c, s, j, "-")
                break
        if not hit:
            print(" %s slot%s : not found inside the contig (%d bp)" % (smp, slot, len(canon)))
            LOG.write(TAB.join([smp, slot, fn, "-", "-", "-", "0", "0",
                                md5(canon), "NOT_FOUND"]) + NL)
            LOG.flush()
            continue
        hap, c, s, i, strand = hit
        # orient the flank to the direction in which the canonical is stored
        if strand == "+":
            lf, rf = s[max(0, i - FL):i], s[i + len(canon):i + len(canon) + FL]
        else:
            lf, rf = rc(s[i + len(canon):i + len(canon) + FL]), rc(s[max(0, i - FL):i])
        new = lf + canon + rf
        assert new[len(lf):len(lf) + len(canon)] == canon, "정본 손상"
        io.open(outp, "w", newline=NL).write(
            ">%s|flank5=%d|flank3=%d|src=HPRC:%s:%s%s%s%s"
            % (fn.replace(".fa", ""), len(lf), len(rf), c, strand, NL, new, NL))
        print("  %s 슬롯%s : 정본 %d bp · 앞 %d · 뒤 %d (hap%s %s)"
              % (smp, slot, len(canon), len(lf), len(rf), hap, strand))
        LOG.write(TAB.join([smp, slot, fn, c, hap, strand, str(len(lf)), str(len(rf)),
                            md5(canon), "OK"]) + NL)
        LOG.flush()
LOG.close()
