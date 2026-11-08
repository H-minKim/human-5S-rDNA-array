# -*- coding: utf-8 -*-
# flank restoration round 2, retry of failures plus the anchor method (2026-09-14)
#
#   two classes not done in round 1
#     (a) contig not downloaded: network, simply download again
#     (b) canonical not found in the contig: repaired slots (the middle is a read consensus, so it is not in the assembly) and the like
#          -> search the contig only for the 5kb at both ends of the canonical and use what lies outside as the flank
#             (even if the middle differs, both ends are as in the assembly, so they serve as anchors)
#
#   Invariant: the canonical array sequence is not changed by a single letter. Flanks are only attached.
#   Arguments: sample [sample ...]
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
ANC = 5000
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

LOG = io.open(os.path.join(W, "FLANK_LOG2_%s.tsv" % os.environ.get("TAG", "x")), "a",
              encoding="utf8", newline=NL)


def fetch(u, ctg, tag):
    fai = os.path.join(TMP, tag + ".fai")
    gzi = os.path.join(TMP, tag + ".gzi")
    for _ in range(3):
        try:
            subprocess.run(["curl", "-sfL", "-o", fai, u + ".fai"], check=True, timeout=1200)
            subprocess.run(["curl", "-sfL", "-o", gzi, u + ".gzi"], check=True, timeout=1200)
            break
        except Exception:
            continue
    else:
        return None
    SZ = {}
    for L in io.open(fai, encoding="utf8"):
        r = L.rstrip(NL).split(TAB)
        if len(r) >= 2:
            SZ[r[0]] = int(r[1])
    if ctg not in SZ:
        return None
    for _ in range(3):
        rr = subprocess.run(["samtools", "faidx", "--fai-idx", fai, "--gzi-idx", gzi,
                             u, "%s:1-%d" % (ctg, SZ[ctg])], capture_output=True, text=True)
        if rr.returncode == 0 and len(rr.stdout) > 10000:
            break
    else:
        return None
    for x in (fai, gzi):
        try:
            os.remove(x)
        except OSError:
            pass
    return "".join(rr.stdout.split(NL)[1:]).upper()


def locate(s, canon):
    """Whole canonical first, then anchors at both ends (5kb) if that fails. Returns (start, end, strand, method)."""
    i = s.find(canon)
    if i >= 0:
        return i, i + len(canon), "+", "FULL"
    j = s.find(rc(canon))
    if j >= 0:
        return j, j + len(canon), "-", "FULL"
    for strand, q in (("+", canon), ("-", rc(canon))):
        a = s.find(q[:ANC])
        b = s.rfind(q[-ANC:])
        if a >= 0 and b >= 0 and b + ANC > a:
            return a, b + ANC, strand, "ANCHOR"
    return None


for smp in sys.argv[1:]:
    if smp not in URL or smp not in SEQF:
        print("  %s : 인덱스/정본 없음" % smp)
        continue
    need = [(slot, fn) for slot, fn in SEQF[smp] if not os.path.exists(os.path.join(W, fn))]
    if not need:
        continue
    CACHE = {}
    for hap, anm, u in URL[smp]:
        c = CTG.get((smp, hap))
        if not c:
            continue
        s = fetch(u, c, "%s_%s_r2" % (smp, hap))
        if s:
            CACHE[hap] = (c, s)
    if not CACHE:
        print("  %s : contig 못 받음(재시도 후에도)" % smp)
        continue
    for slot, fn in need:
        canon = rd(os.path.join(CMPD, fn))
        hit = None
        for hap, (c, s) in sorted(CACHE.items()):
            r = locate(s, canon)
            if r:
                hit = (hap, c, s) + r
                break
        if not hit:
            print(" %s slot%s : not found by the anchor either" % (smp, slot))
            LOG.write(TAB.join([smp, slot, fn, "-", "-", "-", "0", "0", md5(canon), "NOT_FOUND"]) + NL)
            LOG.flush()
            continue
        hap, c, s, a0, a1, strand, how = hit
        if strand == "+":
            lf, rf = s[max(0, a0 - FL):a0], s[a1:a1 + FL]
        else:
            lf, rf = rc(s[a1:a1 + FL]), rc(s[max(0, a0 - FL):a0])
        new = lf + canon + rf
        assert new[len(lf):len(lf) + len(canon)] == canon, "정본 손상"
        io.open(os.path.join(W, fn), "w", newline=NL).write(
            ">%s|flank5=%d|flank3=%d|src=HPRC:%s:%s:%s%s%s%s"
            % (fn.replace(".fa", ""), len(lf), len(rf), c, strand, how, NL, new, NL))
        print("  %s 슬롯%s : 정본 %d · 앞 %d · 뒤 %d (hap%s %s %s)"
              % (smp, slot, len(canon), len(lf), len(rf), hap, strand, how))
        LOG.write(TAB.join([smp, slot, fn, c, hap, strand, str(len(lf)), str(len(rf)),
                            md5(canon), how]) + NL)
        LOG.flush()
LOG.close()
