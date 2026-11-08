# -*- coding: utf-8 -*-
# 5S methylation extraction
#   input: MM/ML of placed BAMs (junction ~/placed_bam, ARRAY_ONLY ~/placed_bam_g8), raw read as is
#          unit boundaries: unit_start/unit_end of the map CONSENSUS_MAP (canonV8)
#   coordinate rule (verified with sequence in the pilot, CpG hit rate 0.986)
#     - MM is relative to the original read orientation. For reverse-aligned records the SAM SEQ is flipped, so
#       original coordinate op -> SEQ coordinate sp = len(SEQ)-1-op, then moved to reference coordinates with CIGAR.
#     - a C of a reverse read lands on the G position of the canonical forward -> subtract 1 to match the C of the CpG.
# (without this correction the hit rate drops to 0.65)
#     - junction BAMs have a [flank5000][canonical][flank5000] frame, so reference coordinate -5000. ARRAY_ONLY as is.
#   outputs (gzip)
#     METH_READUNIT_<node>.tsv.gz : individual, source, setting, slot, unit number, read, strand, CpG count, methylated count, beta, calls
#         calls = list of "within-unit coordinate:value" (value 0-255, ML as is) -> each region of each read can be viewed directly
#     METH_SITE_<node>.tsv.gz     : aggregated by individual, source, slot, unit number and within-unit coordinate (n, n_meth, beta)
#   threshold: methylation call ML >= 128 (probability 0.5)
#   args: output folder [individual list file]
import io
import os
import re
import sys
import gzip
import bisect
import subprocess
import collections

H = "/home/ubuntu"
NL = chr(10)
TAB = chr(9)
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
MLCUT = 128
OUT = sys.argv[1]
LIST = sys.argv[2] if len(sys.argv) > 2 else None
os.makedirs(OUT, exist_ok=True)
NODE = subprocess.run(["hostname"], capture_output=True, text=True).stdout.strip()
MAPF = os.path.join(H, "g8jxmap_v10", "CONSENSUS_MAP_G8_full.tsv")


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


# ── unit boundaries (canonical coordinates, 1-based inclusive)
U = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] != "UNIT":
        continue
    U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]), int(r[ix["array_pos"]])))
for k in U:
    U[k].sort()
STARTS = dict((k, [x[0] for x in v]) for k, v in U.items())

SMPS = None
if LIST:
    SMPS = set(x.strip() for x in io.open(LIST, encoding="utf8") if x.strip())


def bams():
    for d, off, kind in ((os.path.join(H, "placed_bam"), 5000, "junction"),
                         (os.path.join(H, "NONE_placed_bam_g8"), 0, "array_only")):
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".bam"):
                continue
            p = fn.split(".")
            smp, cfg, src = p[0], p[1], p[2]
            if SMPS and smp not in SMPS:
                continue
            yield os.path.join(d, fn), smp, cfg, src, off, kind


def calls_of(seq, flag, mm, ml):
    """C+m probabilities (SEQ coordinate, value). MM is relative to the original read orientation."""
    orig = rc(seq) if (flag & 16) else seq
    vals = [int(x) for x in ml.split(",") if x != ""]
    out = []
    base = 0
    idx = None
    for blk in [p for p in mm.split(";") if p]:
        f2 = blk.split(",")
        deltas = [int(x) for x in f2[1:]]
        if f2[0].startswith("C+m"):
            if idx is None:
                idx = [i for i, c in enumerate(orig) if c == "C"]
            j = -1
            for t, d in enumerate(deltas):
                j += d + 1
                if j >= len(idx):
                    break
                op = idx[j]
                sp = (len(seq) - 1 - op) if (flag & 16) else op
                out.append((sp, vals[base + t] if base + t < len(vals) else 0))
        base += len(deltas)
    return out


RU = gzip.open(os.path.join(OUT, "METH_READUNIT_%s.tsv.gz" % NODE), "wt", encoding="utf8", newline=NL)
RU.write(TAB.join(["sample", "data_source", "cfg", "kind", "slot", "array_pos", "read_id",
                   "strand", "n_cpg", "n_meth", "beta", "calls"]) + NL)
SITE = collections.Counter()
SITEM = collections.Counter()
T = collections.Counter()

for b, smp, cfg, src, off, kind in bams():
    pr = subprocess.Popen(["samtools", "view", b], stdout=subprocess.PIPE, text=True)
    for L in pr.stdout:
        r = L.rstrip(NL).split(TAB)
        if len(r) < 11:
            continue
        flag = int(r[1])
        if flag & 0x900 or r[2] == "*":
            continue
        sn, pos, cig, seq = r[2], int(r[3]), r[5], r[9]
        sl = sn.split("_h")[-1]
        key = (smp, sl)
        if key not in U:
            T["슬롯 지도 없음"] += 1
            continue
        mm = ml = None
        for x in r[11:]:
            if x.startswith("MM:Z:"):
                mm = x[5:]
            elif x.startswith("ML:B:C,"):
                ml = x[7:]
        if not mm or not ml:
            T["MM/ML 없음"] += 1
            continue
        m = {}
        rp = pos - 1
        q = 0
        for g in CIG.finditer(cig):
            k, op = int(g.group(1)), g.group(2)
            if op == "S":
                q += k
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k):
                    m[q + t] = rp + t
                rp += k
                q += k
            elif op == "I":
                q += k
            elif op in "DN":
                rp += k
        st = "-" if (flag & 16) else "+"
        per = collections.defaultdict(list)
        for sp, v in calls_of(seq, flag, mm, ml):
            if sp not in m:
                continue
            cp = m[sp] - off                      # canonical 0-based
            if flag & 16:
                cp -= 1                           # reverse C -> C of the CpG
            p1 = cp + 1                           # canonical 1-based
            if p1 < 1:
                continue
            S = STARTS[key]
            i = bisect.bisect_right(S, p1) - 1
            if i < 0:
                continue
            s0, e0, ap = U[key][i]
            if p1 > e0:
                continue
            per[ap].append((p1 - s0 + 1, v))
        T["read"] += 1
        for ap, lst in per.items():
            lst.sort()
            nm = sum(1 for _, v in lst if v >= MLCUT)
            RU.write(TAB.join([smp, src, cfg, kind, sl, str(ap), r[0], st, str(len(lst)), str(nm),
                               "%.4f" % (nm / float(len(lst))),
                               ",".join("%d:%d" % (a, v) for a, v in lst)]) + NL)
            T["read-유닛"] += 1
            for a, v in lst:
                SITE[(smp, src, sl, ap, a)] += 1
                if v >= MLCUT:
                    SITEM[(smp, src, sl, ap, a)] += 1
    pr.stdout.close()
    pr.wait()
RU.close()

SI = gzip.open(os.path.join(OUT, "METH_SITE_%s.tsv.gz" % NODE), "wt", encoding="utf8", newline=NL)
SI.write(TAB.join(["sample", "data_source", "slot", "array_pos", "unit_pos", "n", "n_meth", "beta"]) + NL)
for k in sorted(SITE):
    n = SITE[k]
    nm = SITEM.get(k, 0)
    SI.write(TAB.join([k[0], k[1], k[2], str(k[3]), str(k[4]), str(n), str(nm),
                       "%.4f" % (nm / float(n))]) + NL)
SI.close()
for k, v in sorted(T.items()):
    print("  %-16s %d" % (k, v))
print("  자리 %d" % len(SITE))
