# -*- coding: utf-8 -*-
# Build tier1 references: use the array segment of spanning reads as per-haplotype references
#
#  haplotype grouping: group by (proximal hap, distal hap) pairs.
#    proximal and distal numbers are independent, so using only one side mixes half of them
#    (measured: only 46.9% of reads have proximal=distal). Spanning reads carry both sides, so
#    the pair fixes the haplotype. Each individual should have at most 2 pair kinds, and 97.2% actually do.
#
#  If a haplotype has several spanning reads, the one with the longest array bp is used.
#    Max bp, not max unit count: in measurement, pairs off by 4 or more units had almost the same bp
#    (ratio 0.10~0.28), i.e. the sequence was fine and it was a counting error where blastn missed gene HSPs.
#    Since the sequence is used as a reference, the criterion must be the length actually covered.
#
#  array boundaries are set by two 50 kb probes outside the CHM13 junctions (same as other outputs).
import io
import os
import sys
import gzip
import subprocess
import collections

W = os.path.expanduser("~/tier1ref")
os.makedirs(W, exist_ok=True)
PROBE = os.path.expanduser("~/asmunits3/probe.fa")
GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
NL = "\n"

# haplotype table (side must be part of the key)
H = collections.defaultdict(dict)
for hp in ("~/v5k/hap_reads_v5.tsv", "~/giabv5/hap_reads.tsv",
           "~/hg008v5/hap_reads.tsv"):
    p = os.path.expanduser(hp)
    if not os.path.exists(p):
        continue
    for i, L in enumerate(io.open(p, encoding="utf8")):
        if i == 0:
            continue
        f = L.rstrip(NL).split("\t")
        if len(f) >= 4:
            q = f[2].split("|")
            if len(q) >= 2:
                H[(f[0], q[1])][f[1]] = int(f[3])

SPAN = collections.defaultdict(list)
for L in io.open(os.path.expanduser("~/v5u/units_all.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) >= 7 and f[2] == "SPAN":
        SPAN[f[0]].append(f[1])


def find_reads(smp):
    for p in ("~/s5s_ont/hits/%s.hits.fastq.gz",
              "~/kgp_ont/out/%s.array.fastq.gz",
              "~/giab_hits/%s.hits.fastq.gz",
              "~/cancer5s/hits/%s.hits.fastq.gz"):
        q = os.path.expanduser(p % smp)
        if os.path.exists(q):
            return q
    return None


OUT = io.open(os.path.join(W, "tier1_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
for smp in sys.argv[1:]:
    if smp not in SPAN:
        continue
    q = find_reads(smp)
    if not q:
        continue
    want = set(SPAN[smp])
    seqs = {}
    with gzip.open(q, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline()
            fh.readline()
            rid = h[1:].rstrip(NL).replace(" ", "\t").split("\t")[0]
            if rid in want:
                seqs[rid] = s
    if not seqs:
        continue
    fa = os.path.join(W, "s_%s.fa" % smp)
    with io.open(fa, "w", newline=NL) as g:
        for r, s in seqs.items():
            g.write(u">%s\n%s\n" % (r, s))
    # probes use -N 20 -p 0.1. With --secondary=no a probe attaches to only one read.
    paf = subprocess.run(["minimap2", "-x", "map-ont", "-t", TH, "-N", "20", "-p", "0.1",
                          fa, PROBE], capture_output=True, text=True)
    FL = collections.defaultdict(dict)
    for L in paf.stdout.splitlines():
        f = L.split("\t")
        if len(f) < 12:
            continue
        pq, t, ts, te, nm = f[0], f[5], int(f[7]), int(f[8]), int(f[9])
        if nm < 3000:
            continue
        cur = FL[t].get(pq)
        if cur is None or nm > cur[2]:
            FL[t][pq] = (ts, te, nm)
    os.remove(fa)

    cand = collections.defaultdict(list)   # (ph,dh) -> [(bp, read, s, e)]
    for r in seqs:
        d = FL.get(r, {})
        if "PROX50k" not in d or "DIST50k" not in d:
            continue
        hh = H.get((smp, r), {})
        if "prox" not in hh or "dist" not in hh:
            continue
        pm = (d["PROX50k"][0] + d["PROX50k"][1]) // 2
        dm = (d["DIST50k"][0] + d["DIST50k"][1]) // 2
        a = sorted([d["PROX50k"][0], d["PROX50k"][1],
                    d["DIST50k"][0], d["DIST50k"][1]])
        s0, e0 = a[1], a[2]
        if e0 - s0 < 4000:
            continue
        side = "PROX_FIRST" if pm < dm else "DIST_FIRST"
        cand[(hh["prox"], hh["dist"])].append((e0 - s0, r, s0, e0, side))
    for k in sorted(cand):
        cand[k].sort(reverse=True)
        bp, r, s0, e0, side = cand[k][0]
        seq = seqs[r][s0:e0]
        # references are unified in proximal->distal orientation (to match the coordinate system of other outputs)
        if side == "DIST_FIRST":
            seq = seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]
        nm = "%s_p%d_d%d" % (smp, k[0], k[1])
        io.open(os.path.join(W, "t1_%s.fa" % nm), "w", newline=NL).write(
            ">%s|tier1|%s|%d%s%s%s" % (nm, r, bp, NL, seq, NL))
        OUT.write("%s\t%d\t%d\t%s\t%d\t%d\t%s%s"
                  % (smp, k[0], k[1], r, bp, len(cand[k]), side, NL))
        OUT.flush()
        sys.stderr.write("  %s (근위h%d,원위h%d) read %s · array %d bp · 후보 %d · %s\n"
                         % (smp, k[0], k[1], r[:12], bp, len(cand[k]), side))
        sys.stderr.flush()
OUT.close()
print("done")
