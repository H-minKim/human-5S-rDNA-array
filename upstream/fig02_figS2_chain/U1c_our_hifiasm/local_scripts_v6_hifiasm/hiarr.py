# -*- coding: utf-8 -*-
# extract the 5S array from our hifiasm assemblies and count units
#   exactly the same method as asm_arr_all.py / hicheck.py (if it differs, the comparison does not hold)
#     (1) align two 50 kb probes outside the CHM13 junctions to the contig (both sides must attach)
#     (2) take the span between the two probes as the array
#     (3) count as units the 2,100-2,400 bp intervals between 5S gene blastn HSP starts
#   results are streamed to standard output one line at a time (the calling script appends them to the summary table).
import io
import os
import sys
import subprocess
import collections

GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
PP = os.path.expanduser("~/asmunits3/probe.fa")
TH = "4"
NL = "\n"
SMP, D, ARR = sys.argv[1], sys.argv[2], sys.argv[3]


def units_of(seq, tmp):
    fa = os.path.join(tmp, "u.fa")
    io.open(fa, "w", newline=NL).write(">arr" + NL + seq + NL)
    db = os.path.join(tmp, "ud")
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                        "-evalue", "1e-5", "-max_target_seqs", "100000000",
                        "-num_threads", TH], capture_output=True, text=True, check=True)
    st = set()
    for L in r.stdout.splitlines():
        f = L.split("\t")
        if len(f) >= 2:
            st.add(min(int(f[0]), int(f[1])) - 1)
    K = sorted(st)
    return sum(1 for i in range(len(K) - 1) if 2100 <= K[i + 1] - K[i] <= 2400)


# all p_ctg into one fasta (hap1/hap2/primary all examined, duplicates cleaned up later)
fa = os.path.join(D, "all_ctg.fa")
seqs = {}
for g in sorted(os.listdir(D)):
    if not g.endswith(".p_ctg.gfa") or "noseq" in g:
        continue
    src = "hap1" if "hap1" in g else ("hap2" if "hap2" in g else "pri")
    for L in io.open(os.path.join(D, g), encoding="utf8"):
        if L[0] != "S":
            continue
        f = L.rstrip(NL).split("\t")
        if len(f) < 3 or f[2] == "*":
            continue
        seqs["%s|%s" % (src, f[1])] = f[2]
if not seqs:
    sys.stdout.write("%s\tNO_CTG\t-\t-\t-\t-%s" % (SMP, NL))
    sys.exit()
with io.open(fa, "w", newline=NL) as g:
    for k, v in seqs.items():
        g.write(u">%s\n%s\n" % (k, v))

paf = subprocess.run(["minimap2", "-x", "asm20", "-t", TH, "-N", "20", "-p", "0.1",
                      fa, PP], capture_output=True, text=True)
hit = collections.defaultdict(dict)
for L in paf.stdout.splitlines():
    f = L.split("\t")
    if len(f) < 12:
        continue
    pq, t, ts, te, nm = f[0], f[5], int(f[7]), int(f[8]), int(f[9])
    if nm < 5000:
        continue
    cur = hit[t].get(pq)
    if cur is None or nm > cur[2]:
        hit[t][pq] = (ts, te, nm)

got = []
for t, d in hit.items():
    if "PROX50k" not in d or "DIST50k" not in d:
        continue
    a = sorted([d["PROX50k"][0], d["PROX50k"][1], d["DIST50k"][0], d["DIST50k"][1]])
    seg = seqs[t][a[1]:a[2]]
    if len(seg) < 4000:
        continue
    got.append((t, a[2] - a[1], seg, d["PROX50k"][2] + d["DIST50k"][2]))

if not got:
    # also record contigs where only one probe attached (needed to know how much was lost)
    half = sum(1 for t, d in hit.items() if len(d) == 1)
    sys.stdout.write("%s\tNO_BOTH_PROBE\t%d\t-\t-\t-%s" % (SMP, half, NL))
    sys.exit()

# collapse duplicate contigs of the same length into one (hap1/hap2/pri emit the same one several times)
got.sort(key=lambda z: (-z[1], -z[3]))
keep = []
for t, L, seg, nm in got:
    if any(abs(L - k[1]) <= 200 for k in keep):
        continue
    keep.append((t, L, seg, nm))

for i, (t, L, seg, nm) in enumerate(keep[:4]):
    u = units_of(seg, D)
    io.open(os.path.join(ARR, "hiarr_%s_%d.fa" % (SMP, i + 1)), "w",
            newline=NL).write(">%s|%s|%d|%d%s%s%s" % (SMP, t, L, u, NL, seg, NL))
    sys.stdout.write("%s\tOK\t%d\t%s\t%d\t%d%s" % (SMP, len(keep), t, L, u, NL))
