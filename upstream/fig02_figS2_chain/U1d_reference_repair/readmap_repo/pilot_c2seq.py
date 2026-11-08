# -*- coding: utf-8 -*-
# Can consensus2 be made into a sequence, a trial (the canonical files are not touched)
#
#  Procedure (the user-specified order, at the sequence level)
#    reads = reads placed in this slot in round 3 (consensus2 map)
#    ref0  = current representative sequence of the slot (spanning read array)
#    it1   : reads -> align to ref0 -> samtools consensus -> seq1
#    it2   : reads -> align to seq1 -> samtools consensus -> seq2   (= consensus2 sequence)
#  Validation
#    cut seq2 into units (same as G6: gene HSP start), call variants against the canonical unit,
#    and compare with the consensus2 variant list and with the original spanning-read map of the unit with the same index.
#    If seq2 matches consensus2 much better, "the sequence consensus reproduces consensus2".
import io
import os
import re
import sys
import gzip
import subprocess
import collections

NL = "\n"
FL = 5000
TH = "8"
SMP, SL = sys.argv[1], sys.argv[2]
W = os.path.expanduser("~/c2seq/%s_%s" % (SMP, SL))
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
C = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CS = re.compile(r"([:*+-])([A-Za-z0-9]+)")


def rc(s):
    return s.translate(C)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write("FAIL: %s\n%s\n" % (cmd, r.stderr[-400:]))
        raise SystemExit(1)
    return r.stdout


# Representative sequence
SEQF = {}
for x in io.open(os.path.expanduser("~/hap_seq_node.tsv"), encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
for i, x in enumerate(io.open(os.path.expanduser("~/CONSENSUS_SUM_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
a0 = rd(os.path.expanduser("~/cmp4d/" + SEQF[(SMP, SL)]))
if RCF.get((SMP, SL)):
    a0 = rc(a0)

# Reads placed in round 3
want = set()
for fn in os.listdir(os.path.expanduser("~/aojx7")):
    if not fn.startswith("aojv_"):
        continue
    for L in io.open(os.path.expanduser("~/aojx7/" + fn), encoding="utf8"):
        r = L.rstrip(NL).split("\t")
        if len(r) >= 3 and r[0] == SMP and r[2] == SL:
            want.add(r[1])
fq = None
for p in ("~/s5s_ont/hits/%s.hits.fastq.gz" % SMP, "~/s5s_new/hits/%s.hits.fastq.gz" % SMP):
    q = os.path.expanduser(p)
    if os.path.exists(q) and os.path.getsize(q) > 1000:
        fq = q
        break
rfq = os.path.join(W, "reads.fq")
n = 0
with io.open(rfq, "w", newline=NL) as g, gzip.open(fq, "rt") as fh:
    rec = []
    for L in fh:
        rec.append(L)
        if len(rec) == 4:
            rid = rec[0][1:].split()[0]
            if rid[:16] in want:
                g.write(rec[0] + rec[1] + "+" + NL + rec[3])
                n += 1
            rec = []
print("배치 read %d (id %d)" % (n, len(want)))


def polish(ref_arr, tag):
    ref = os.path.join(W, "%s.ref.fa" % tag)
    io.open(ref, "w", newline=NL).write(">r" + NL + rc(FD) + ref_arr + rc(FP) + NL)
    bam = os.path.join(W, "%s.bam" % tag)
    run("minimap2 -ax map-ont -t %s --secondary=no %s %s 2>/dev/null | samtools sort -@4 -o %s - && samtools index %s"
        % (TH, ref, rfq, bam, bam))
    out = os.path.join(W, "%s.cons.fa" % tag)
    reg = "r:%d-%d" % (FL + 1, FL + len(ref_arr))
    run("samtools consensus -a --show-del no --show-ins yes -r %s -f fasta -o %s %s" % (reg, out, bam))
    s = rd(out)
    dep = run("samtools depth -a -r %s %s | awk '{s+=$3; n++; if($3<3) z++} END{printf \"%%.1f %%d %%d\", s/n, z, n}'"
              % (reg, bam)).split()
    return s, float(dep[0]), int(dep[1]), int(dep[2])


FD = rd(os.path.join(R, "FLANK_distal_50000.fa"))[:FL]
FP = rd(os.path.join(R, "FLANK_proximal_50000.fa"))[-FL:]
s1, d1, z1, n1 = polish(a0, "it1")
print("it1 : 길이 %d (원래 %d) · N %d · 평균깊이 %.1f · 깊이<3 %d/%d" % (len(s1), len(a0), s1.count("N"), d1, z1, n1))
s2, d2, z2, n2 = polish(s1.replace("N", ""), "it2")
print("it2 : 길이 %d · N %d · 평균깊이 %.1f · 깊이<3 %d/%d" % (len(s2), s2.count("N"), d2, z2, n2))
io.open(os.path.join(W, "cons2.fa"), "w", newline=NL).write(">%s_%s_cons2%s%s%s" % (SMP, SL, NL, s2, NL))


def units(seq, tag):
    fa = os.path.join(W, "%s.u.fa" % tag)
    io.open(fa, "w", newline=NL).write(">s" + NL + seq + NL)
    run("makeblastdb -in %s -dbtype nucl -out %s.db >/dev/null" % (fa, fa))
    out = run("blastn -task blastn -query %s -db %s.db -outfmt '6 sstart send' -dust no -word_size 9 -evalue 1e-5 -max_target_seqs 100000000"
              % (GENE, fa))
    st = sorted(set(min(int(a), int(b)) - 1 for a, b in (x.split("\t") for x in out.splitlines() if x)))
    return [seq[st[i]:st[i + 1]] for i in range(len(st) - 1) if 2100 <= st[i + 1] - st[i] <= 2400]


def varof(us, tag):
    q = os.path.join(W, "%s.v.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(us):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    out = run("minimap2 -x asm5 -c --cs -t %s --secondary=no %s %s" % (TH, CANON, q))
    V = {}
    for L in out.splitlines():
        f = L.split("\t")
        if len(f) < 12:
            continue
        i = int(f[0][1:])
        if i in V:
            continue
        cs = ""
        for x in f[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        tp = int(f[7])
        o = set()
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == ":":
                tp += int(v)
            elif op == "*":
                o.add("%d%s" % (tp + 1, v[1].upper()))
                tp += 1
            elif op == "-":
                o.add("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                o.add("%dI%d" % (tp + 1, len(v)))
        V[i] = o
    return [V.get(i) for i in range(len(us))]


V2 = varof(units(s2, "c2"), "c2")
# Comparison targets: consensus2 (MAP_C2) and the original map (G6)
def loadmap(fn):
    M = {}
    f = io.open(os.path.expanduser(fn), encoding="utf8")
    h = f.readline().rstrip(NL).split("\t")
    ix = dict((x, i) for i, x in enumerate(h))
    for L in f:
        r = L.rstrip(NL).split("\t")
        if len(r) < 15 or r[0] != SMP or r[1] != SL or r[ix["seg_type"]] != "UNIT":
            continue
        M[int(r[ix["array_pos"]])] = set() if r[ix["varlist"]] == "-" else set(r[ix["varlist"]].split(","))
    return [M[k] for k in sorted(M)]


C2 = loadmap("~/MAP_C2.tsv")
G6 = loadmap("~/CONSENSUS_MAP_G6.tsv")
print("유닛 수 : 서열합의2 %d · 합의2 지도 %d · 원래 지도 %d" % (len(V2), len(C2), len(G6)))


def best_off(A, B):
    best = None
    for off in range(-3, 4):
        tp = fn = fp = 0
        for i in range(len(A)):
            j = i + off
            if 0 <= j < len(B) and A[i] is not None:
                tp += len(A[i] & B[j]); fn += len(B[j] - A[i]); fp += len(A[i] - B[j])
        t = tp + fn + fp
        if t and (best is None or tp / float(t) > best[0]):
            best = (tp / float(t), off, tp, fn, fp)
    return best


for nm, B in (("합의2 변이목록", C2), ("원래 관통read 지도", G6)):
    b = best_off(V2, B)
    print("  서열합의2 vs %-16s : 일치 %.3f%% (오프셋 %+d · TP %d FN %d FP %d)" % (nm, 100 * b[0], b[1], b[2], b[3], b[4]))
