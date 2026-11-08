# -*- coding: utf-8 -*-
# [v11-fix1 copy 2026-09-16: map g8canonv10, canonical cmp4d_v11, positions = anchors given directly (HG00133 h2 119/129, HG03270 h1 40/47). Otherwise unchanged]
# Backbone polishing consensus (2026-09-12), replaces the mafft MSA version, which made chimeras in tandem arrays.
#   Method (same principle as racon, implemented directly because racon/medaka are not on the node)
#     (1) pick the read of median length among the fragments as the backbone
#     (2) map the other fragments to the backbone with minimap2 (map-ont, no secondary alignments)
#     (3) majority vote of read bases at each backbone position (substitutions, deletions and insertions all applied) -> polished sequence
#     (4) use the polished sequence as the new backbone and repeat once more (2 rounds)
#   Since no MSA is used, reads slipped by one unit do not contaminate the consensus.
#   Test: compare the anchor unit intervals (2kb before, 2kb after) with the assembly to measure base accuracy.
#   Arguments: fragment folder, output folder
import io
import os
import re
import sys
import glob
import subprocess
import collections
import statistics

H = "/home/ubuntu"
NL = chr(10)
TAB = chr(9)
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
SEGD, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
W = os.path.join(OUT, "tmp")
os.makedirs(W, exist_ok=True)
R = os.path.join(H, "s5s_tuc", "ref")
GENE = os.path.join(R, "5S_gene_query.fa")
CMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


def rdfa(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()


GLEN = len(rdfa(GENE))
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
SITE = {}
for (smp_, sl_, aL_, aR_) in [("HG00133", "2", 119, 129), ("HG03270", "1", 40, 47)]:
    v = dict((x[0], x) for x in U[(smp_, sl_)])
    gsum = sum(v[a + 1][1] - v[a][2] - 1 for a in range(aL_, aR_))      # total gap between anchors (for the record)
    SITE[(smp_, sl_, aL_)] = (v[aL_][1], v[aR_][2], gsum)


def canon(smp, sl):
    a = rdfa(os.path.join(H, "cmp4d_v11", SEQF[(smp, sl)]))
    return rc(a) if FLIP.get((smp, sl)) else a


def gene_ends(seq, tag):
    fa = os.path.join(W, "g_%s.fa" % tag)
    db = os.path.join(W, "gd_%s" % tag)
    io.open(fa, "w", newline=NL).write(">x" + NL + seq + NL)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db], capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send length pident qstart qend", "-dust", "no",
                        "-word_size", "9", "-evalue", "1e-5", "-max_target_seqs", "100000000"],
                       capture_output=True, text=True)
    hits = []
    for L in r.stdout.splitlines():
        g = L.split(TAB)
        if len(g) < 6:
            continue
        a, b, qs, qe = int(g[0]), int(g[1]), int(g[4]), int(g[5])
        if a <= b:
            s_, e_ = a - (qs - 1), b + (GLEN - qe)
        else:
            s_, e_ = b - (GLEN - qe), a + (qs - 1)
        hits.append((max(1, s_), min(len(seq), e_), int(g[2]),
                     (qs - 1) + (GLEN - qe)))   # fourth = extension length
    # F method: interval overlap plus fragment calls, in place of a 1 kb grid
    FRAG = 5
    LOx = 2100
    keep = []
    for h in sorted(hits):
        if keep and h[0] <= keep[-1][1]:
            if h[2] > keep[-1][2]:
                keep[-1] = h
        else:
            keep.append(h)
    keep.sort(key=lambda x: x[1])
    i = 0
    while i + 1 < len(keep):
        if keep[i + 1][1] - keep[i][1] < LOx:
            fa_, fb_ = keep[i][3] > FRAG, keep[i + 1][3] > FRAG
            if fa_ and not fb_:
                keep.pop(i)
                i = max(0, i - 1)
                continue
            if fb_ and not fa_:
                keep.pop(i + 1)
                i = max(0, i - 1)
                continue
        i += 1
    return [x[1] for x in keep]


def brk(ends):
    return sum(1 for i in range(len(ends) - 1) if not (2100 <= ends[i + 1] - ends[i] <= 2400))


def polish(backbone, reads, tag):
    """Attaches reads to the backbone and polishes it by per-position majority vote (substitutions, deletions, insertions)"""
    bf = os.path.join(W, "b_%s.fa" % tag)
    io.open(bf, "w", newline=NL).write(">b" + NL + backbone + NL)
    qf = os.path.join(W, "q_%s.fa" % tag)
    with io.open(qf, "w", newline=NL) as g:
        for i, s in enumerate(reads):
            g.write(">r%d%s%s%s" % (i, NL, s, NL))
    sam = subprocess.run(["minimap2", "-a", "-x", "map-ont", "-t", "8", "--secondary=no", bf, qf],
                         capture_output=True, text=True).stdout
    base = [collections.Counter() for _ in range(len(backbone))]
    ins = [collections.Counter() for _ in range(len(backbone) + 1)]
    n = 0
    for L in sam.splitlines():
        if L.startswith("@"):
            continue
        r = L.split(TAB)
        if len(r) < 10 or r[2] == "*" or (int(r[1]) & 0x900):
            continue
        n += 1
        rp = int(r[3]) - 1
        q = 0
        sq = r[9]
        for g2 in CIG.finditer(r[5]):
            k2, op = int(g2.group(1)), g2.group(2)
            if op == "S":
                q += k2
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k2):
                    if rp + t < len(base):
                        base[rp + t][sq[q + t]] += 1
                rp += k2
                q += k2
            elif op == "I":
                if rp <= len(backbone):
                    ins[rp][sq[q:q + k2]] += 1
                q += k2
            elif op in "DN":
                for t in range(k2):
                    if rp + t < len(base):
                        base[rp + t]["-"] += 1
                rp += k2
    if n == 0:
        return backbone, 0
    out = []
    for i in range(len(backbone)):
        c = ins[i]
        if c:
            s2, v = c.most_common(1)[0]
            if v * 2 > n:
                out.append(s2)
        c = base[i]
        if not c:
            out.append(backbone[i])
            continue
        b, v = c.most_common(1)[0]
        if b != "-":
            out.append(b)
    return "".join(out), n


files = collections.defaultdict(list)
for p in glob.glob(os.path.join(SEGD, "seg_*.fa")):
    core = os.path.basename(p)[4:].rsplit(".", 2)[0]
    smp, sl, ap = core.rsplit("_", 2)
    files[(smp, sl, int(ap))].append(p)

RES = io.open(os.path.join(OUT, "GAP_POLISH_RESULT.tsv"), "w", encoding="utf8", newline=NL)
RES.write(TAB.join(["sample", "slot", "after_unit", "gap_bp", "n_reads", "anchor_bp", "cons_bp",
                    "delta_bp", "acc5", "acc3", "genes_asm", "genes_cons", "breaks_asm",
                    "breaks_cons", "verdict"]) + NL)
print("%-9s %-3s %-6s %-7s %-5s %-8s %-8s %-7s %-7s %s"
      % ("개체", "h", "앞유닛", "간격", "read", "합의bp", "차이", "앞정확", "뒤정확", "판정"))


def idpair(q, ref, tag):
    qf = os.path.join(W, "iq_%s.fa" % tag)
    io.open(qf, "w", newline=NL).write(">q" + NL + q + NL)
    rf = os.path.join(W, "ir_%s.fa" % tag)
    io.open(rf, "w", newline=NL).write(">r" + NL + ref + NL)
    paf = subprocess.run(["minimap2", "-x", "asm20", "-t", "8", "--secondary=no", rf, qf],
                         capture_output=True, text=True)
    best = None
    for L in paf.stdout.splitlines():
        f2 = L.split(TAB)
        if len(f2) < 12:
            continue
        if best is None or int(f2[10]) > best[1]:
            best = (int(f2[9]), int(f2[10]))
    return (100.0 * best[0] / best[1]) if best else 0.0


for k in sorted(files):
    smp, sl, ap = k
    if k not in SITE:
        continue
    a0, a1, g = SITE[k]
    seqs = []
    for p in files[k]:
        cur = []
        for L in io.open(p, encoding="utf8"):
            if L.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(L.strip().upper())
        if cur:
            seqs.append("".join(cur))
    if len(seqs) < 2:
        continue
    tag = "%s_%s_%d" % (smp, sl, ap)
    med = statistics.median(len(s) for s in seqs)
    backbone = min(seqs, key=lambda s: abs(len(s) - med))
    cons = backbone
    for it in range(2):
        cons, used = polish(cons, [s for s in seqs if s is not backbone], tag)
    A = canon(smp, sl)
    old = A[a0 - 1:a1]
    acc5 = idpair(cons[:2000], old[:2400], "5" + tag)
    acc3 = idpair(cons[-2000:], old[-2400:], "3" + tag)
    e_o, e_c = gene_ends(old, tag + "_o"), gene_ends(cons, tag + "_c")
    b_o, b_c = brk(e_o), brk(e_c)
    if acc5 < 99 or acc3 < 99:
        v = "합의 품질 미달"
    elif b_c == 0 and b_o > 0:
        v = "조립본 오류(합의가 고침)"
    elif b_c > 0 and abs(len(cons) - len(old)) < 200:
        v = "실제 변이(합의도 같음)"
    else:
        v = "부분 개선" if b_c < b_o else "끊김 남음"
    print("%-9s %-3s %-6d %-7d %-5d %-8d %+7d %-7.2f %-7.2f 유전자 %d->%d · 끊김 %d->%d · %s"
          % (smp, sl, ap, g, len(seqs), len(cons), len(cons) - len(old), acc5, acc3,
             len(e_o), len(e_c), b_o, b_c, v))
    RES.write(TAB.join([smp, sl, str(ap), str(g), str(len(seqs)), str(len(old)), str(len(cons)),
                        str(len(cons) - len(old)), "%.2f" % acc5, "%.2f" % acc3, str(len(e_o)),
                        str(len(e_c)), str(b_o), str(b_c), v]) + NL)
    io.open(os.path.join(OUT, "pol_%s.fa" % tag), "w", newline=NL).write(
        ">pol|%s|h%s|after_unit%d|anchor%d-%d|reads%d%s%s%s" % (smp, sl, ap, a0, a1, len(seqs), NL, cons, NL))
RES.close()
