# -*- coding: utf-8 -*-
# Testing the canonical flanks: do the assembly flanks match the two ends of the canonical sequence?
#   Step: validating the reference set.
#   Produces the evidence for the reference-set quality in Supplementary Figure 1.

# -*- coding: utf-8 -*-
# Compare a consensus built from the flank reads alone against the map sequence.
#
#  The question:
#    build a consensus map from the sequences of the assigned flank reads alone,
#     and compare that read-derived sequence with the map sequence.
#
#  Why this matters:
#    the map sequence comes from assemblies (64.3% public, 33.2% our hifiasm, 2.5% spanning reads),
#    while the reads are an independent dataset. Agreement means the map holds that individual's real molecule;
#    disagreement means either the map, that is the assembly, or the read assignment is wrong.
#
#  method
#    1) Assign the reads to slots from HAP_READ_FINAL.tsv, using only the flank reads already resolved.
#    2) Cut each read into units and compute the variant set per unit.
#    3) Place each read by chain alignment to determine which map unit it starts at
#       (the chainplace.py rule: score |S intersect V| - |S - V| over K consecutive units; leave-one-out accuracy 100% at k = 8).
#    4) Per unit, pool the variants of the reads covering that position and take a majority consensus.
#    5) Compare against the map's variant list, as recall and precision.
#
#  ONT reads are about 99% accurate, so an individual variant carries noise; the result is therefore stratified by frequency
#    (measured: precision is 99.3% for high-frequency variants and 23.5% for low-frequency ones).
#  Reads whose placement is uncertain are excluded, by a lower bound on the margin.
import io
import os
import re
import sys
import gzip
import subprocess
import collections

W = os.path.expanduser("~/consflank")
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
NL = "\n"
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
K = int(os.environ.get("K", "8"))
MINMG = int(os.environ.get("MINMG", "1"))     # the minimum placement margin
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(CMP)[::-1]


def score(S, V):
    return len(S & V) - len(S - V)


MAP = collections.defaultdict(dict)
for i, L in enumerate(io.open(os.path.expanduser("~/CONSENSUS_MAP_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 15 or f[4] == "PARTIAL_PROX":
        continue
    MAP[(f[0], f[1])][int(f[3])] = set() if f[14] == "-" else set(f[14].split(","))

ASG = collections.defaultdict(list)          # smp -> [(read, slot)]
seen = set()
for i, L in enumerate(io.open(os.path.expanduser("~/HAP_READ_FINAL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 7 or f[4] == "-":
        continue
    if (f[0], f[1]) in seen:
        continue
    seen.add((f[0], f[1]))
    ASG[f[0]].append((f[1], f[4]))


def find_reads(s):
    for p in ("~/s5s_ont/hits/%s.hits.fastq.gz" % s,
              "~/s5s_new/hits/%s.hits.fastq.gz" % s,
              "~/kgp_ont/out/%s.array.fastq.gz" % s):
        q = os.path.expanduser(p)
        if os.path.exists(q) and os.path.getsize(q) > 1000:
            return q
    return None


def units_of(seq, tag):
    fa = os.path.join(W, "r_%s.fa" % tag)
    io.open(fa, "w", newline=NL).write(">r" + NL + seq + NL)
    db = os.path.join(W, "rd_%s" % tag)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                        "-evalue", "1e-5", "-max_target_seqs", "100000000",
                        "-num_threads", TH], capture_output=True, text=True)
    plus = minus = 0
    st = set()
    for L in r.stdout.splitlines():
        g = L.split("\t")
        if len(g) < 2:
            continue
        a, b = int(g[0]), int(g[1])
        if a <= b:
            plus += 1
        else:
            minus += 1
        st.add(min(a, b) - 1)
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    if minus > plus:
        return None, True
    Kp = sorted(st)
    U = []
    for i in range(len(Kp) - 1):
        if 2100 <= Kp[i + 1] - Kp[i] <= 2400:
            U.append(seq[Kp[i]:Kp[i + 1]])
    return U, False


def variants(units, tag):
    if not units:
        return []
    q = os.path.join(W, "u_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(units):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm5", "-c", "--cs", "-t", TH,
                        "--secondary=no", CANON, q], capture_output=True, text=True)
    V = {}
    for L in r.stdout.splitlines():
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
        out = set()
        for m in re.finditer(r"([:*+-])([A-Za-z0-9]+)", cs):
            op, v = m.group(1), m.group(2)
            if op == ":":
                tp += int(v)
            elif op == "*":
                out.add("%d%s" % (tp + 1, v[1].upper()))
                tp += 1
            elif op == "-":
                out.add("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                out.add("%dI%d" % (tp + 1, len(v)))
        V[i] = out
    os.remove(q)
    return [V.get(i, set()) for i in range(len(units))]


def place(chain, smp, slot):
    M = MAP.get((smp, slot))
    if not M or len(chain) < K:
        return None
    keys = sorted(M)
    best = []
    for i0 in keys:
        if not all((i0 + t) in M for t in range(K)):
            continue
        n = min(len(chain), max(keys) - i0 + 1)
        if n < K:
            continue
        s = sum(score(M[i0 + t], chain[t]) for t in range(n) if (i0 + t) in M)
        best.append((s, i0))
    if len(best) < 2:
        return None
    best.sort(reverse=True)
    return best[0][1], best[0][0] - best[1][0]


OUT = io.open(os.path.join(W, "cf_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
smps = [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]
DONE = os.path.join(W, "done")
os.makedirs(DONE, exist_ok=True)

for smp in smps:
    mk = os.path.join(DONE, smp + ".done")
    if os.path.exists(mk):
        continue
    q = find_reads(smp)
    if not q or smp not in ASG:
        continue
    want = dict(ASG[smp])
    seqs = {}
    with gzip.open(q, "rt") as fh:
        j = 0
        rid = None
        for L in fh:
            m = j & 3
            j += 1
            if m == 0:
                # The read IDs in HAP_READ_FINAL.tsv are full 36-character UUIDs.
                #   Truncating them makes the matching fail entirely.
                rid = L[1:].split()[0]
            elif m == 1 and rid in want:
                seqs[rid] = L.strip()
    # the read variants are pooled per unit position
    COV = collections.defaultdict(list)      # (slot, unit) -> list of variant sets
    nplace = nskip = 0
    for rid, slot in sorted(want.items()):
        sq = seqs.get(rid)
        if not sq:
            continue
        tag = "%s_%s_%s" % (TAG, smp, rid[:8])
        U, isrc = units_of(sq, tag)
        if isrc:
            U, _ = units_of(rc(sq), tag)
        if not U or len(U) < K:
            nskip += 1
            continue
        ch = variants(U, tag)
        p = place(ch, smp, slot)
        if not p or p[1] < MINMG:
            nskip += 1
            continue
        i0, mg = p
        nplace += 1
        for t, v in enumerate(ch):
            COV[(slot, i0 + t)].append(v)
    # a majority consensus is built per unit and compared with the map
    for (slot, u), lst in sorted(COV.items()):
        M = MAP.get((smp, slot), {}).get(u)
        if M is None:
            continue
        n = len(lst)
        c = collections.Counter()
        for v in lst:
            for x in v:
                c[x] += 1
        for th, lab in ((0.5, "maj"), (0.8, "strict")):
            cons = set(x for x, k in c.items() if k / float(n) >= th)
            OUT.write("%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d%s" %
                      (smp, slot, u, lab, n, len(M), len(cons),
                       len(M & cons), len(cons - M), NL))
    OUT.flush()
    io.open(mk, "w").write("ok")
    sys.stderr.write("  %s 배치 %d · 제외 %d · 유닛자리 %d\n" %
                     (smp, nplace, nskip, len(COV)))
    sys.stderr.flush()
OUT.close()
print("done")
