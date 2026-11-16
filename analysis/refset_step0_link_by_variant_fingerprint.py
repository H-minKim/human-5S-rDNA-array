# -*- coding: utf-8 -*-
# Link the remaining 72 individuals by a chain of variant fingerprints.
#   Step 1: fix the haplotypes.

# -*- coding: utf-8 -*-
# Linking the proximal and distal haplotypes of the remaining 72 individuals, by a chain of variant fingerprints.
#
#  The problem is a narrow one.
#    For all 72, the haplotypes are already separated within each side (hap_reads_v5.tsv:
#    all four groups prox1, prox2, dist1 and dist2 exist in 72 of 72, with a median read-count balance of 0.80).
#    What is missing is only the link between the proximal and the distal groups:
#    per individual, whether it is prox1-dist1 or prox1-dist2, a choice between two options.
#    It is not a problem of re-placing the reads.
#
#  The method is the one already validated here, in chainplace.py and loo2.py:
#    the score is |S intersect V| - |S - V|, that is matches minus misses,
#    where S is the map's variant set for that unit and V the read's,
#    over a chain of K consecutive units. Leave-one-out accuracy measured at 49% for k = 1, 70% for k = 2, 94% for k = 4 and 100% for k = 8.
#    minimap2 is not used: measured, it is not a tool for assigning haplotypes
#      (accuracy measured against spanning reads is 42.1% at MAPQ 1-19, and 44% of proximal flank reads fall below MAPQ 20).
#
#  The map's variant list serves as the fingerprint, as in chainplace.py.
#    The map is a consensus over many reads, so it plays the part of the shared variants.
#
#  A shared fingerprint throws information away, so this proceeds in two stages:
#    stage 1 decides from the group's shared variants, those carried by a majority of its reads;
#    stage 2, where stage 1 leaves no margin, sums the scores of the individual reads
#             and so discards nothing.
import io
import os
import sys
import gzip
import subprocess
import collections

W = os.path.expanduser("~/linkhap")
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
NL = "\n"
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
K = int(os.environ.get("K", "8"))
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(CMP)[::-1]


def score(S, V):
    return len(S & V) - len(S - V)


# ---- the map: (individual, slot) -> {unit index: variant set} ----
MAP = collections.defaultdict(dict)
for i, L in enumerate(io.open(os.path.expanduser("~/CONSENSUS_MAP_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 15 or f[4] == "PARTIAL_PROX":
        continue
    vs = set()
    if f[14] != "-":
        for v in f[14].split(","):
            vs.add(v)
    MAP[(f[0], f[1])][int(f[3])] = vs

# ---- the side label of each read ----
SIDE = collections.defaultdict(dict)
for i, L in enumerate(io.open(os.path.expanduser("~/hap_reads_v5.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 4:
        continue
    p = f[2].split("|")
    rid = p[1] if len(p) > 1 else f[2]
    SIDE[f[0]][(rid, f[1])] = f[3]


def find_reads(s):
    for p in ("~/s5s_ont/hits/%s.hits.fastq.gz" % s,
              "~/s5s_new/hits/%s.hits.fastq.gz" % s,
              "~/kgp_ont/out/%s.array.fastq.gz" % s):
        q = os.path.expanduser(p)
        if os.path.exists(q) and os.path.getsize(q) > 1000:
            return q
    return None


def unit_chain(seq, tag):
    """cut the read into units and give, for each unit, its variant set against the standard unit."""
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
    import re
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
    """place the chain on the map of that slot and return (best score, margin)."""
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
    return best[0][0], best[0][0] - best[1][0]


OUT = io.open(os.path.join(W, "link_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
smps = [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]
DONE = os.path.join(W, "done")
os.makedirs(DONE, exist_ok=True)

for smp in smps:
    mk = os.path.join(DONE, smp + ".done")
    if os.path.exists(mk):
        continue
    q = find_reads(smp)
    slots = sorted(s for (a, s) in MAP if a == smp)
    if not q or len(slots) < 2:
        sys.stderr.write("  건너뜀 %s (read %s · 슬롯 %d)\n" % (smp, bool(q), len(slots)))
        continue
    want = set(r for (r, sd) in SIDE.get(smp, {}))
    seqs = {}
    with gzip.open(q, "rt") as fh:
        j = 0
        rid = None
        for L in fh:
            m = j & 3
            j += 1
            if m == 0:
                rid = L[1:].split()[0]
            elif m == 1 and rid in want:
                seqs[rid] = L.strip()
    # the score summed per group
    SC = collections.defaultdict(lambda: collections.Counter())
    NR = collections.Counter()
    for (rid, sd), h in sorted(SIDE.get(smp, {}).items()):
        sq = seqs.get(rid)
        if not sq:
            continue
        tag = "%s_%s_%s" % (TAG, smp, rid[:8])
        U, isrc = unit_chain(sq, tag)
        if isrc:
            U, _ = unit_chain(rc(sq), tag)
        if not U or len(U) < K:
            continue
        ch = variants(U, tag)
        NR[(sd, h)] += 1
        for slot in slots:
            p = place(ch, smp, slot)
            if p:
                SC[(sd, h)][slot] += p[0]
    # group -> slot
    call = {}
    for k in SC:
        c = SC[k]
        if len(c) < 2:
            continue
        a = c.most_common()
        call[k] = (a[0][0], a[0][1] - a[1][1], NR[k])
    for (sd, h), (slot, mg, nr) in sorted(call.items()):
        OUT.write("%s\t%s\t%s\t%s\t%d\t%d\n" % (smp, sd, h, slot, mg, nr))
    OUT.flush()
    io.open(mk, "w").write("ok")
    sys.stderr.write("  %s 군 %d 판정 (read %d)\n" % (smp, len(call), sum(NR.values())))
    sys.stderr.flush()
OUT.close()
print("done")
