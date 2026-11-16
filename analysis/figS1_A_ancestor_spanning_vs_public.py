# -*- coding: utf-8 -*-
# Comparing array-spanning reads against the public assemblies at base-pair resolution.
#   Step: S1 A upstream.

# -*- coding: utf-8 -*-
# Do the array-spanning reads and the public assemblies agree at the level of individual positions?
#
#  The question:
#    beyond agreeing unit by unit, do the positions of the distinctive variants agree?
#
#  What is already known:
#    minimap2 overall identity: public assembly against our hifiasm assembly is 100.0000% (median),
#                          and public assembly against an array-spanning read is 99.4282% (median).
#    A spanning read therefore differs about once every 175 bp, and the nature of those differences is what matters.
#
#  method
#    For the 129 pairs matched as the same haplotype, to within one unit:
#    1) each sequence is divided into units, by the chain of 5S genes, start to start within 2,100-2,400 bp;
#    2) a variant list is built per unit against the standard 2,240 bp unit;
#    3) the variant positions are compared between units of the same index.
#
#    The unit indices are referenced to the proximal junction on both sides, so the phase matches
#    (measured: median off_prox 1,516 bp, and 98.9% of the 616 haplotypes within +-5).
#
#  Substitutions and indels are examined separately.
#    The dominant ONT error is the indel, so counting them together hides the substitution agreement.
#  The variants are stratified by frequency.
#    A variant present in most units of a haplotype is a genuine lineage marker,
#    while one present in only one or two units may be an error. The two must be separated.
import io
import os
import re
import sys
import subprocess
import collections

W = os.path.expanduser("~/psw")
os.makedirs(W, exist_ok=True)
D = os.path.expanduser("~/cmp4d")
R = os.path.expanduser("~/s5s_tuc/ref")
GENE = os.path.join(R, "5S_gene_query.fa")
CANON = os.path.join(R, "unit_canon.fa")
NL = "\n"
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(CMP)[::-1]


def units_of(p, tag):
    """Splits the sequence into units, the orientation unification included."""
    sq = "".join(x.strip() for x in
                 io.open(p, encoding="utf8").read().split(NL)[1:])
    for _ in range(2):
        fa = os.path.join(W, "g_%s.fa" % tag)
        db = os.path.join(W, "gd_%s" % tag)
        io.open(fa, "w", newline=NL).write(">x" + NL + sq + NL)
        subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
        r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                            "-outfmt", "6 sstart send", "-dust", "no",
                            "-word_size", "9", "-evalue", "1e-5",
                            "-max_target_seqs", "100000000", "-num_threads", TH],
                           capture_output=True, text=True)
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
            sq = rc(sq)
            continue
        break
    K = sorted(st)
    U = []
    for i in range(len(K) - 1):
        if 2100 <= K[i + 1] - K[i] <= 2400:
            U.append((len(U) + 1, sq[K[i]:K[i + 1]]))
    if K:
        last = min(K[-1] + 2240, len(sq))
        if last - K[-1] >= 2000:
            U.append((len(U) + 1, sq[K[-1]:last]))
    return U, (K[0] if K else -1)


def variants(U, tag):
    """Per unit variants against the canonical unit, substitutions and indels apart."""
    q = os.path.join(W, "u_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in U:
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
        sub, ind = set(), set()
        for m in re.finditer(r"([:*+-])([A-Za-z0-9]+)", cs):
            op, v = m.group(1), m.group(2)
            if op == ":":
                tp += int(v)
            elif op == "*":
                sub.add("%d%s" % (tp + 1, v[1].upper()))
                tp += 1
            elif op == "-":
                ind.add("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                ind.add("%dI%d" % (tp + 1, len(v)))
        V[i] = (sub, ind)
    os.remove(q)
    return V


OUT = io.open(os.path.join(W, "ps_%s.tsv" % TAG), "w", encoding="utf8", newline=NL)
OUT.write("\t".join(["sample", "units", "n_unit_cmp", "class", "kind",
                     "n_pub", "n_span", "n_both", "recall", "precision"]) + NL)
SITE = io.open(os.path.join(W, "pssite_%s.tsv" % TAG), "w", encoding="utf8", newline=NL)
SITE.write("sample\tsite\tkind\tclass\tpub_units\tspan_units\tboth_units\n")

rows = [L.rstrip(NL).split("\t") for L in
        io.open(os.path.expanduser("~/pubspan_pairs.tsv"), encoding="utf8") if L.strip()]
mine = [r for i, r in enumerate(rows) if i % int(sys.argv[1]) == int(sys.argv[2])]

for n, r in enumerate(mine):
    smp, nu, fpub, fspan = r[0], int(r[1]), r[2], r[3]
    pp, ps = os.path.join(D, fpub), os.path.join(D, fspan)
    if not (os.path.exists(pp) and os.path.exists(ps)):
        continue
    tag = "%s_%s_%d" % (TAG, smp, n)
    UA, ka = units_of(pp, tag + "a")
    UB, kb = units_of(ps, tag + "b")
    if not UA or not UB:
        continue
    VA = variants(UA, tag + "a")
    VB = variants(UB, tag + "b")
    m = min(len(UA), len(UB))
    if m < 5:
        continue
    # frequency stratification: in how many units of that haplotype the variant appears
    cntA = collections.Counter()
    cntB = collections.Counter()
    for i in range(1, m + 1):
        for x in VA.get(i, (set(), set()))[0] | VA.get(i, (set(), set()))[1]:
            cntA[x] += 1
        for x in VB.get(i, (set(), set()))[0] | VB.get(i, (set(), set()))[1]:
            cntB[x] += 1
    for kind, gi in (("substitution", 0), ("indel", 1)):
        for cls, lo, hi in (("high (>=50% of units)", 0.5, 1.01),
                            ("middle (10-50%)", 0.1, 0.5),
                            ("low (<10%)", 0.0, 0.1)):
            A = set()
            B = set()
            for i in range(1, m + 1):
                for x in VA.get(i, (set(), set()))[gi]:
                    if lo <= cntA[x] / float(m) < hi:
                        A.add((i, x))
                for x in VB.get(i, (set(), set()))[gi]:
                    if lo <= cntB[x] / float(m) < hi:
                        B.add((i, x))
            both = A & B
            OUT.write("\t".join([smp, str(nu), str(m), cls, kind,
                                 str(len(A)), str(len(B)), str(len(both)),
                                 "%.4f" % (len(both) / float(len(A))) if A else "-",
                                 "%.4f" % (len(both) / float(len(B))) if B else "-"]) + NL)
    # per position: in how many units the position appears in each of the two sequences
    for x in set(cntA) | set(cntB):
        a = cntA.get(x, 0)
        b = cntB.get(x, 0)
        bo = sum(1 for i in range(1, m + 1)
                 if x in (VA.get(i, (set(), set()))[0] | VA.get(i, (set(), set()))[1])
                 and x in (VB.get(i, (set(), set()))[0] | VB.get(i, (set(), set()))[1]))
        fr = max(a, b) / float(m)
        cls = ("high" if fr >= 0.5 else "middle" if fr >= 0.1 else "low")
        SITE.write("%s\t%s\t%s\t%s\t%d\t%d\t%d%s" %
                   (smp, x, "indel" if ("D" in x or "I" in x) else "substitution",
                    cls, a, b, bo, NL))
    OUT.flush()
    SITE.flush()
    sys.stderr.write("  %d/%d %s units %d\n" % (n + 1, len(mine), smp, m))
    sys.stderr.flush()
OUT.close()
SITE.close()
print("done")
