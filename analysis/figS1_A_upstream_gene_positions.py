# -*- coding: utf-8 -*-
# makeblastdb + blastn → genepos.tsv
#   Step: S1 A upstream.
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# The gene start positions are recorded as they are, which separates a real difference in unit count from a boundary that was never drawn.
#   Comparing whole fasta lengths mixes in how much flank each sequence happens to carry.
#        To see the array alone, the distance from the first gene to the last must be used.
#   Reported: the file, the number of gene starts, K[0], K[-1], the interior distance K[-1] - K[0], the total length,
#             the leading flank K[0], the trailing flank len - K[-1],
#             and the gap distribution: normal (2100-2400), large (over 2400, a missed boundary) and small (under 2100).
#             The large gaps are also expressed in units: 4480 is about 2 units, so one boundary was missed.
import io, os, sys, glob, subprocess
W = os.path.expanduser("~/psw3"); os.makedirs(W, exist_ok=True)
D = os.path.expanduser("~/cmp4d")
R = os.path.expanduser("~/s5s_tuc/ref")
GENE = os.path.join(R, "5S_gene_query.fa")
NL = "\n"; TH = os.environ.get("TH", "3")
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(CMP)[::-1]


def starts(p, tag):
    sq = "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])
    for _ in range(2):
        fa = os.path.join(W, "g_%s.fa" % tag); db = os.path.join(W, "gd_%s" % tag)
        io.open(fa, "w", newline=NL).write(">x" + NL + sq + NL)
        subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
        r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                            "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                            "-evalue", "1e-5", "-max_target_seqs", "100000000",
                            "-num_threads", TH], capture_output=True, text=True)
        plus = minus = 0; st = set()
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
            sq = rc(sq); continue
        break
    return sorted(st), len(sq)


FL = [os.path.join(D, x.strip()) for x in
      io.open(os.path.expanduser("~/genepos_list.txt")) if x.strip()]
out = io.open(os.path.expanduser("~/genepos.tsv"), "w", encoding="utf8", newline=NL)
out.write("file\tn_gene\tk0\tklast\tinner_bp\ttotal_bp\tflank5\tflank3"
          "\tgap_ok\tgap_big\tgap_small\tbig_gaps\n")
for i, p in enumerate(FL):
    b = os.path.basename(p)
    try:
        K, L = starts(p, "gp%d" % i)
    except Exception:
        out.write("%s\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t0\t0\t0\t-\n" % b); out.flush(); continue
    if not K:
        out.write("%s\t0\t-1\t-1\t-1\t%d\t-1\t-1\t0\t0\t0\t-\n" % (b, L)); out.flush(); continue
    g = [K[j + 1] - K[j] for j in range(len(K) - 1)]
    ok = sum(1 for x in g if 2100 <= x <= 2400)
    big = [x for x in g if x > 2400]
    small = sum(1 for x in g if x < 2100)
    out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
              % (b, len(K), K[0], K[-1], K[-1] - K[0], L, K[0], L - K[-1],
                 ok, len(big), small,
                 ",".join("%.2f" % (x / 2240.0) for x in sorted(big, reverse=True)[:8]) or "-"))
    out.flush()
    sys.stderr.write("%d/%d %s\n" % (i + 1, len(FL), b)); sys.stderr.flush()
print("done")
