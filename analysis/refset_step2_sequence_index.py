# -*- coding: utf-8 -*-
# The unit count per file, written to seqidx_all.tsv.
#   Step 2: the sequence index.
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# The sequence index: for every .fa in cmp4d, the individual, the source, the bp and the unit count are measured exhaustively.
#
#  The unit count is in the unit frame, that is the number of normal 2,100-2,400 bp gaps plus 1.
#    This is the same frame as hap_final2.tsv. Dividing bp by 2240 is never used.
#  Sorting is done only with Python's sorted(); awk's string sort scrambles the coordinates.
import io
import os
import re
import sys
import subprocess

D = os.path.expanduser("~/cmp4d")
W = os.path.expanduser("~/seqidxw")
os.makedirs(W, exist_ok=True)
GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
NL = "\n"

PUB_D = {"asmunits3", "hgsvcunits"}
SPAN_D = {"tier1ref"}
HIFI_D = {"arr", "resub", "resub2", "hicheck", "refix", "fill60", "sweep31", "nopurge"}


def base_sample(b, kind):
    if kind == "span":
        s = b[3:] if b.startswith("t1_") else b
        return s.split("_p")[0]
    s = b[6:] if b.startswith("hiarr_") else (b[4:] if b.startswith("arr_") else b)
    s = re.sub(r"_\d+$", "", s)                 # strip the trailing _N first
    for pat in (r"_(sw|f)[2-6]$", r"_l0$", r"_sub250$", r"_noSPAN$"):
        s = re.sub(pat, "", s)
    return re.sub(r"_\d+$", "", s)


def units_of(sq, tag):
    fa = os.path.join(W, "u_%s.fa" % tag)
    io.open(fa, "w", newline=NL).write(">x" + NL + sq + NL)
    db = os.path.join(W, "ud_%s" % tag)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                        "-evalue", "1e-5", "-max_target_seqs", "100000000",
                        "-num_threads", TH], capture_output=True, text=True)
    st = set()
    for L in r.stdout.splitlines():
        g = L.split("\t")
        if len(g) >= 2:
            st.add(min(int(g[0]), int(g[1])) - 1)
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    K = sorted(st)                                  # python sorted
    gaps = [K[i + 1] - K[i] for i in range(len(K) - 1)]
    good = sum(1 for x in gaps if 2100 <= x <= 2400)
    frac = 100.0 * good / len(gaps) if gaps else 0.0
    return good + 1, frac, len(K)


files = sorted(x for x in os.listdir(D) if x.endswith(".fa"))
mine = [f for i, f in enumerate(files) if i % int(sys.argv[1]) == int(sys.argv[2])]
OUT = io.open(os.path.join(W, "seqidx_%s.tsv" % TAG), "w", encoding="utf8", newline=NL)
OUT.write("file\tsample\tkind\tdir\tbp\tunits\tchain_pct\tn_hsp\n")
for n, f in enumerate(mine):
    p = f.split("__")
    if len(p) < 3:
        continue
    d, b = p[0], p[2][:-3]
    kind = ("pub" if d in PUB_D else "span" if d in SPAN_D
            else "hifi" if d in HIFI_D else None)
    if not kind:
        continue
    L = io.open(os.path.join(D, f), encoding="utf8").read().split(NL)
    sq = "".join(x.strip() for x in L[1:])
    if len(sq) < 4000:
        continue
    u, fr, nh = units_of(sq, "%s_%d" % (TAG, n))
    OUT.write("%s\t%s\t%s\t%s\t%d\t%d\t%.1f\t%d%s"
              % (f, base_sample(b, kind), kind, d, len(sq), u, fr, nh, NL))
    OUT.flush()
    if n % 25 == 0:
        sys.stderr.write("  %d / %d\n" % (n, len(mine)))
        sys.stderr.flush()
OUT.close()
print("done %d" % len(mine))
