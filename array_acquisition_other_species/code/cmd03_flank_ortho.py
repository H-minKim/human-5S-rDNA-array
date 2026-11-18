# -*- coding: utf-8 -*-
# Which human flank is the flank on the array_pos 1 side homologous to? Decided by shared k-mers.
#   The canonical fasta is [flank5 20kb][array][flank3 20kb] and is stored in the gene, that is transcription, orientation.
#   flank5 is therefore the array_pos 1 side, as confirmed by canon_start in 06_units_from_canonical, where position 1 lies at about 22 kb.
#   The test: the k-mer set of a species' flank5 is compared against the human flank5 and flank3 sets, on both strands, by Jaccard and by hit count.
# A call needs one side clearly larger than the other; otherwise it is recorded as unknown.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import glob
import io
import os
import sys

K = 21
FL = 20000
S = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1]          # the folder holding the copies of the canonical fasta files
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fa(p):
    s = []
    for L in io.open(p, encoding="utf8"):
        if not L.startswith(">"):
            s.append(L.strip().upper())
    return "".join(s)


def kset(s):
    out = set()
    for i in range(len(s) - K + 1):
        w = s[i:i + K]
        if "N" in w:
            continue
        out.add(min(w, rc(w)))        # strand-neutral
    return out


seqs = {}
for p in sorted(glob.glob(SRC + "/*.fa")):
    sp = os.path.basename(p).split("__")[0]
    s = read_fa(p)
    seqs[sp] = (kset(s[:FL]), kset(s[-FL:]))
    print("%-26s len %d" % (sp, len(s)))

h5, h3 = seqs["human_CHM13v2"]
o = io.open(S + "/FLANK_ORTHOLOGY.tsv", "w", encoding="utf8", newline="\n")
o.write("\t".join(["species", "side_tested", "n_kmer", "hit_human_flank5(array_pos1쪽)", "hit_human_flank3(반대쪽)",
                   "jaccard_f5", "jaccard_f3", "ratio_f5_over_f3", "call"]) + "\n")
for sp in sorted(seqs):
    if sp == "human_CHM13v2":
        continue
    for side, ks in (("flank5(array_pos 1 쪽)", seqs[sp][0]), ("flank3(반대쪽)", seqs[sp][1])):
        a, b = len(ks & h5), len(ks & h3)
        ja = a / max(1, len(ks | h5))
        jb = b / max(1, len(ks | h3))
        if max(a, b) < 300:
            call = "판정불가(공유 k-mer 너무 적음 · 300 미만)"
        elif a >= 3 * max(1, b):
            call = "사람 flank5 쪽(= 사람 array_pos 1 쪽)"
        elif b >= 3 * max(1, a):
            call = "사람 flank3 쪽(= 사람 반대쪽)"
        else:
            call = "모름(양쪽이 비슷)"
        o.write("\t".join([sp, side, str(len(ks)), str(a), str(b), "%.4f" % ja, "%.4f" % jb,
                           "%.2f" % (a / max(1, b)), call]) + "\n")
o.close()
print(io.open(S + "/FLANK_ORTHOLOGY.tsv", encoding="utf8").read())
