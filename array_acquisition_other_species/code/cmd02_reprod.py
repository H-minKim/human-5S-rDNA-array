# -*- coding: utf-8 -*-
# Does an array-level value reproduce between individuals? Zebrafish ONT, 6 individuals, the 30 chr4 arrays.
#   Input: 13_array_level_diff/ARRAY_LEVEL_TABLE.tsv, the output of cmd02_stats.py.
#   Outputs: ZF_CHR4_ARRAY_REPRO.tsv, the whole-array methylation per array and individual with the read count, plus the mean over the 6 individuals;
#          ZF_CHR4_ARRAY_REPRO_SPEARMAN.tsv, a Spearman correlation per pair of individuals, over the arrays with at least 20 reads in both.
#   Each individual's value is measured within that individual alone; the individuals are never pooled, and only whether the ranking agrees is examined.
# MIT License · Copyright (c) 2026 Hungmin Kim
import csv
import io
import itertools
import numpy as np
from scipy.stats import spearmanr

D = "$RDNA_OTHERSP_DIR/data_v2/13_array_level_diff/"
T = "\t"
S = ["Control%d" % i for i in range(1, 7)]
d = {}
for r in csv.DictReader(io.open(D + "ARRAY_LEVEL_TABLE.tsv", encoding="utf8"), delimiter=T):
    if r["species"] == "zebrafish" and r["chrom"] == "chr4":
        d.setdefault(r["locus"], {})[r["sample"]] = (float(r["whole"]), int(r["n_read"]))
L = sorted(d, key=lambda l: np.mean([d[l][s][0] for s in S if s in d[l]]))
o = io.open(D + "ZF_CHR4_ARRAY_REPRO.tsv", "w", encoding="utf8", newline="\n")
o.write(T.join(["locus"] + ["%s_whole" % s for s in S] + ["%s_nread" % s for s in S] + ["mean_whole"]) + "\n")
for l in L:
    o.write(T.join([l] + ["%.3f" % d[l][s][0] if s in d[l] else "-" for s in S] +
                   [str(d[l][s][1]) if s in d[l] else "0" for s in S] +
                   ["%.3f" % np.mean([d[l][s][0] for s in S if s in d[l]])]) + "\n")
o.close()
M = np.array([[d[l][s][0] if s in d[l] and d[l][s][1] >= 20 else np.nan for s in S] for l in L])
o = io.open(D + "ZF_CHR4_ARRAY_REPRO_SPEARMAN.tsv", "w", encoding="utf8", newline="\n")
o.write(T.join(["a", "b", "n_array", "spearman"]) + "\n")
rs = []
for i, j in itertools.combinations(range(6), 2):
    ok = ~np.isnan(M[:, i]) & ~np.isnan(M[:, j])
    r = spearmanr(M[ok, i], M[ok, j])[0]
    rs.append(r)
    o.write(T.join([S[i], S[j], str(int(ok.sum())), "%.3f" % r]) + "\n")
o.write(T.join(["중앙값", "", "", "%.3f" % np.median(rs)]) + "\n")
o.write(T.join(["최소", "", "", "%.3f" % min(rs)]) + "\n")
o.write(T.join(["최대", "", "", "%.3f" % max(rs)]) + "\n")
o.close()
print("완료 array", len(L), "쌍", len(rs), "중앙 %.3f" % np.median(rs))
