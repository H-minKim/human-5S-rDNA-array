# -*- coding: utf-8 -*-
# Builds the input table for Figure 3G, FIG3_G_TYPEDIST.tsv.
#
#   What it produces:
#     one point per haplotype, with the mean distance to the members of its own array-type group
#     and the mean distance to the other groups written on the same row.
#
#   Cohort: H3 = H1 intersect H2, because the claim uses both the variant-based grouping and the methylation.
#   Gate: L2_HAP_UNIT filtered to n_read >= 5, then only haplotypes with a value in all five zones.
#   Distance: the mean |delta| between two haplotypes over the five zone betas.
#   Excluded: pairs from the same individual, which are not independent.
#   Groups: cluster size at least 5, and the cluster index not 0.
#   Random numbers: default_rng(20261104), with 2,000 permutations of the group labels, printed to the terminal and not written to the table.
#
#   Four inputs; the absolute paths are assembled from ROOT and GRP below.
#     1) <ROOT>/02_tables/TIERS_H0H3_SLOTS_2026-10-05.tsv, the tiers, used to select H3
#     ② <GRP>   5S-rDNA-paper-en/tables/X12_HAP_GROUPS.tsv
#               the groups of Figure 2M (single linkage, Jaccard >= 0.90, cluster size >= 5)
#     3) <ROOT>/02_tables/L2_HAP_UNIT.parquet, haplotype x unit methylation; only seven columns are read
#     4) <ROOT>/02_tables/POP_COHORT.tsv, the superpopulation and the unit count
#
#   Written to <FIG>/FIG3_G_TYPEDIST.tsv.
#             Setting FIG3G_OUT writes it elsewhere, so the table can be reproduced without touching the published input.
#
#   Read by panel G of fig3_main_draw_en.py and the G_type_distance sheet of fig3_sourcedata_en.py.
#   The p value printed in the figure is a paired Wilcoxon, not the permutation p.
import hashlib
import io
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.dont_write_bytecode = True
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
S = ""
NL = chr(10)
ROOT = os.environ["DATA_ROOT"]
FIG = ROOT + "/03_figures/paper_Fig"
GRP = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tables", "X12_HAP_GROUPS.tsv")
md5 = lambda p: hashlib.md5(io.open(p, "rb").read()).hexdigest()[:12]
rng = np.random.default_rng(20261104)

TI = pd.read_csv(ROOT + "/02_tables/TIERS_H0H3_SLOTS_2026-10-05.tsv", sep="\t", dtype=str)
TI["hap"] = TI["sample"].astype(str) + "_" + TI.slot.astype(str)
H3 = set(TI[TI.tier == "H3"].hap)

G = pd.read_csv(GRP, sep="\t")
G.columns = ["clu", "hap", "sample", "slot", "clu_size", "n_token"]
big = G[(G.clu_size >= 5) & (G.clu != 0)]

L2 = pd.read_parquet(ROOT + "/02_tables/L2_HAP_UNIT.parquet",
                     columns=["sample", "slot", "ap", "ad", "N", "n_read", "b_mean"])
L2 = L2[L2.n_read >= 5].copy()
L2["hap"] = L2["sample"].astype(str) + "_" + L2.slot.astype(str)
Z = [("P1_3", (L2.ap <= L2.ad) & (L2.ap <= 3)),
     ("P4_10", (L2.ap <= L2.ad) & (L2.ap >= 4) & (L2.ap <= 10)),
     ("INT", ((L2.ap <= L2.ad) & (L2.ap > 10)) | ((L2.ad < L2.ap) & (L2.ad > 21))),
     ("D11_21", (L2.ad < L2.ap) & (L2.ad >= 11) & (L2.ad <= 21)),
     ("D1_10", (L2.ad < L2.ap) & (L2.ad <= 10))]
V = pd.DataFrame(index=sorted(L2.hap.unique()))
for nm, msk in Z:
    V[nm] = L2[msk].groupby("hap").b_mean.mean()
V = V.dropna()
P = pd.read_csv(ROOT + "/02_tables/POP_COHORT.tsv", sep="\t")
P["hap"] = P["sample"].astype(str) + "_" + P.slot.astype(str)
P = P.set_index("hap")[["sample", "sp", "N"]]
D = big.set_index("hap").join(V, how="inner").join(P[["sp", "N"]], how="inner")
D = D.rename(columns={"sample": "ind"})
D = D[[h in H3 for h in D.index]]
COLS = ["P1_3", "P4_10", "INT", "D11_21", "D1_10"]
X = D[COLS].to_numpy(float)
clu = D.clu.astype(int).to_numpy()
ind = D.ind.astype(str).to_numpy()
hap = np.array(D.index)
n = len(hap)
print("%scohort H3, haplotypes %d, individuals %d, groups %d%s"
      % (S, n, len(set(ind)), len(set(clu)), S))

DIST = np.abs(X[:, None, :] - X[None, :, :]).mean(axis=2)
OK = (ind[:, None] != ind[None, :])
np.fill_diagonal(OK, False)
same = clu[:, None] == clu[None, :]
rows = []
for i in range(n):
    w = DIST[i][same[i] & OK[i]]
    b = DIST[i][(~same[i]) & OK[i]]
    if len(w) == 0 or len(b) == 0:
        continue
    rows.append(dict(hap=hap[i], 개체=ind[i], 무리=clu[i], 무리크기=int((clu == clu[i]).sum()),
                     초집단=D.sp.iloc[i], 유닛수=int(D.N.iloc[i]),
                     n_군안=len(w), 군안거리=round(float(w.mean()), 5),
                     n_군밖=len(b), 군밖거리=round(float(b.mean()), 5),
                     차=round(float(w.mean() - b.mean()), 5)))
T = pd.DataFrame(rows)
print("  points (haplotypes) %d, closer within group %d (%.1f%%)"
      % (len(T), int((T.차 < 0).sum()), 100.0 * (T.차 < 0).mean()))
st, pw = wilcoxon(T.군안거리, T.군밖거리)
print("  paired Wilcoxon p = %.3g, median within-pair difference %+.4f"
      % (pw, float(np.median(T.차))))
print("  median: within group %.4f, outside group %.4f"
      % (T.군안거리.median(), T.군밖거리.median()))

# the null from permuting group labels, on the haplotype-level statistic
obs = float(np.median(T.차))
cnt = 0
NP = 2000
idx = {h: k for k, h in enumerate(hap)}
for _ in range(NP):
    p = clu.copy()
    rng.shuffle(p)
    sm = p[:, None] == p[None, :]
    dif = []
    for i in range(n):
        w = DIST[i][sm[i] & OK[i]]
        b = DIST[i][(~sm[i]) & OK[i]]
        if len(w) and len(b):
            dif.append(w.mean() - b.mean())
    if np.median(dif) <= obs:
        cnt += 1
print("  group label permutations %d: p = %.4f (floor %.4f)" % (NP, (cnt + 1) / (NP + 1), 1 / (NP + 1)))

OUT = os.environ.get("FIG3G_OUT") or os.path.join(FIG, "FIG3_G_TYPEDIST.tsv")
HDR = [
    "#  Fig 3 (G) 입력 — array-type 군 안 ↔ 군 밖 메틸 거리",
    "#  코호트 : H3 = H1 ∩ H2 (변이 참조 ∩ 메틸 참조) · 하플로타입 %d · 개체 %d · 군 14" % (
        len(T), T["개체"].nunique()),
    "#  군 : X12_HAP_GROUPS.tsv (Fig 2M · single linkage Jaccard >= 0.90 · 무리 크기 >= 5)",
    "#  메틸 : L2_HAP_UNIT.parquet (n_read >= 5) 의 다섯 구역 평균 beta",
    "#  거리 : 두 하플로타입의 다섯 구역 |delta beta| 평균 · 같은 개체 짝은 뺀다",
    "#  점 하나 = 한 하플로타입 (군 안 평균 거리 · 군 밖 평균 거리)",
]
#  No provenance header is written inside the table. The four figure input tables in this folder (FIG3_G_TYPEDIST,
#  FIG5V2_PAIRZONE5, FIG8_PAIRS and NORM_VS_CANCER_SAMPLE) carry no # lines, and
#  the two scripts that read them do so without comment="#"
#  (fig3_main_draw_en.py:103 and fig3_sourcedata_en.py:36 call pd.read_csv(sep="\t") plainly).
#  A # line would therefore be read as the column names and break the figure.
#  The header is written to a side file instead, and the table begins with its column names.
with io.open(OUT, "w", encoding="utf-8", newline=NL) as f:
    T.to_csv(f, sep="\t", index=False, lineterminator=NL)
with io.open(OUT + ".about.txt", "w", encoding="utf-8", newline=NL) as f:
    f.write(NL.join(HDR).replace("#  ", "") + NL)
print()
print("%stable %s, %s, %d lines%s" % (S, os.path.basename(OUT), md5(OUT), len(T), S))
print("   sidecar file: %s.about.txt" % os.path.basename(OUT))
print(S + "done" + S)
