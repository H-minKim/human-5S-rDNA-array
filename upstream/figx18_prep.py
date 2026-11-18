# -*- coding: utf-8 -*-
# FigX18 data: what remains once the data source is controlled for.
#   The pooled comparison, 8 cancer against 30 normal, mixes data sources, so it is recounted using only pairs from the same patient and the same source.
import io, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = io.open(os.path.join(HERE, "FIGX18_PREP.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
COLS = ["고_gene", "고_sp", "고_gene_sp"]

R = pd.read_csv(os.path.join(HERE, "NORM_VS_CANCER_SAMPLE.tsv"), sep="\t")
S = pd.read_csv(os.path.join(HERE, "ZSS_SAMPLE.tsv"), sep="\t")
R["군"] = R["sample"].map(dict(zip(S["sample"], S["군"])))
src = {}
for f in ("CR_READUNIT.parquet", "CELL_READUNIT.parquet"):
    d = pd.read_parquet(os.path.join(HERE, f), columns=["sample", "src"])
    for s, g in d.groupby("sample"): src[s] = "/".join(sorted(g.src.unique()))
for s in R["sample"]:
    if s.endswith("_LCL"): src[s] = "HPRC_Y7_LCL"
    elif s.endswith("_PBMC"): src[s] = "HPRC_Y7_PBMC"
R["src"] = R["sample"].map(src)

# ---------- A: the map of data sources ----------
A = (R.pivot_table(index="src", columns="덩어리", values="sample", aggfunc="count")
     .fillna(0).astype(int).reset_index())
for k in ("암", "정상"):
    if k not in A.columns: A[k] = 0
A["짝가능"] = np.where((A.암 > 0) & (A.정상 > 0), "짝 가능", np.where(A.암 > 0, "암만", "정상만"))
A["정상_gene_sp"] = A.src.map(R[R.덩어리 == "정상"].groupby("src").고_gene_sp.median())
A = A.sort_values(["짝가능", "src"])
P("[A] data sources %d" % len(A))
P(A.to_string(index=False))

# ---------- B and C: the pairs ----------
PAIRS = [("COLO829-Tr", "COLO829-BLr", "COLO829_ONT_2024", "같음"),
         ("H1437-T", "H1437-N", "CASTLE_ONT_R10", "같음"),
         ("H2009-T", "H2009-N", "CASTLE_ONT_R10", "같음"),
         ("HCC1937-T", "HCC1937-N", "CASTLE_ONT_R10", "같음"),
         ("HCC1954-T", "HCC1954-N", "CASTLE_ONT_R10", "같음"),
         ("Hs578T-T", "Hs578T-N", "CASTLE_ONT_R10", "다름"),
         ("HG008-T", "HG008-N-P", "GIAB_HG008_ONT", "같음")]
V = R.set_index("sample")
rows = []
for t, n, sr, mode in PAIRS:
    r = dict(짝=t.rsplit("-", 1)[0].replace("-Tr", ""), 종양=t, 정상=n, 자료원=sr, 호출방식=mode)
    for c in COLS:
        k = c.replace("고_", "")
        r["정상_" + k] = V.loc[n, c]; r["종양_" + k] = V.loc[t, c]
        r["Δ" + k] = V.loc[t, c] - V.loc[n, c]
    rows.append(r)
B = pd.DataFrame(rows)
P()
P("[B, C] pairs %d (HCC1395 has only 1 tumor read, so it forms no pair)" % len(B))
P(B[["짝", "자료원", "호출방식", "정상_gene_sp", "종양_gene_sp", "Δgene_sp", "Δgene", "Δsp"]].to_string(index=False))

# ---------- D: pooled against paired ----------
def sgn(v): return stats.binomtest(int((np.array(v) < 0).sum()), len(v), 0.5).pvalue
can = R[R.덩어리 == "암"]; nor = R[R.덩어리 == "정상"]
cas = R[R.src == "CASTLE_ONT_R10"]
D = []
u = stats.mannwhitneyu(can.고_gene_sp, nor.고_gene_sp)
D.append(dict(판="pooled 암 8 vs 정상 30", 자료원="섞임", n="8 vs 30",
              Δ=can.고_gene_sp.median() - nor.고_gene_sp.median(), p=u.pvalue, 검정="Mann-Whitney"))
u = stats.mannwhitneyu(cas[cas.덩어리 == "암"].고_gene_sp, cas[cas.덩어리 == "정상"].고_gene_sp)
D.append(dict(판="CASTLE 안 5 vs 5", 자료원="한 자료원", n="5 vs 5",
              Δ=cas[cas.덩어리 == "암"].고_gene_sp.median() - cas[cas.덩어리 == "정상"].고_gene_sp.median(),
              p=u.pvalue, 검정="Mann-Whitney"))
D.append(dict(판="짝 7 (같은 환자·같은 자료원)", 자료원="통제", n="7 짝",
              Δ=B.Δgene_sp.median(), p=stats.wilcoxon(B.Δgene_sp).pvalue, 검정="윌콕슨"))
b6 = B[B.호출방식 == "같음"]
D.append(dict(판="짝 6 (호출방식까지 같은 것)", 자료원="통제", n="6 짝",
              Δ=b6.Δgene_sp.median(), p=stats.wilcoxon(b6.Δgene_sp).pvalue, 검정="윌콕슨"))
D = pd.DataFrame(D)
D["부호검정p"] = [np.nan, np.nan, sgn(B.Δgene_sp), sgn(b6.Δgene_sp)]
P()
P("[D] the same metric (gene-sp within high-methylation units) counted with only the version changed")
P(D.to_string(index=False))

# ---------- E: the cancer pairs against the immortalisation pairs ----------
IM = pd.read_csv(os.path.join(HERE, "FIGX17_IMM.tsv"), sep="\t")
IM["Δ"] = IM.LCL - IM.PBMC
P()
P("[E] cancer pairs 7 (data source controlled) vs immortalization pairs 10 (same person, same data source)")
P("   cancer           median Δ %+.4f, decreased %d/%d, Wilcoxon p %.4f, sign test p %.4f"
  % (B.Δgene_sp.median(), int((B.Δgene_sp < 0).sum()), len(B),
     stats.wilcoxon(B.Δgene_sp).pvalue, sgn(B.Δgene_sp)))
P("   immortalization median Δ %+.4f, decreased %d/%d, Wilcoxon p %.4f, sign test p %.4f"
  % (IM.Δ.median(), int((IM.Δ < 0).sum()), len(IM),
     stats.wilcoxon(IM.Δ).pvalue, sgn(IM.Δ)))
P("   the constraint (immortalization) has a smaller p than the evidence (cancer) = it is firmer")
IM.to_csv(os.path.join(HERE, "FIGX18_IMM.tsv"), sep="\t", index=False, float_format="%.5f")

# ---------- F: the ceiling set by the number of pairs ----------
rows = []
for n in range(3, 21):
    rows.append(dict(n=n, 부호검정_최소p=stats.binomtest(0, n, 0.5).pvalue,
                     윌콕슨_최소p=stats.wilcoxon(-np.arange(1., n + 1)).pvalue))
F = pd.DataFrame(rows)
P()
P("[F] smallest attainable p with n pairs (assuming all in the same direction)")
P(F.to_string(index=False))
P("   -> with 7 pairs the sign test minimum is %.4f, Wilcoxon minimum %.4f; the 39-cell correction threshold 0.0013 is unattainable"
  % (F[F.n == 7].부호검정_최소p.iloc[0], F[F.n == 7].윌콕슨_최소p.iloc[0]))

for nm, df in (("SOURCE", A), ("PAIRS", B), ("COMPARE", D), ("MINP", F)):
    df.to_csv(os.path.join(HERE, "FIGX18_%s.tsv" % nm), sep="\t", index=False, float_format="%.6g")
P()
P("[saved] FIGX18_SOURCE/PAIRS/COMPARE/MINP/IMM.tsv, FIGX18_PREP.txt")
LOG.close()
