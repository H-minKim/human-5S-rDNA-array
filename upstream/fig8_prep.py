# -*- coding: utf-8 -*-
# Fig 8 자료 — 자료원을 맞춘 두 짝 대조
#   (1) 암 : 같은 환자 · 같은 자료원 (종양 vs 그 환자의 정상)          7 짝
#   (2) 불멸화 : 같은 사람 · 같은 코호트(HPRC Year7) (LCL vs PBMC)   10 짝
#   지표 : 고메틸(ON) 유닛 안 ① 유전자부 β ② 스페이서 β ③ gene - sp
import io, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = io.open(os.path.join(HERE, "FIG8_PREP.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()

R = pd.read_csv(os.path.join(HERE, "NORM_VS_CANCER_SAMPLE.tsv"), sep="\t")
S = pd.read_csv(os.path.join(HERE, "ZSS_SAMPLE.tsv"), sep="\t")
R["군"] = R["sample"].map(dict(zip(S["sample"], S["군"])))
PW = pd.read_csv(os.path.join(HERE, "FIGX18_PAIRS.tsv"), sep="\t")

# ---------- (1) 암 짝 ----------
C = PW[PW.정상 != "HG008-N-D"].copy()      # HG008 은 정상이 둘이라 췌장(N-P)만
C["짝이름"] = C.짝
C["대조"] = "암"
C["a_이름"] = C.정상; C["b_이름"] = C.종양
for k in ("gene", "sp", "gene_sp"):
    C["a_" + k] = C["정상_" + k]; C["b_" + k] = C["종양_" + k]
P("[1] 암 짝 %d (같은 환자 · 같은 자료원)" % len(C))
P(C[["짝이름", "자료원", "호출방식", "a_gene", "b_gene", "a_gene_sp", "b_gene_sp"]].to_string(index=False))

# ---------- (2) 불멸화 짝 ----------
Y = R[R.군.isin(["정상 1차세포(PBMC)", "비암 세포주(LCL)"])].copy()
Y["person"] = Y["sample"].str.replace(r"_(LCL|PBMC)$", "", regex=True)
Y["t"] = np.where(Y["sample"].str.endswith("LCL"), "LCL", "PBMC")
rows = []
for p_, g in Y.groupby("person"):
    if set(g.t) != {"LCL", "PBMC"}: continue
    a = g[g.t == "PBMC"].iloc[0]; b = g[g.t == "LCL"].iloc[0]
    rows.append(dict(짝이름=p_, 대조="불멸화", 자료원="HPRC_Y7", 호출방식="같음",
                     a_이름=a["sample"], b_이름=b["sample"],
                     a_gene=a.고_gene, b_gene=b.고_gene,
                     a_sp=a.고_sp, b_sp=b.고_sp,
                     a_gene_sp=a.고_gene_sp, b_gene_sp=b.고_gene_sp))
I = pd.DataFrame(rows)
P()
P("[2] 불멸화 짝 %d (같은 사람 · HPRC Year7 · PBMC → LCL)" % len(I))
P(I[["짝이름", "a_gene", "b_gene", "a_gene_sp", "b_gene_sp"]].to_string(index=False))

PAIR = pd.concat([C[["짝이름", "대조", "자료원", "호출방식", "a_이름", "b_이름",
                     "a_gene", "b_gene", "a_sp", "b_sp", "a_gene_sp", "b_gene_sp"]], I],
                 ignore_index=True)
PAIR.to_csv(os.path.join(HERE, "FIG8_PAIRS.tsv"), sep="\t", index=False, float_format="%.5f")

# ---------- 검정 ----------
P()
P("[검정] 짝마다 b - a (암 = 종양 - 정상 · 불멸화 = LCL - PBMC)")
rows = []
for nm, g in PAIR.groupby("대조"):
    for k, lab in (("gene", "고메틸 유닛 안 유전자부 β"), ("sp", "고메틸 유닛 안 스페이서 β"),
                   ("gene_sp", "gene - sp")):
        d = (g["b_" + k] - g["a_" + k]).dropna()
        w = stats.wilcoxon(d)
        sg = stats.binomtest(int((d < 0).sum()), len(d), 0.5).pvalue
        rows.append(dict(대조=nm, 지표=lab, n=len(d), 중앙Δ=float(np.median(d)),
                         내려간짝=int((d < 0).sum()), 윌콕슨p=w.pvalue, 부호검정p=sg))
        P("   %-6s %-22s n %2d · 중앙 Δ %+.4f · 내려간 %d/%d · 윌콕슨 p %.4f · 부호검정 p %.4f"
          % (nm, lab, len(d), np.median(d), int((d < 0).sum()), len(d), w.pvalue, sg))
T = pd.DataFrame(rows)
T.to_csv(os.path.join(HERE, "FIG8_TESTS.tsv"), sep="\t", index=False, float_format="%.6g")

# ---------- 호출방식 같은 것만 (암 6 짝) ----------
P()
P("[민감도] 암에서 호출방식까지 같은 6 짝만")
c6 = PAIR[(PAIR.대조 == "암") & (PAIR.호출방식 == "같음")]
for k, lab in (("gene", "유전자부"), ("gene_sp", "gene - sp")):
    d = (c6["b_" + k] - c6["a_" + k]).dropna()
    P("   %-8s n %d · 중앙 Δ %+.4f · 내려간 %d/%d · 윌콕슨 p %.4f"
      % (lab, len(d), np.median(d), int((d < 0).sum()), len(d), stats.wilcoxon(d).pvalue))

# ---------- 자료원 차이(참고선) ----------
P()
P("[참고] 자료원 차이 — Y7 LCL 10 vs HPRC 본 코호트 280 (둘 다 LCL 유래)")
H = pd.read_csv(os.path.join(HERE, "HPRC280_SAMPLE.tsv"), sep="\t")
lcl = R[R.군 == "비암 세포주(LCL)"]
ref = []
for k, col in (("gene", "고_gene"), ("sp", "고_sp"), ("gene_sp", "고_gene_sp")):
    d = lcl[col].median() - H[col].median()
    u = stats.mannwhitneyu(lcl[col].dropna(), H[col].dropna())
    ref.append(dict(지표=k, Δ=float(d), p=u.pvalue))
    P("   %-8s Δ %+.4f · p %.4g" % (k, d, u.pvalue))
pd.DataFrame(ref).to_csv(os.path.join(HERE, "FIG8_SOURCEREF.tsv"), sep="\t", index=False, float_format="%.6g")
P()
P("[저장] FIG8_PAIRS.tsv · FIG8_TESTS.tsv · FIG8_SOURCEREF.tsv")
LOG.close()
