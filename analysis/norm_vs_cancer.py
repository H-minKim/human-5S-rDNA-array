# -*- coding: utf-8 -*-
# 물음 정상(조직+1차세포+세포주) vs 암(세포주+조직) 을 두 덩어리로 합쳐 보면 결과가 뭔가
#   정상 = 정상 조직(2) + PBMC(10) + LCL(10) + 비암 세포주(8) = 30 시료
#   암   = 암 세포주(8) = 8 시료
#   정정 : 전에 「암 세포주 7 + 암 조직 1」 로 적었는데 그 1(HG008-T)도
#     세포주 다 (환자유래 PDAC · ATCC CRL-3734 · PMC12267650). 시료 집합은 그대로 8 이라
#     이 글이 내는 수는 하나도 안 바뀐다 — 바뀐 것은 이름뿐이다.
#   판정 : 시료 단위 Mann-Whitney + AUC · BH q
import io, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats
# 표를 읽고 쓰는 자리 — 저장소에서는 표가 코드 옆에 없으므로 환경변수로 받는다 (2026-10-21)
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
LOG = io.open(os.path.join(HERE, "NORM_VS_CANCER.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
S = pd.read_csv(os.path.join(HERE, "ZSS_SAMPLE.tsv"), sep="\t")
ZONES = ["Z1 근1-3", "Z2 근4-10", "Z3 근11-17", "ZI 내부", "Z4 원11-22", "Z5 원9-10", "Z6 원3-8", "Z7 원1-2"]
CAN = ["암 세포주"]   # 옛 라벨 "암 조직(환자)" 는 없어졌다 (시료 집합은 같은 8)
S["덩어리"] = np.where(S.군.isin(CAN), "암", "정상")
P("[군 합치기] 암 = 암 세포주(8) · 정상 = 정상 조직 + PBMC + LCL + 비암 세포주")
P("   시료 수 : %s" % S.drop_duplicates("sample").덩어리.value_counts().to_dict())

def bh(ps):
    ps = np.asarray(ps, float); ok = ~np.isnan(ps); q = np.full(len(ps), np.nan)
    p = ps[ok]; n = len(p)
    if n == 0: return q
    o = np.argsort(p); r = np.empty(n); r[o] = np.arange(1, n + 1)
    qq = np.minimum.accumulate((p * n / r)[o][::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(qq, 1); q[ok] = out
    return q

rows = []
for state in ("고메틸", "저메틸"):
    for col in ("gene", "sp", "gene_sp"):
        for z in ZONES:
            d = S[(S.상태 == state) & (S.구역 == z)]
            a = d[d.덩어리 == "암"][col].dropna(); b = d[d.덩어리 == "정상"][col].dropna()
            r = dict(상태=state, 지표=col, 구역=z, n암=len(a), n정상=len(b),
                     암=a.median() if len(a) else np.nan, 정상=b.median() if len(b) else np.nan)
            if len(a) >= 3 and len(b) >= 5:
                u = stats.mannwhitneyu(a, b)
                r["p"] = u.pvalue; r["AUC"] = u.statistic / (len(a) * len(b))
            else:
                r["p"] = np.nan; r["AUC"] = np.nan
            rows.append(r)
T = pd.DataFrame(rows)
T["q"] = bh(T.p.values)

for state in ("고메틸", "저메틸"):
    for col, nm in (("sp", "스페이서 β"), ("gene", "유전자부 β"), ("gene_sp", "gene - sp")):
        P()
        P("=" * 92)
        P("[%s 유닛 안 · %s]" % (state, nm))
        P("=" * 92)
        P("   %-11s %5s %6s %9s %9s %8s %8s %8s" % ("구역", "n암", "n정상", "암", "정상", "AUC", "p", "BH q"))
        for _, r in T[(T.상태 == state) & (T.지표 == col)].iterrows():
            P("   %-11s %5d %6d %9s %9s %8s %8s %8s"
              % (r.구역, r.n암, r.n정상,
                 "%.4f" % r.암 if r.암 == r.암 else "-",
                 "%.4f" % r.정상 if r.정상 == r.정상 else "-",
                 "%.2f" % r.AUC if r.AUC == r.AUC else "-",
                 "%.3f" % r.p if r.p == r.p else "불가",
                 "%.3f" % r.q if r.q == r.q else "-"))

P()
P("=" * 92)
P("[한 줄 요약]")
P("=" * 92)
sig = T[T.p < 0.05].copy()
P("   검정한 칸 %d · nominal p<0.05 %d · BH q<0.05 %d"
  % (int(T.p.notna().sum()), len(sig), int((T.q < 0.05).sum())))
P()
P("   유의한 칸(nominal)")
for _, r in sig.sort_values("p").iterrows():
    P("      %-6s %-9s %-11s 암 %+.4f vs 정상 %+.4f · AUC %.2f · p %.4f · q %.3f"
      % (r.상태, r.지표, r.구역, r.암, r.정상, r.AUC, r.p, r.q))
P()
P("   q<0.05 로 살아남는 칸")
sq = T[T.q < 0.05]
if len(sq):
    for _, r in sq.sort_values("q").iterrows():
        P("      %-6s %-9s %-11s 암 %+.4f vs 정상 %+.4f · AUC %.2f · q %.4f"
          % (r.상태, r.지표, r.구역, r.암, r.정상, r.AUC, r.q))
else:
    P("      없음")
T.to_csv(os.path.join(HERE, "NORM_VS_CANCER.tsv"), sep="\t", index=False, float_format="%.5f")

# array 전체(구역 안 가름)
P()
P("=" * 92)
P("[구역을 안 가르고 array 전체로] 시료 단위")
P("=" * 92)
A = pd.read_parquet(os.path.join(HERE, "CR_READUNIT.parquet"))
B = pd.read_parquet(os.path.join(HERE, "CELL_READUNIT.parquet"))
C = pd.read_parquet(os.path.join(HERE, "TISSUE_READUNIT.parquet"))
cols = [c for c in A.columns if c in B.columns and c in C.columns]
D = pd.concat([A[cols], B[cols], C[cols]], ignore_index=True)
GRP = dict(zip(S["sample"], S.군)); D = D[D["sample"].isin(GRP)].copy()
D["덩어리"] = np.where(D["sample"].map(GRP).isin(CAN), "암", "정상")
D["state"] = np.where(D.b_call >= 0.508, "고메틸", "저메틸")
def wb(g, k):
    n = g["n_" + k].sum()
    return ((g["b_" + k] * g["n_" + k]).sum() / n) if n else np.nan
rows2 = []
for s, g in D.groupby("sample"):
    hi = g[g.state == "고메틸"]
    rows2.append(dict(sample=s, 덩어리=g.덩어리.iloc[0],
                      전체b=g.n_meth.sum() / g.n_cpg.sum(),
                      저메틸몫=(g.state == "저메틸").mean(),
                      고_sp=wb(hi, "sp"), 고_gene=wb(hi, "gene"),
                      고_gene_sp=wb(hi, "gene") - wb(hi, "sp")))
R = pd.DataFrame(rows2)
P("   %-14s %10s %10s %8s %8s" % ("지표", "암(n8)", "정상(n30)", "AUC", "p"))
for c, nm in (("전체b", "전체 β"), ("저메틸몫", "저메틸 유닛 몫"), ("고_sp", "고메틸 sp"),
              ("고_gene", "고메틸 gene"), ("고_gene_sp", "고메틸 gene-sp")):
    a = R[R.덩어리 == "암"][c].dropna(); b = R[R.덩어리 == "정상"][c].dropna()
    u = stats.mannwhitneyu(a, b)
    P("   %-14s %10.4f %10.4f %8.2f %8.4f %s"
      % (nm, a.median(), b.median(), u.statistic / (len(a) * len(b)), u.pvalue,
         "" if u.pvalue < 0.05 else ""))
R.to_csv(os.path.join(HERE, "NORM_VS_CANCER_SAMPLE.tsv"), sep="\t", index=False, float_format="%.5f")
P()
P("   ※ 전체 β 는 호출 척도(ML)를 타므로 결론에 쓰지 않는다(FigX15 B).")
P()
P("[저장] NORM_VS_CANCER.tsv · NORM_VS_CANCER_SAMPLE.tsv")
LOG.close()
