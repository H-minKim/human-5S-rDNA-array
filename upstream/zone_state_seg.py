# -*- coding: utf-8 -*-
# 군 × array 구역 × 유닛 상태(저/고) × 구획(스페이서·유전자부·차) 한 표로
# zone = the canonical zone definition, the zone() of 01_code/v36_cancer_by_zone.py
#   군   = 정상 조직(환자) · 정상 1차세포(PBMC) · 비암 세포주(LCL·GM12878·짝지은 정상) · 암 조직 · 암 세포주
#   묻는 것 : 암 조직과 암 세포주가 함께 정상과 다른 자리가 있나
import io, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
LOG = io.open(os.path.join(HERE, "ZONE_STATE_SEG.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
THR = 0.508

CR = pd.read_parquet(os.path.join(HERE, "CR_READUNIT.parquet"))
CL = pd.read_parquet(os.path.join(HERE, "CELL_READUNIT.parquet"))
TI = pd.read_parquet(os.path.join(HERE, "TISSUE_READUNIT.parquet"))
cols = [c for c in CR.columns if c in CL.columns and c in TI.columns]
D = pd.concat([CR[cols], CL[cols], TI[cols]], ignore_index=True)

GRP = {}
for s in ("HG008-N-P", "HG008-N-D"): GRP[s] = "정상 조직(환자)"
GRP["HG008-T"] = "암 조직(환자)"
for s in ("COLO829-Tr", "H1437-T", "H2009-T", "HCC1937-T", "HCC1954-T", "Hs578T-T", "K562r"):
    GRP[s] = "암 세포주"
for s in ("COLO829-BLr", "H1437-N", "H2009-N", "HCC1937-N", "HCC1954-N", "Hs578T-N",
          "GM12878", "GM12878FS"): GRP[s] = "비암 세포주"
for s in TI[TI.tissue == "LCL"]["sample"].unique(): GRP[s] = "비암 세포주(LCL)"
for s in TI[TI.tissue == "PBMC"]["sample"].unique(): GRP[s] = "정상 1차세포(PBMC)"
D = D[D["sample"].isin(GRP)].copy()
D["군"] = D["sample"].map(GRP)
D["state"] = np.where(D.b_call >= THR, "고메틸", "저메틸")
D["apx"] = D.N - D.ad + 1 if "ad" in D.columns else D.N - D.ap + 1
D["adx"] = D.ap            # array_pos 1 = 원위 (정본 규칙)
D["apx"] = D.N - D.ap + 1  # 근위 index

def zone(ap, ad):
    if ap <= ad:
        if ap <= 3:  return "Z1 근1-3"
        if ap <= 10: return "Z2 근4-10"
        if ap <= 17: return "Z3 근11-17"
        return "ZI 내부"
    if ad <= 2:  return "Z7 원1-2"
    if ad <= 8:  return "Z6 원3-8"
    if ad <= 10: return "Z5 원9-10"
    if ad <= 22: return "Z4 원11-22"
    return "ZI 내부"
D["Z"] = [zone(a, d) for a, d in zip(D.apx.values, D.adx.values)]
ZONES = ["Z1 근1-3", "Z2 근4-10", "Z3 근11-17", "ZI 내부", "Z4 원11-22", "Z5 원9-10", "Z6 원3-8", "Z7 원1-2"]
ORD = ["정상 조직(환자)", "정상 1차세포(PBMC)", "비암 세포주(LCL)", "비암 세포주", "암 조직(환자)", "암 세포주"]
P("read×유닛 %s · 시료 %d" % (format(len(D), ","), D["sample"].nunique()))
for g in ORD:
    v = D[D.군 == g]
    P("   %-18s 시료 %2d · 유닛 %s" % (g, v["sample"].nunique(), format(len(v), ",")))

def wb(x, k):
    n = x["n_" + k].sum()
    return ((x["b_" + k] * x["n_" + k]).sum() / n) if n else np.nan

# ---------------- 시료 × 구역 × 상태 ----------------
rows = []
for (s, z, st), g in D.groupby(["sample", "Z", "state"]):
    if len(g) < 20: continue
    sp = wb(g, "sp"); ge = wb(g, "gene")
    rows.append(dict(sample=s, 군=GRP[s], 구역=z, 상태=st, n유닛=len(g),
                     sp=sp, gene=ge, gene_sp=ge - sp, ups=wb(g, "ups")))
S = pd.DataFrame(rows)
S.to_csv(os.path.join(HERE, "ZSS_SAMPLE.tsv"), sep="\t", index=False, float_format="%.5f")

# ---------------- 군 요약 (시료 중앙값) ----------------
P()
P("[표] 군 × 구역 × 상태 — 값은 시료 중앙값 (n = 그 칸에 들어간 시료 수)")
for st in ("저메틸", "고메틸"):
    P()
    P("=== %s 유닛 안 ===" % st)
    P("%-10s %-18s %5s %8s %8s %9s" % ("구역", "군", "시료", "sp", "gene", "gene-sp"))
    for z in ZONES:
        for g in ORD:
            v = S[(S.구역 == z) & (S.군 == g) & (S.상태 == st)]
            if not len(v): continue
            P("%-10s %-18s %5d %8.4f %8.4f %+9.4f"
              % (z, g, len(v), v.sp.median(), v.gene.median(), v.gene_sp.median()))
        P("-" * 62)

# ---------------- 암 공통 : 두 암 군이 함께 정상과 다른가 ----------------
P()
P("[검정] 암 조직(환자) 과 암 세포주 가 함께 정상 쪽과 다른 자리")
P("   정상 쪽 = 정상 조직(환자) + 정상 1차세포(PBMC) + 비암 세포주(LCL) + 비암 세포주")
NORM = ["정상 조직(환자)", "정상 1차세포(PBMC)", "비암 세포주(LCL)", "비암 세포주"]
res = []
for st in ("저메틸", "고메틸"):
    for z in ZONES:
        for col in ("sp", "gene", "gene_sp"):
            n = S[(S.구역 == z) & (S.상태 == st) & (S.군.isin(NORM))][col].dropna()
            c = S[(S.구역 == z) & (S.상태 == st) & (S.군 == "암 세포주")][col].dropna()
            t = S[(S.구역 == z) & (S.상태 == st) & (S.군 == "암 조직(환자)")][col].dropna()
            if len(n) < 8 or len(c) < 4 or len(t) < 1: continue
            u = stats.mannwhitneyu(c, n)
            zt = (t.iloc[0] - n.mean()) / n.std(ddof=1) if n.std(ddof=1) > 0 else np.nan
            res.append(dict(구역=z, 상태=st, 지표=col, n정상=len(n), n암세포주=len(c),
                            정상중앙=n.median(), 암세포주중앙=c.median(), 암조직=t.iloc[0],
                            Δ암세포주=c.median() - n.median(), 암조직z=zt, p세포주=u.pvalue))
R = pd.DataFrame(res)
R["q"] = np.minimum(1, R.p세포주 * len(R))
# 같은 방향 + 세포주 유의 + 조직도 |z|>1
R["암공통"] = ((R.p세포주 < 0.05) & (np.sign(R.Δ암세포주) == np.sign(R.암조직z)) & (R.암조직z.abs() > 1))
R = R.sort_values("p세포주")
P()
P(R.round(4).to_string(index=False))
R.to_csv(os.path.join(HERE, "ZSS_TEST.tsv"), sep="\t", index=False, float_format="%.5f")
P()
P("암 조직과 암 세포주가 같은 방향으로 벗어난 칸 (세포주 p<0.05 · 조직 |z|>1 · 방향 일치)")
hit = R[R.암공통]
if len(hit):
    P(hit[["구역", "상태", "지표", "정상중앙", "암세포주중앙", "암조직", "Δ암세포주", "암조직z", "p세포주", "q"]]
      .round(4).to_string(index=False))
else:
    P("   없음")
P("   (검정 %d 칸 · Bonferroni q 열 참고 — q<0.05 인 칸 %d)" % (len(R), int((R.q < 0.05).sum())))

# ---------------- 조직 비교 (같은 사람 LCL vs PBMC) ----------------
P()
P("[tissue] LCL and PBMC of the same 10 individuals : paired test")
T = TI.copy()
T["state"] = np.where(T.b_call >= THR, "저메틸", "고메틸")
T["state"] = np.where(T.b_call >= THR, "고메틸", "저메틸")
per = []
for (pid, ti), g in T.groupby(["person", "tissue"]):
    sp = wb(g, "sp"); ge = wb(g, "gene")
    per.append(dict(person=pid, tissue=ti, n=len(g), b=g.n_meth.sum() / g.n_cpg.sum(),
                    sp=sp, gene=ge, gene_sp=ge - sp, 저메틸몫=(g.state == "저메틸").mean()))
PT = pd.DataFrame(per).pivot(index="person", columns="tissue")
PT.to_csv(os.path.join(HERE, "ZSS_TISSUE_PERSON.tsv"), sep="\t", float_format="%.5f")
for col, nm in (("b", "전체 β"), ("gene_sp", "gene−sp"), ("gene", "유전자부 β"),
                ("sp", "스페이서 β"), ("저메틸몫", "저메틸 유닛 몫")):
    a = PT[(col, "LCL")]; b = PT[(col, "PBMC")]
    d = a - b
    w = stats.wilcoxon(a, b)
    P("   %-14s LCL %+.4f · PBMC %+.4f · Δ(LCL−PBMC) 중앙 %+.4f · LCL 더 낮음 %d/%d · 윌콕슨 p %.4f"
      % (nm, a.median(), b.median(), d.median(), int((d < 0).sum()), len(d), w.pvalue))
P()
P("   ※ LCL 은 EBV 로 불멸화한 세포주 이고 PBMC 는 1 차세포다. 같은 사람에서 비교하므로")
P("     개체 차이는 지워진다. 조직 차이인지 불멸화 효과 인지는 이 자료만으로 못 가른다.")
P()
P("[저장] ZSS_SAMPLE.tsv · ZSS_TEST.tsv · ZSS_TISSUE_PERSON.tsv")
LOG.close()
