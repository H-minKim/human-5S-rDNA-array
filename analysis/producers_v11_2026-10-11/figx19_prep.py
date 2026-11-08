# -*- coding: utf-8 -*-
# FigX19 자료 — 근위 변이집합이 유닛 메틸에 미치는 영향 (유닛 단위)
#   자료 : L2_HAP_UNIT.parquet (하플로타입 × 유닛 · HPRC_ONT_submissions)
#   규칙 : ap = 근위 끝에서 센 자리(ap 1 = 근위 말단) · ad = 원위 끝에서 센 자리 (실측 확인함)
#   위치가 메틸을 지배하므로 모든 비교는 위치 띠 안에서 · 개체 색조는 하플로타입 안에서 지운다
import io, os, sys, itertools
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import itertools
import numpy as np, pandas as pd
from scipy import stats

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    이 글이 읽는 자료의 뿌리는 README.md 의
    "Running this code" 절에 적혀 있다."""
    import os as _o, sys as _s
    _x = _o.environ.get(_v)
    if not _x:
        _s.exit(
            "environment variable %s is not set; it is required by %s"
            % (_v, _o.path.basename(__file__))
            + "\n(a module stops the script that imports it)."
            + "\nThe deposited copy has no default; see README.md"
            + " section \"Running this code\"."
        )
    return _x

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = io.open(os.path.join(HERE, "FIGX19_PREP.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
rng = np.random.default_rng(19)
T = os.path.join(_need("RDNA_V11_ROOT"), "02_tables")
MINREAD = 5

#  원본은 `T = "$DATA_ROOT/02_tables/"` 로 끝에 `/` 가 있었고
#    여기서 문자열을 이어붙였다. 경로치환이 그 `/` 를 없애 `…02_tablesL2_HAP_UNIT.parquet` 가 됐다.
L = pd.read_parquet(os.path.join(T, "L2_HAP_UNIT.parquet"), columns=[
    "sample", "slot", "ap", "ad", "N", "side", "f_hi", "n_read", "varlist", "n_var",
    "n_var_common", "n_var_rare", "n_var_gene", "nv_coding", "vsig", "vsig_n",
    "b_gene_mean", "b_sp_mean", "b_mean", "data_source", "pop", "super_pop",
    "f_hi_thr"])
P("[입력] L2_HAP_UNIT.parquet · %s 행 · 자료원 %s" % (format(len(L), ","), set(L.data_source)))
L["hap"] = L["sample"].astype(str) + "|" + L.slot.astype(str)
P("   개체 %d · 하플로타입 %d · 유닛 %s"
  % (L["sample"].nunique(), L.hap.nunique(), format(len(L), ",")))
L = L[L.n_read >= MINREAD].copy()
L["vl"] = L.varlist.fillna("-").astype(str)
P("   read >= %d 만 : 유닛 %s · 개체 %d · 하플로 %d"
  % (MINREAD, format(len(L), ","), L["sample"].nunique(), L.hap.nunique()))
PROV = [dict(항목="바탕 자료", 자료="L2_HAP_UNIT.parquet", 자료원="HPRC_ONT_submissions",
             개체=L["sample"].nunique(), 하플로=L.hap.nunique(), 유닛=len(L),
             비고="read>=%d 인 유닛만" % MINREAD)]

# ---------- A : 위치가 먼저다 ----------
P()
P("[A] 위치가 메틸을 지배한다 — 근위 자리별 f_hi (그래서 아래 비교는 전부 위치 띠 안에서 한다)")
A = (L[L.ap <= 20].groupby("ap")
     .agg(유닛=("f_hi", "size"), 하플로=("hap", "nunique"), f_hi중앙=("f_hi", "median"),
          f_hi평균=("f_hi", "mean"), 변이수평균=("n_var", "mean")).reset_index())
P(A.to_string(index=False))
A.to_csv(os.path.join(HERE, "FIGX19_A_POS.tsv"), sep="\t", index=False, float_format="%.4f")

# 위치 띠 : 근위 1-3 / 4-10 (정본 구역 Z1 · Z2) · 견줌용 내부 11-30
def band(ap):
    if ap <= 3: return "근위 1-3"
    if ap <= 10: return "근위 4-10"
    if ap <= 30: return "근위 11-30"
    return "그 밖"
L["띠"] = [band(a) for a in L.ap]
PROX = L[L.띠.isin(["근위 1-3", "근위 4-10", "근위 11-30"])].copy()
P()
P("   근위 30 유닛 안 : 유닛 %s · 하플로 %d · 개체 %d"
  % (format(len(PROX), ","), PROX.hap.nunique(), PROX["sample"].nunique()))

# ---------- B : 변이 개수 ↔ 메틸 (위치 통제 전/후) ----------
P()
P("[B] 변이 개수 와 f_hi — 통제 전/후")
rows = []
r0, p0 = stats.spearmanr(PROX.n_var, PROX.f_hi)
rows.append(dict(구분="통제 전 (근위 30 전부)", n=len(PROX), rho=r0, p=p0))
P("   통제 전 : rho %+.3f (p %.3g · 유닛 %d)" % (r0, p0, len(PROX)))
for b, g in PROX.groupby("띠"):
    r, p = stats.spearmanr(g.n_var, g.f_hi)
    rows.append(dict(구분="위치 띠 안 : " + b, n=len(g), rho=r, p=p))
    P("   %-12s rho %+.3f (p %.3g · 유닛 %d · 하플로 %d)" % (b, r, p, len(g), g.hap.nunique()))
# 하플로타입 안에서 (개체 색조 제거) : 하플로마다 rho 를 구해 부호검정
hr = []
for h, g in PROX.groupby("hap"):
    if g.n_var.nunique() > 2 and len(g) >= 15:
        hr.append(stats.spearmanr(g.n_var, g.f_hi)[0])
hr = np.array([x for x in hr if x == x])
rows.append(dict(구분="하플로타입 안 rho 의 중앙", n=len(hr), rho=np.median(hr),
                 p=stats.wilcoxon(hr).pvalue))
P("   하플로타입 안 rho 중앙 %+.3f (하플로 %d · 윌콕슨 p %.3g)"
  % (np.median(hr), len(hr), stats.wilcoxon(hr).pvalue))
B = pd.DataFrame(rows)
B.to_csv(os.path.join(HERE, "FIGX19_B_NVAR.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="B 변이 개수↔메틸", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=PROX["sample"].nunique(), 하플로=len(hr), 유닛=len(PROX),
                 비고="근위 1-30 · 하플로 안 rho 는 유닛 15 개 이상인 하플로만"))

# ---------- C : 토큰 하나하나 (같은 하플로 · 같은 띠 안) ----------
P()
P("[C] 변이 토큰마다 — 같은 하플로타입·같은 위치 띠 안에서 '가진 유닛 vs 안 가진 유닛'")
ex = PROX[["hap", "띠", "ap", "f_hi", "vl"]].copy()
ex["toks"] = ex.vl.map(lambda v: [] if v == "-" else v.split(","))
long = ex.explode("toks").dropna(subset=["toks"])
cnt = long.toks.value_counts()
P("   근위 30 안에 나타난 토큰 %d 종 · 유닛 %s 개에 붙어 있다" % (len(cnt), format(len(long), ",")))
rows = []
for tok, n_unit in cnt.items():
    if n_unit < 200: continue
    has = long[long.toks == tok][["hap", "띠", "ap"]].assign(has=1).drop_duplicates(["hap", "띠", "ap"])
    m = ex.merge(has, on=["hap", "띠", "ap"], how="left")
    m["has"] = m.has.fillna(0)
    diffs = []
    for (h, b), g in m.groupby(["hap", "띠"]):
        a = g[g.has == 1].f_hi; c = g[g.has == 0].f_hi
        if len(a) >= 2 and len(c) >= 2: diffs.append(a.median() - c.median())
    if len(diffs) < 20: continue
    diffs = np.array(diffs)
    rows.append(dict(토큰=tok, 유닛=int(n_unit), 짝=len(diffs), Δ중앙=np.median(diffs),
                     내려간짝=int((diffs < 0).sum()), p=stats.wilcoxon(diffs).pvalue))
C = pd.DataFrame(rows).sort_values("p")
if len(C):
    C["q_BH"] = np.minimum.accumulate(
        (C.p.values * len(C) / np.arange(1, len(C) + 1))[::-1])[::-1]
    P("   검정한 토큰 %d 종 (근위 유닛 200 개 이상 · 짝 20 쌍 이상)" % len(C))
    P("   %-10s %7s %6s %9s %9s %10s %10s" % ("토큰", "유닛", "짝", "Δ중앙", "내려간짝", "p", "q(BH)"))
    for _, r in C.head(14).iterrows():
        P("   %-10s %7d %6d %+9.4f %6d/%-4d %10.3g %10.3g"
          % (r.토큰, r.유닛, r.짝, r.Δ중앙, r.내려간짝, r.짝, r.p, r.q_BH))
    P("   q<0.05 인 토큰 %d 종 / %d" % (int((C.q_BH < .05).sum()), len(C)))
C.to_csv(os.path.join(HERE, "FIGX19_C_TOKEN.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="C 토큰별 효과", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=PROX["sample"].nunique(), 하플로=PROX.hap.nunique(), 유닛=len(PROX),
                 비고="근위 1-30 · 토큰 %d 종 검정 · 짝=하플로×띠" % len(C)))

# ---------- D : 변이집합(유형) 단위 ----------
P()
P("[D] 변이집합(유형) 단위 — 같은 띠 안에서 유형마다 f_hi 가 다른가")
rows = []
for b, g in PROX.groupby("띠"):
    vc = g.vl.value_counts()
    top = [v for v in vc.index[:10] if vc[v] >= 30]
    groups = [g[g.vl == v].f_hi.values for v in top]
    if len(groups) < 3: continue
    kw = stats.kruskal(*groups)
    P("   %-12s 유형 %d 개(유닛 %d 이상) · Kruskal-Wallis p %.3g · f_hi 중앙 %.3f ~ %.3f"
      % (b, len(top), 30, kw.pvalue, min(np.median(x) for x in groups),
         max(np.median(x) for x in groups)))
    for v, x in zip(top, groups):
        rows.append(dict(띠=b, 유형=v[:60] if v != "-" else "변이 없음", 유닛=len(x),
                         하플로=g[g.vl == v].hap.nunique(), f_hi중앙=np.median(x),
                         f_hi평균=np.mean(x), KW_p=kw.pvalue, 변이수=0 if v == "-" else len(v.split(","))))
D = pd.DataFrame(rows)
D.to_csv(os.path.join(HERE, "FIGX19_D_TYPE.tsv"), sep="\t", index=False, float_format="%.4f")
PROV.append(dict(항목="D 변이집합 유형", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=PROX["sample"].nunique(), 하플로=PROX.hap.nunique(), 유닛=int(D.유닛.sum()),
                 비고="띠마다 상위 10 유형(유닛 30 이상)"))

# ---------- E : 같은 집합이면 메틸도 닮는가 (다른 하플로타입끼리) ----------
P()
P("[E] 같은 변이집합을 가진 유닛은 메틸도 닮는가 — 같은 띠 · 다른 하플로타입 유닛쌍")
# 자를 둘로 낸다 — 본표는 f_hi_thr (컷 0.508), 대조는 f_hi (컷 0.65).
#   까닭 : Fig4 묶음이 0.508 로 통일돼 있다 ( 요청).
def _epair(col, seed=19):
    _rng = np.random.default_rng(seed)
    rows = []
    for b, g in PROX.groupby("띠"):
        g2 = g[g.vl != "-"]
        same, diff = [], []
        vc = g2.vl.value_counts()
        for v in vc.index[:40]:
            s = g2[g2.vl == v]
            if len(s) < 4:
                continue
            idx = list(s.index)
            for i2, j2 in itertools.combinations(_rng.permutation(idx)[:12], 2):
                if s.at[i2, "hap"] != s.at[j2, "hap"]:
                    same.append(abs(s.at[i2, col] - s.at[j2, col]))
        pool = g2.sample(min(len(g2), 3000), random_state=19)
        ar = pool.sample(1200, replace=True, random_state=1)
        br = pool.sample(1200, replace=True, random_state=2)
        for (i2, r1), (j2, r2) in zip(ar.iterrows(), br.iterrows()):
            if r1.hap != r2.hap and r1.vl != r2.vl:
                diff.append(abs(r1[col] - r2[col]))
        if len(same) < 30 or len(diff) < 30:
            continue
        u = stats.mannwhitneyu(same, diff)
        rows.append(dict(띠=b, 같은집합쌍=len(same), 같은중앙=np.median(same),
                         다른집합쌍=len(diff), 다른중앙=np.median(diff),
                         Δ=np.median(same) - np.median(diff), p=u.pvalue))
    return pd.DataFrame(rows)


E = _epair("f_hi_thr")
E65 = _epair("f_hi")
for _, r in E.iterrows():
    P("   %-12s [컷 0.508] 같은 집합 쌍 %d (중앙 %.3f) vs 다른 집합 쌍 %d (중앙 %.3f) · Δ %+.3f · p %.3g"
      % (r.띠, r.같은집합쌍, r.같은중앙, r.다른집합쌍, r.다른중앙, r.Δ, r.p))
for _, r in E65.iterrows():
    P("   %-12s [컷 0.65 ] Δ %+.3f · p %.3g" % (r.띠, r.Δ, r.p))
E65.to_csv(os.path.join(HERE, "FIGX19_E_PAIR_cut065.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="E 같은 집합끼리 닮는가", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=PROX["sample"].nunique(), 하플로=PROX.hap.nunique(), 유닛=len(PROX),
                 비고="다른 하플로타입끼리만 짝지음 · 같은 위치 띠 · 컷 0.508(f_hi_thr) · 대조표 cut065"))
E.to_csv(os.path.join(HERE, "FIGX19_E_PAIR.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="E 같은 집합끼리 닮는가", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=PROX["sample"].nunique(), 하플로=PROX.hap.nunique(), 유닛=len(PROX),
                 비고="다른 하플로타입끼리만 짝지음 · 같은 위치 띠"))

# ---------- F : 변이가 어디에 있나 (구획) ----------
P()
P("[F] 변이가 유전자부에 있을 때 — 같은 띠 안에서 b_gene / b_sp")
rows = []
for b, g in PROX.groupby("띠"):
    g = g.dropna(subset=["b_gene_mean", "b_sp_mean"])
    a = g[g.n_var_gene > 0]; c = g[g.n_var_gene == 0]
    if len(a) < 30 or len(c) < 30: continue
    for col, nm in (("b_gene_mean", "유전자부 β"), ("b_sp_mean", "스페이서 β"), ("f_hi", "f_hi")):
        u = stats.mannwhitneyu(a[col], c[col])
        rows.append(dict(띠=b, 지표=nm, 유전자부변이있음=len(a), 없음=len(c),
                         중앙있음=a[col].median(), 중앙없음=c[col].median(),
                         Δ=a[col].median() - c[col].median(), p=u.pvalue))
        P("   %-12s %-8s 있음 %.3f (n=%d) vs 없음 %.3f (n=%d) · Δ %+.4f · p %.3g"
          % (b, nm, a[col].median(), len(a), c[col].median(), len(c),
             a[col].median() - c[col].median(), u.pvalue))
F = pd.DataFrame(rows)
#  FIGX19_F_SEG.tsv is written in English (column names, band and metric values). F itself
#    keeps its names, because the provenance row below still reads them. Reader : figS4_A_unit_beta_three_sources.py
F_EN = F.rename(columns={"띠": "band", "지표": "metric", "유전자부변이있음": "n_with_gene_body_variant",
                         "없음": "n_without", "중앙있음": "median_with", "중앙없음": "median_without"})
F_EN["band"] = F_EN["band"].map({"근위 1-3": "proximal 1-3", "근위 4-10": "proximal 4-10",
                                 "근위 11-30": "proximal 11-30"})
F_EN["metric"] = F_EN["metric"].map({"유전자부 β": "gene body β", "스페이서 β": "spacer β", "f_hi": "f_hi"})
F_EN.to_csv(os.path.join(HERE, "FIGX19_F_SEG.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="F 변이 위치(구획)", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=PROX["sample"].nunique(), 하플로=PROX.hap.nunique(),
                 유닛=int(F.유전자부변이있음.sum() + F.없음.sum()) if len(F) else 0,
                 비고="b_gene_mean·b_sp_mean 이 있는 유닛만"))

pd.DataFrame(PROV).to_csv(os.path.join(HERE, "FIGX19_PROVENANCE.tsv"), sep="\t", index=False)
P()
P("[자료 명세]")
P(pd.DataFrame(PROV).to_string(index=False))
P()
P("[저장] FIGX19_A_POS · B_NVAR · C_TOKEN · D_TYPE · E_PAIR · F_SEG · PROVENANCE")
LOG.close()
