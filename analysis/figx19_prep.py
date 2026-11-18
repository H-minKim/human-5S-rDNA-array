# -*- coding: utf-8 -*-
# FigX19 data: the effect of the proximal variant set on unit methylation, per unit.
# A script of the same name is also in the repository:
#   analysis/producers_v11_2026-10-11/figx19_prep.py.
#   The name is the same but the table differs: FIGX19_F_SEG.tsv from this script has 15 rows and 9 columns including n_hap,
#   while the other has 9 rows and 8 columns without n_hap. They are told apart by row and column count, not by md5.
#   Figure 3 and Supplementary Figure 3 read the table this script writes.
#   Data: L2_HAP_UNIT.parquet (haplotype x unit, HPRC_ONT_submissions).
#   Convention: ap is the position counted from the proximal end (ap 1 = the proximal terminus) and ad is counted from the distal end; both confirmed by measurement.
#   Position dominates methylation, so every comparison is made within a position band, and the individual offset is removed within haplotype.
import io, os, sys, itertools
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import itertools
import numpy as np, pandas as pd
from scipy import stats
# The data roots come from environment variables.
#   They have no default: the script stops with a message if one is not set.
DATA_ROOT = os.environ["DATA_ROOT"]
# Where the tables are read and written. The tables are not deposited next to the code, so these paths come from the environment.
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
LOG = io.open(os.path.join(HERE, "FIGX19_PREP.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
rng = np.random.default_rng(19)
T = DATA_ROOT + "/02_tables/"
MINREAD = 5

L = pd.read_parquet(T + "L2_HAP_UNIT.parquet", columns=[
    "sample", "slot", "ap", "ad", "N", "side", "f_hi", "n_read", "varlist", "n_var",
    "n_var_common", "n_var_rare", "n_var_gene", "nv_coding", "vsig", "vsig_n",
    "b_gene_mean", "b_sp_mean", "b_mean", "data_source", "pop", "super_pop",
    "f_hi_thr"])
P("[input] L2_HAP_UNIT.parquet, %s rows, data source %s" % (format(len(L), ","), set(L.data_source)))
L["hap"] = L["sample"].astype(str) + "|" + L.slot.astype(str)
P("   individuals %d, haplotypes %d, units %s"
  % (L["sample"].nunique(), L.hap.nunique(), format(len(L), ",")))
L = L[L.n_read >= MINREAD].copy()
L["vl"] = L.varlist.fillna("-").astype(str)
P("   read >= %d only: units %s, individuals %d, haplotypes %d"
  % (MINREAD, format(len(L), ","), L["sample"].nunique(), L.hap.nunique()))
PROV = [dict(항목="바탕 자료", 자료="L2_HAP_UNIT.parquet", 자료원="HPRC_ONT_submissions",
             개체=L["sample"].nunique(), 하플로=L.hap.nunique(), 유닛=len(L),
             비고="read>=%d 인 유닛만" % MINREAD)]

# ---------- A: position comes first ----------
P()
P("[A] position dominates methylation: f_hi by proximal position (so every comparison below is made within a position band)")
A = (L[L.ap <= 20].groupby("ap")
     .agg(유닛=("f_hi", "size"), 하플로=("hap", "nunique"), f_hi중앙=("f_hi", "median"),
          f_hi평균=("f_hi", "mean"), 변이수평균=("n_var", "mean")).reset_index())
P(A.to_string(index=False))
A.to_csv(os.path.join(HERE, "FIGX19_A_POS.tsv"), sep="\t", index=False, float_format="%.4f")

# Positional zones: the common five.
#   proximal 1-3, proximal 4-17, interior, distal 11-22, distal 1-10
#   They group the canonical eight zones (zone_state_seg.py:37): Z1 / Z2+Z3 / ZI / Z4 / Z5+Z6+Z7.
#   No new boundary is introduced. In a short array one unit can fall at both ends, so the split is made on ap <= ad.
def band(ap, ad):
    if ap <= ad:
        if ap <= 3:  return "근위 1-3"
        if ap <= 17: return "근위 4-17"
        return "내부"
    if ad <= 10: return "원위 1-10"
    if ad <= 22: return "원위 11-22"
    return "내부"
ZORD = ["근위 1-3", "근위 4-17", "내부", "원위 11-22", "원위 1-10"]
L["띠"] = [band(a, d) for a, d in zip(L.ap.values, L.ad.values)]
assert L.띠.isin(ZORD).all(), "구역 배정에서 빠진 유닛이 있다"
ZU = L.copy()                      # the zone units, which cover more than the proximal end
# Also reported: the number of haplotypes with an observed unit in each band, so the denominator can be read straight from the table.
#   This is a count within the same set.
NHAP = ZU.groupby("띠").hap.nunique().to_dict()
P("   haplotypes by zone: " + " · ".join("%s %d" % (k, NHAP[k]) for k in sorted(NHAP)))
P()
P("   all 5 zones: units %s, haplotypes %d, individuals %d"
  % (format(len(ZU), ","), ZU.hap.nunique(), ZU["sample"].nunique()))

# ---------- B: variant count against methylation, before and after controlling for position ----------
P()
P("[B] variant count and f_hi: before/after control")
rows = []
r0, p0 = stats.spearmanr(ZU.n_var, ZU.f_hi)
rows.append(dict(구분="통제 전 (5 구역 전부)", n=len(ZU), rho=r0, p=p0))
P("   before control: rho %+.3f (p %.3g, units %d)" % (r0, p0, len(ZU)))
for b, g in ZU.groupby("띠"):
    r, p = stats.spearmanr(g.n_var, g.f_hi)
    rows.append(dict(구분="위치 띠 안 : " + b, n=len(g), rho=r, p=p))
    P("   %-12s rho %+.3f (p %.3g, units %d, haplotypes %d)" % (b, r, p, len(g), g.hap.nunique()))
# Within haplotype, removing the individual offset: rho is computed per haplotype and a sign test applied.
hr = []
for h, g in ZU.groupby("hap"):
    if g.n_var.nunique() > 2 and len(g) >= 15:
        hr.append(stats.spearmanr(g.n_var, g.f_hi)[0])
hr = np.array([x for x in hr if x == x])
rows.append(dict(구분="하플로타입 안 rho 의 중앙", n=len(hr), rho=np.median(hr),
                 p=stats.wilcoxon(hr).pvalue))
P("   within-haplotype rho median %+.3f (haplotypes %d, Wilcoxon p %.3g)"
  % (np.median(hr), len(hr), stats.wilcoxon(hr).pvalue))
B = pd.DataFrame(rows)
B.to_csv(os.path.join(HERE, "FIGX19_B_NVAR.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="B 변이 개수↔메틸", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=ZU["sample"].nunique(), 하플로=len(hr), 유닛=len(ZU),
                 비고="5 구역 전부 · 하플로 안 rho 는 유닛 15 개 이상인 하플로만"))

# ---------- C: token by token, within the same haplotype and the same band ----------
P()
P("[C] per variant token: within the same haplotype and the same position band, 'units with it vs units without it'")
ex = ZU[["hap", "띠", "ap", "f_hi", "vl"]].copy()
ex["toks"] = ex.vl.map(lambda v: [] if v == "-" else v.split(","))
long = ex.explode("toks").dropna(subset=["toks"])
cnt = long.toks.value_counts()
P("   tokens seen within proximal 30: %d kinds, carried by %s units" % (len(cnt), format(len(long), ",")))
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
    P("   tokens tested: %d kinds (200 or more proximal units, 20 or more pairs)" % len(C))
    P("   %-10s %7s %6s %9s %9s %10s %10s" % ("토큰", "유닛", "짝", "Δ중앙", "내려간짝", "p", "q(BH)"))
    for _, r in C.head(14).iterrows():
        P("   %-10s %7d %6d %+9.4f %6d/%-4d %10.3g %10.3g"
          % (r.토큰, r.유닛, r.짝, r.Δ중앙, r.내려간짝, r.짝, r.p, r.q_BH))
    P("   tokens with q<0.05: %d kinds / %d" % (int((C.q_BH < .05).sum()), len(C)))
C.to_csv(os.path.join(HERE, "FIGX19_C_TOKEN.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="C 토큰별 효과", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=ZU["sample"].nunique(), 하플로=ZU.hap.nunique(), 유닛=len(ZU),
                 비고="5 구역 · 토큰 %d 종 검정 · 짝=하플로×구역" % len(C)))

# ---------- D: per variant set (haplotype type) ----------
P()
P("[D] by variant set (type): does f_hi differ between types within the same band")
rows = []
for b, g in ZU.groupby("띠"):
    vc = g.vl.value_counts()
    top = [v for v in vc.index[:10] if vc[v] >= 30]
    groups = [g[g.vl == v].f_hi.values for v in top]
    if len(groups) < 3: continue
    kw = stats.kruskal(*groups)
    P("   %-12s types %d (%d or more units), Kruskal-Wallis p %.3g, f_hi median %.3f ~ %.3f"
      % (b, len(top), 30, kw.pvalue, min(np.median(x) for x in groups),
         max(np.median(x) for x in groups)))
    for v, x in zip(top, groups):
        rows.append(dict(띠=b, 유형=v[:60] if v != "-" else "변이 없음", 유닛=len(x),
                         하플로=g[g.vl == v].hap.nunique(), f_hi중앙=np.median(x),
                         f_hi평균=np.mean(x), KW_p=kw.pvalue, 변이수=0 if v == "-" else len(v.split(","))))
D = pd.DataFrame(rows)
D.to_csv(os.path.join(HERE, "FIGX19_D_TYPE.tsv"), sep="\t", index=False, float_format="%.4f")
PROV.append(dict(항목="D 변이집합 유형", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=ZU["sample"].nunique(), 하플로=ZU.hap.nunique(), 유닛=int(D.유닛.sum()),
                 비고="구역마다 상위 10 유형(유닛 30 이상)"))

# ---------- E: do units with the same variant set also resemble each other in methylation, across haplotypes ----------
P()
P("[E] do units with the same variant set also resemble each other in methylation: same band, unit pairs from different haplotypes")
# Two rulers are reported: the main table uses f_hi_thr (cut 0.508) and the control uses f_hi (cut 0.65).
#   The Figure 4 set uses 0.508 throughout.
def _epair(col, seed=19):
    _rng = np.random.default_rng(seed)
    rows = []
    for b, g in ZU.groupby("띠"):
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
        rows.append(dict(띠=b, n_hap=NHAP.get(b, 0),
                         같은집합쌍=len(same), 같은중앙=np.median(same),
                         다른집합쌍=len(diff), 다른중앙=np.median(diff),
                         Δ=np.median(same) - np.median(diff), p=u.pvalue))
    return pd.DataFrame(rows)


E = _epair("f_hi_thr")
E65 = _epair("f_hi")
for _, r in E.iterrows():
    P("   %-12s [cut 0.508] same-set pairs %d (median %.3f) vs different-set pairs %d (median %.3f), Δ %+.3f, p %.3g"
      % (r.띠, r.같은집합쌍, r.같은중앙, r.다른집합쌍, r.다른중앙, r.Δ, r.p))
for _, r in E65.iterrows():
    P("   %-12s [cut 0.65 ] Δ %+.3f, p %.3g" % (r.띠, r.Δ, r.p))
E65.to_csv(os.path.join(HERE, "FIGX19_E_PAIR_cut065.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="E 같은 집합끼리 닮는가", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=ZU["sample"].nunique(), 하플로=ZU.hap.nunique(), 유닛=len(ZU),
                 비고="다른 하플로타입끼리만 짝지음 · 같은 위치 띠 · 컷 0.508(f_hi_thr) · 대조표 cut065"))
E.to_csv(os.path.join(HERE, "FIGX19_E_PAIR.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="E 같은 집합끼리 닮는가", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=ZU["sample"].nunique(), 하플로=ZU.hap.nunique(), 유닛=len(ZU),
                 비고="다른 하플로타입끼리만 짝지음 · 같은 위치 띠"))

# ---------- F: where the variants lie (compartment) ----------
P()
P("[F] when the variant is in the gene region: b_gene / b_sp within the same band")
rows = []
for b, g in ZU.groupby("띠"):
    g = g.dropna(subset=["b_gene_mean", "b_sp_mean"])
    a = g[g.n_var_gene > 0]; c = g[g.n_var_gene == 0]
    if len(a) < 30 or len(c) < 30: continue
    for col, nm in (("b_gene_mean", "유전자부 β"), ("b_sp_mean", "스페이서 β"), ("f_hi", "f_hi")):
        u = stats.mannwhitneyu(a[col], c[col])
        rows.append(dict(띠=b, 지표=nm, n_hap=NHAP.get(b, 0),
                         유전자부변이있음=len(a), 없음=len(c),
                         중앙있음=a[col].median(), 중앙없음=c[col].median(),
                         Δ=a[col].median() - c[col].median(), p=u.pvalue))
        P("   %-12s %-8s present %.3f (n=%d) vs none %.3f (n=%d), Δ %+.4f, p %.3g"
          % (b, nm, a[col].median(), len(a), c[col].median(), len(c),
             a[col].median() - c[col].median(), u.pvalue))
F = pd.DataFrame(rows)
F.to_csv(os.path.join(HERE, "FIGX19_F_SEG.tsv"), sep="\t", index=False, float_format="%.6g")
PROV.append(dict(항목="F 변이 위치(구획)", 자료="L2_HAP_UNIT", 자료원="HPRC_ONT_submissions",
                 개체=ZU["sample"].nunique(), 하플로=ZU.hap.nunique(),
                 # The three measures emit the same unit on three rows, so summing them would treble the count.
                 # (44,835 against the correct 14,945.)
                 유닛=int(F[F.지표 == "유전자부 β"].유전자부변이있음.sum()
                        + F[F.지표 == "유전자부 β"].없음.sum()) if len(F) else 0,
                 비고="b_gene_mean·b_sp_mean 이 있는 유닛만"))

pd.DataFrame(PROV).to_csv(os.path.join(HERE, "FIGX19_PROVENANCE.tsv"), sep="\t", index=False)
P()
P("[data specification]")
P(pd.DataFrame(PROV).to_string(index=False))
P()
P("[saved] FIGX19_A_POS, B_NVAR, C_TOKEN, D_TYPE, E_PAIR, F_SEG, PROVENANCE")
LOG.close()
