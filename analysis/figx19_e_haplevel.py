# -*- coding: utf-8 -*-
# ① 을 하플로타입 하나에 값 하나 로 — 하위표집 없음 · 난수 없음 · 모든 교차쌍 정확 집계
#
#   하플로타입 h 마다 두 값을 낸다 (같은 띠 안, 다른 하플로타입의 유닛과만 견줌) :
#     같은집합거리(h) = h 의 유닛 ↔ 같은 변이집합 을 가진 다른 하플로타입 유닛 의 평균 |Δ|
#     다른집합거리(h) = h 의 유닛 ↔ 다른 변이집합 을 가진 다른 하플로타입 유닛 의 평균 |Δ|
#   → 하플로타입당 값 둘 → 짝지은 윌콕슨 부호순위 검정 (관측이 독립)
#
#   합 |a-b| 는 정렬 + 누적합으로 전수 계산한다 (표집 없음 ⇒ 난수 의존 없음).
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats


def sumabs_within(values, group):
    """같은 group 안에서 각 원소가 그 group 의 모든 원소와 갖는 |차| 의 합, 그리고 group 크기."""
    v = np.asarray(values, float)
    g = np.asarray(group)
    out = np.zeros(len(v)); cnt = np.zeros(len(v))
    order = np.lexsort((v, g))
    gs = g[order]; vs = v[order]
    bounds = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1], True])
    for b0, b1 in zip(bounds[:-1], bounds[1:]):
        a = vs[b0:b1]; m = b1 - b0
        P = np.concatenate(([0.0], np.cumsum(a)))
        r = np.arange(m)
        s = r * a - P[r] + (P[m] - P[r + 1]) - (m - 1 - r) * a
        out[order[b0:b1]] = s; cnt[order[b0:b1]] = m - 1
    return out, cnt


def run(col):
    L = pd.read_parquet(__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/02_tables/L2_HAP_UNIT.parquet",
                        columns=["sample", "slot", "ap", "ad", "n_read", "varlist", "f_hi", "f_hi_thr"])
    L = L[L.n_read >= 5].copy()
    L["hap"] = L["sample"].astype(str) + "|" + L.slot.astype(str)
    L["vl"] = L.varlist.fillna("-").astype(str)
    # 공통 5 구역 (2026-10-16) = 정본 8 구역을 묶은 것 Z1 / Z2+Z3 / ZI / Z4 / Z5+Z6+Z7.
    #   main (F)(G)(H) · FigX26 · supp (G) 와 같은 조건식 이다.
    def band(ap, ad):
        if ap <= ad:
            if ap <= 3:  return "근위 1-3"
            if ap <= 17: return "근위 4-17"
            return "내부"
        if ad <= 10: return "원위 1-10"
        if ad <= 22: return "원위 11-22"
        return "내부"
    L["띠"] = [band(a, d) for a, d in zip(L.ap.values, L.ad.values)]
    rows = []
    for b in ("근위 1-3", "근위 4-17", "내부", "원위 11-22", "원위 1-10"):
        g = L[(L.띠 == b) & (L.vl != "-")].reset_index(drop=True)
        v = g[col].values.astype(float)
        hap = pd.factorize(g.hap)[0]
        vl = pd.factorize(g.vl)[0]
        hv = hap.astype(np.int64) * (vl.max() + 1) + vl          # (hap, set) 합성 키
        one = np.zeros(len(v), dtype=np.int64)
        S_all, C_all = sumabs_within(v, one)                      # 전체
        S_hap, C_hap = sumabs_within(v, hap)                      # 같은 하플로 안
        S_set, C_set = sumabs_within(v, vl)                       # 같은 집합 안
        S_hs,  C_hs  = sumabs_within(v, hv)                       # 같은 하플로 & 같은 집합
        S_same = S_set - S_hs;  C_same = C_set - C_hs             # 같은 집합 · 다른 하플로
        S_cross = S_all - S_hap; C_cross = C_all - C_hap          # 다른 하플로 전부
        S_diff = S_cross - S_same; C_diff = C_cross - C_same      # 다른 집합 · 다른 하플로
        d = pd.DataFrame(dict(hap=g.hap.values, Ss=S_same, Cs=C_same, Sd=S_diff, Cd=C_diff))
        a = d.groupby("hap").sum()
        a = a[(a.Cs > 0) & (a.Cd > 0)]
        same = a.Ss / a.Cs; diff = a.Sd / a.Cd
        w = stats.wilcoxon(same, diff)
        rows.append(dict(띠=b, 하플로타입=len(a), 유닛=len(g),
                         같은집합거리중앙=float(np.median(same)),
                         다른집합거리중앙=float(np.median(diff)),
                         Δ=float(np.median(same - diff)),
                         더가까운하플로=int((same < diff).sum()),
                         p=float(w.pvalue)))
    return pd.DataFrame(rows)


for col, lab in (("f_hi_thr", "컷 0.508 (f_hi_thr) 본표"), ("f_hi", "컷 0.65 (f_hi) 대조")):
    R = run(col)
    print("\n%s" % lab)
    for _, r in R.iterrows():
        print("   %-12s 하플로타입 %3d (유닛 %5d) · 같은집합 %.4f ↔ 다른집합 %.4f · "
              "짝 안 Δ 중앙 %+.4f · 더 가까운 하플로 %d/%d · 윌콕슨 짝 p %.4g"
              % (r.띠, r.하플로타입, r.유닛, r.같은집합거리중앙, r.다른집합거리중앙,
                 r.Δ, r.더가까운하플로, r.하플로타입, r.p))
    R.to_csv(__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/03_figures/paper_Fig/FIGX19_E_HAP_%s.tsv"
             % ("cut0508" if col == "f_hi_thr" else "cut065"),
             sep="\t", index=False, float_format="%.6g")
print("\n[저장] FIGX19_E_HAP_cut0508.tsv · FIGX19_E_HAP_cut065.tsv")
