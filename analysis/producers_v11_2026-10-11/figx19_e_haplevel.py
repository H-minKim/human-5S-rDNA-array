# -*- coding: utf-8 -*-
# Reduce 1) to one value per haplotype: no subsampling, no random numbers, every cross pair counted exactly.
#
#   Two values are produced per haplotype h, comparing only against units of other haplotypes within the same position band:
#     same-set distance(h) = the mean |delta| between h's units and units of other haplotypes carrying the same variant set
#     other-set distance(h) = the mean |delta| between h's units and units of other haplotypes carrying a different variant set
#   -> two values per haplotype -> a paired Wilcoxon signed-rank test, so the observations are independent.
#
#   The sum of |a-b| is computed exhaustively by sorting and a cumulative sum, so there is no sampling and no dependence on random numbers.
import os          # os.path.join is used below, so os must be imported here.
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    The roots of the data this script reads are listed in the
    "Running this code" section of README.md."""
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



def sumabs_within(values, group):
    """within one group, for each element the sum of |difference| to every element of that group, and the group size."""
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
    L = pd.read_parquet(os.path.join(_need("RDNA_V11_ROOT"), "02_tables", "L2_HAP_UNIT.parquet"),
                        columns=["sample", "slot", "ap", "n_read", "varlist", "f_hi", "f_hi_thr"])
    L = L[L.n_read >= 5].copy()
    L["hap"] = L["sample"].astype(str) + "|" + L.slot.astype(str)
    L["vl"] = L.varlist.fillna("-").astype(str)
    def band(ap):
        return "근위 1-3" if ap <= 3 else ("근위 4-10" if ap <= 10 else ("근위 11-30" if ap <= 30 else "그 밖"))
    L["띠"] = [band(a) for a in L.ap]
    rows = []
    for b in ("근위 1-3", "근위 4-10", "근위 11-30"):
        g = L[(L.띠 == b) & (L.vl != "-")].reset_index(drop=True)
        v = g[col].values.astype(float)
        hap = pd.factorize(g.hap)[0]
        vl = pd.factorize(g.vl)[0]
        hv = hap.astype(np.int64) * (vl.max() + 1) + vl          # composite (hap, set) key
        one = np.zeros(len(v), dtype=np.int64)
        S_all, C_all = sumabs_within(v, one)                      # all
        S_hap, C_hap = sumabs_within(v, hap)                      # within the same haplotype
        S_set, C_set = sumabs_within(v, vl)                       # within the same variant set
        S_hs,  C_hs  = sumabs_within(v, hv)                       # same haplotype and same variant set
        S_same = S_set - S_hs;  C_same = C_set - C_hs             # same variant set, different haplotype
        S_cross = S_all - S_hap; C_cross = C_all - C_hap          # all different haplotypes
        S_diff = S_cross - S_same; C_diff = C_cross - C_same      # different variant set, different haplotype
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
    #  the table is written in English (column names and the band values). Only the written
    #    copy is renamed. Reader : figS4_A_unit_beta_three_sources.py
    R_EN = R.rename(columns={"띠": "band", "하플로타입": "haplotypes", "유닛": "units",
                             "같은집합거리중앙": "same_set_distance_median",
                             "다른집합거리중앙": "different_set_distance_median",
                             "더가까운하플로": "haplotypes_closer_to_same_set"})
    R_EN["band"] = R_EN["band"].map({"근위 1-3": "proximal 1-3", "근위 4-10": "proximal 4-10",
                                     "근위 11-30": "proximal 11-30"})
    R_EN.to_csv(os.path.join(_need("RDNA_V11_ROOT"), "03_figures", "paper_Fig", "FIGX19_E_HAP_%s.tsv")
                % ("cut0508" if col == "f_hi_thr" else "cut065"),
                sep="\t", index=False, float_format="%.6g")
print("\n[저장] FIGX19_E_HAP_cut0508.tsv · FIGX19_E_HAP_cut065.tsv")
