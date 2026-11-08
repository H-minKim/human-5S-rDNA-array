# -*- coding: utf-8 -*-
# true n-by-n exhaustive scan, not only methylation but every item x every item
#
#   The question: "which item of a haplotype is related to which item at which position
#              of that haplotype", e.g. array length N and a specific variant in proximal units 4-7
#                                  e.g. a specific variant rate along positions and methylation
#
#   [A] build the haplotype feature matrix (one row = one slot)
#        columns = (position window x variable) + (position window x specific variant carrier fraction) + haplotype constants
#   [B] exhaustive pairwise Spearman + BH-FDR.  n = 510 slots.
#        two slots of the same individual are included, so reproducibility is checked with one slot per individual (n=280).
#   [C] trivial pairs (neighboring windows of the same variable and the like) are flagged and filtered out.
#   [D] unit level: the methylation difference between units with and without a specific variant
#        is compared within the same slot and the same position window (removing position effects) plus a within-stratum permutation null.
#
#  matplotlib glyph trap: ASCII only.
import io, os, warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"
T = OUT + "/02_tables"
LOG = io.open(OUT + "/04_logs/V20_NXN.txt", "w", encoding="utf8", newline="\n")
SEED = 20260914
rng = np.random.RandomState(SEED)


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


P("=" * 112)
P("n 대 n 전수 주사 — 하플로타입의 모든 항목 x 모든 항목")
P("=" * 112)
P()

D = pd.read_parquet(os.path.join(T, "L2_HAP_UNIT.parquet"))
D = D[D.n_read >= 5].copy()
D["sl"] = D["sample"].astype(str) + "_" + D.slot.astype(str)
P("  L2 자리 %s · 슬롯 %d · 개체 %d"
  % (format(len(D), ","), D.sl.nunique(), D["sample"].nunique()))
P()

# ---------------------------------------------------------------- [A] position windows
WIN = [("근1", (D.ap == 1)), ("근2", (D.ap == 2)), ("근3", (D.ap == 3)),
       ("근4_7", (D.ap >= 4) & (D.ap <= 7)),
       ("근8_12", (D.ap >= 8) & (D.ap <= 12)),
       ("근13_17", (D.ap >= 13) & (D.ap <= 17)),
       ("근18_25", (D.ap >= 18) & (D.ap <= 25)),
       ("중앙", (D.ap > 25) & (D.ad > 25)),
       ("원18_25", (D.ad >= 18) & (D.ad <= 25)),
       ("원13_17", (D.ad >= 13) & (D.ad <= 17)),
       ("원8_12", (D.ad >= 8) & (D.ad <= 12)),
       ("원4_7", (D.ad >= 4) & (D.ad <= 7)),
       ("원3", (D.ad == 3)), ("원2", (D.ad == 2)), ("원1", (D.ad == 1))]

VARS = ["f_hi", "b_lo_mean", "b_hi_mean", "b_var_true", "p_switch",
        "b_gene_mean", "b_sp_mean", "r_lo_gene_sp", "r_hi_gene_sp",
        "blk_len_mean", "agree_next",
        "n_var_rare", "n_var_cpg", "n_var_coding", "n_ins", "n_del",
        "nv_AluY", "nv_NTS1", "nv_GArich", "nv_NTS2", "nv_coding",
        "unit_len", "vsig_minfreq"]
VARS = [v for v in VARS if v in D.columns]

slots = sorted(D.sl.unique())
SIDX = {s: i for i, s in enumerate(slots)}
F = pd.DataFrame(index=slots)
MINU = 2
for wn, m in WIN:
    s = D[m]
    g = s.groupby("sl", observed=True)
    cnt = g.size()
    keep = cnt[cnt >= MINU].index
    for v in VARS:
        mu = g[v].mean().reindex(slots)
        mu[~mu.index.isin(keep)] = np.nan
        F["%s|%s" % (wn, v)] = mu.values
P("  [A] 자리창 %d x 변수 %d = %d 열" % (len(WIN), len(VARS), len(WIN) * len(VARS)))

# ---------------------------------------------------------------- specific variant carrier fraction
VW = [("근1_3", (D.ap <= 3)), ("근4_7", (D.ap >= 4) & (D.ap <= 7)),
      ("근8_12", (D.ap >= 8) & (D.ap <= 12)),
      ("근13_25", (D.ap >= 13) & (D.ap <= 25)),
      ("중앙", (D.ap > 25) & (D.ad > 25)),
      ("원13_25", (D.ad >= 13) & (D.ad <= 25)),
      ("원8_12", (D.ad >= 8) & (D.ad <= 12)),
      ("원4_7", (D.ad >= 4) & (D.ad <= 7)),
      ("원1_3", (D.ad <= 3))]

tok_lists = D.varlist.fillna("-").astype(str).str.split(",")
from collections import Counter
cnt_all = Counter()
for tl in tok_lists:
    for t in tl:
        if t and t != "-":
            cnt_all[t] += 1
cand = [t for t, c in cnt_all.items() if 0.02 * len(D) <= c <= 0.98 * len(D)]
P("  후보 변이 토큰(유닛 보유율 2-98%%) %d 개" % len(cand))

TOKM = pd.DataFrame({t: tok_lists.map(lambda L, t=t: 1.0 if t in L else 0.0).values
                     for t in cand}, index=D.index)
# only tokens with large between-haplotype variation
slot_mean = TOKM.groupby(D.sl.values).mean()
sd = slot_mean.std()
TOP = sd.sort_values(ascending=False).head(40).index.tolist()
P("  하플로타입 사이 변동이 가장 큰 변이 토큰 40 개를 쓴다")
P("    " + " ".join(TOP[:20]))
P("    " + " ".join(TOP[20:]))
for wn, m in VW:
    sub = TOKM[m.values]
    key = D.sl.values[m.values]
    g = sub.groupby(key)
    cntw = g.size()
    keep = cntw[cntw >= MINU].index
    mu = g.mean().reindex(slots)
    mu.loc[~mu.index.isin(keep), :] = np.nan
    for t in TOP:
        F["%s|VAR_%s" % (wn, t)] = mu[t].values
P("  변이 창 %d x 토큰 %d = %d 열" % (len(VW), len(TOP), len(VW) * len(TOP)))

# ---------------------------------------------------------------- haplotype constants
HC = (D.drop_duplicates("sl").set_index("sl")
      [["N", "f_hap", "Bp", "Bd", "n_low", "depth", "gate_rate", "cov_frac",
        "n_read_slot", "N_other", "dN", "n_break"]].reindex(slots))
for c in HC.columns:
    F["HAP|" + c] = HC[c].values
P("  하플로타입 상수 %d 열" % HC.shape[1])

ok = F.notna().sum()
F = F.loc[:, ok >= 300]
F = F.loc[:, F.std() > 0]
P("  관측 300 슬롯 이상 · 분산>0 인 열만 남김 : %d 열" % F.shape[1])
F.to_parquet(os.path.join(T, "HAPFEAT.parquet"))
P("  저장 : 02_tables/HAPFEAT.parquet  (%d 슬롯 x %d 열)" % F.shape)
P()

# ---------------------------------------------------------------- [B] exhaustive pairs
def pair_spearman(Fx):
    """Rank each column within its observed values, then compute pairwise-complete Pearson for every pair with matrix products."""
    cols = list(Fx.columns)
    V = Fx.values.astype(float)
    M = np.isfinite(V).astype(float)
    Z = np.zeros_like(V)
    for j in range(V.shape[1]):
        m = np.isfinite(V[:, j])
        Z[m, j] = stats.rankdata(V[m, j])
    X = Z * M
    n = M.T @ M
    Sx = X.T @ M            # [j,k] = sum_i x_j (i observed in both j and k)
    Sy = Sx.T
    Sxy = X.T @ X
    Sxx = (X * X).T @ M
    Syy = Sxx.T
    num = n * Sxy - Sx * Sy
    den = np.sqrt(np.maximum(n * Sxx - Sx ** 2, 1e-12) *
                  np.maximum(n * Syy - Sy ** 2, 1e-12))
    r = num / den
    np.fill_diagonal(r, 1.0)
    return cols, r, n


P("=" * 112)
P("[B] 전수 쌍 Spearman — %d 열 -> %s 쌍"
  % (F.shape[1], format(F.shape[1] * (F.shape[1] - 1) // 2, ",")))
P("=" * 112)
P()
P("  열마다 관측값 안에서 한 번 순위화한 뒤 쌍별 결측 제외 Pearson 을 쓴다")
P("    (쌍마다 다시 순위를 매기는 것과 미세하게 다르지만 표준 방식이다)")
P()
cols, R, NN = pair_spearman(F)
iu = np.triu_indices(len(cols), 1)
rv = R[iu]
nv = NN[iu]
good = np.isfinite(rv) & (nv >= 200)
tt = rv * np.sqrt(np.maximum(nv - 2, 1) / np.maximum(1 - rv ** 2, 1e-12))
pv = 2 * stats.t.sf(np.abs(tt), np.maximum(nv - 2, 1))
pv = np.where(good, pv, 1.0)
o = np.argsort(pv)
q = np.empty_like(pv)
m_ = good.sum()
ranked = pv[o]
qv = ranked * m_ / (np.arange(len(ranked)) + 1)
qv = np.minimum.accumulate(qv[::-1])[::-1]
q[o] = np.clip(qv, 0, 1)

A = np.array(cols)[iu[0]]
B = np.array(cols)[iu[1]]


def fam(c):
    w, v = c.split("|", 1)
    return w, ("VAR" if v.startswith("VAR_") else v)


wa, va = zip(*[fam(c) for c in A])
wb, vb = zip(*[fam(c) for c in B])
PR = pd.DataFrame({"A": A, "B": B, "winA": wa, "varA": va, "winB": wb, "varB": vb,
                   "r": rv, "n": nv, "p": pv, "q": q})
PR = PR[good]
PR["same_var"] = (PR.varA == PR.varB)
PR["same_win"] = (PR.winA == PR.winB)
P("  검정한 쌍 %s · FDR 5%% 통과 %s (%.3f)"
  % (format(len(PR), ","), format(int((PR.q < 0.05).sum()), ","),
     (PR.q < 0.05).mean()))
P("  그중 같은 변수의 다른 창(뻔함) %s · 같은 창 안(뻔함) %s"
  % (format(int(((PR.q < 0.05) & PR.same_var).sum()), ","),
     format(int(((PR.q < 0.05) & PR.same_win).sum()), ",")))
P()

# ---------------------------------------------------------------- reproduction with 1 slot per individual
one = (D.drop_duplicates("sl")[["sample", "sl"]]
       .groupby("sample", observed=True).first()["sl"].values)
F1 = F.loc[[s for s in one if s in F.index]]
_, R1, N1 = pair_spearman(F1)
r1 = R1[iu]
PR["r_1slot"] = r1[good]
P("  개체당 한 슬롯만(n=%d) 으로 다시 재면 r 이 얼마나 유지되나 : "
  "전체 쌍 상관 %.4f" % (len(F1), np.corrcoef(PR.r[np.isfinite(PR.r_1slot)],
                                              PR.r_1slot[np.isfinite(PR.r_1slot)])[0, 1]))
P()

PR.sort_values("q").to_csv(os.path.join(T, "NXN_PAIRS.tsv"), sep="\t", index=False,
                           float_format="%.5g")
P("  저장 : 02_tables/NXN_PAIRS.tsv (모든 쌍)")
P()


def show(title, sel, k=14, note=""):
    P("-" * 112)
    P(title)
    if note:
        P("   " + note)
    P("-" * 112)
    s = sel.reindex(sel.r.abs().sort_values(ascending=False).index)
    s = s[s.q < 0.05].head(k)
    if not len(s):
        P("   (FDR 5% 를 넘는 쌍 없음)")
        P()
        return
    P("   %-26s %-26s %7s %6s %9s %8s" % ("항목 A", "항목 B", "r", "n", "q", "1슬롯r"))
    for t in s.itertuples(index=False):
        P("   %-26s %-26s %+.3f %6d %9.2e %+8.3f"
          % (t.A, t.B, t.r, t.n, t.q, t.r_1slot))
    P()


NONMETH = {"n_var_rare", "n_var_cpg", "n_var_coding", "n_ins", "n_del", "nv_AluY",
           "nv_NTS1", "nv_GArich", "nv_NTS2", "nv_coding", "unit_len",
           "vsig_minfreq", "VAR"}
METH = {"f_hi", "b_lo_mean", "b_hi_mean", "b_var_true", "p_switch", "b_gene_mean",
        "b_sp_mean", "r_lo_gene_sp", "r_hi_gene_sp", "blk_len_mean", "agree_next"}
TECH = {"depth", "gate_rate", "cov_frac", "n_read_slot", "n_break"}

P("=" * 112)
P("[C] 무엇이 무엇과 붙어 있나 — 뻔한 쌍을 걷어내고 갈래별로")
P("=" * 112)
P()

isHAP = (PR.winA == "HAP") | (PR.winB == "HAP")
show("(1) array 길이 N x 특정 자리의 특정 변이  (물어보신 바로 그 조합)",
     PR[((PR.varA == "N") & (PR.varB == "VAR")) |
        ((PR.varB == "N") & (PR.varA == "VAR"))],
     note="한쪽은 하플로타입 길이, 다른 쪽은 그 창 유닛들의 그 변이 보유율")

show("(2) 특정 변이 x 같은 창의 메틸",
     PR[(~PR.same_win) & (False) |
        (PR.same_win & (((PR.varA == "VAR") & PR.varB.isin(METH)) |
                        ((PR.varB == "VAR") & PR.varA.isin(METH))))],
     note="같은 자리창 안에서 변이 보유율과 메틸이 같이 움직이나")

show("(3) 변이 요약(희귀 변이 수·CpG 변이 등) x 메틸 — 같은 창",
     PR[PR.same_win & (((PR.varA in NONMETH if False else PR.varA.isin(NONMETH)) &
                        PR.varB.isin(METH)) |
                       (PR.varB.isin(NONMETH) & PR.varA.isin(METH)))])

show("(4) 멀리 떨어진 두 창 사이 — 한쪽 끝 값이 반대쪽 끝을 아나",
     PR[(~PR.same_win) & (PR.winA != "HAP") & (PR.winB != "HAP") &
        (~PR.same_var) & PR.varA.isin(METH) & PR.varB.isin(METH)])

show("(5) 하플로타입 상수(길이·경계) x 자리별 항목",
     PR[isHAP & (~PR.same_var) &
        ((PR.varA.isin({"N", "Bp", "Bd", "n_low", "f_hap", "dN", "N_other"})) |
         (PR.varB.isin({"N", "Bp", "Bd", "n_low", "f_hap", "dN", "N_other"})))],
     k=16)

show("(6) 기술 지표 x 무엇이든 — 크면 인공물 의심",
     PR[PR.varA.isin(TECH) | PR.varB.isin(TECH)], k=12,
     note="깊이·관문 통과율 같은 것이 생물학 항목과 붙으면 그 항목은 의심해야 한다")

LOG.close()
print("done")
