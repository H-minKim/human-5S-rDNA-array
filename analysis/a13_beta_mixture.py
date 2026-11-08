# -*- coding: utf-8 -*-
"""A13 — 유닛 CpG 메틸 상태를 두 성분 beta-binomial 혼합으로 판정 (β<0.5 임시 기준 대체).

자료: METH_READUNIT in_final=Y · n_cpg>=10 · (n_meth, n_cpg) 를 그대로 사용 (CpG 수가 적으면 판정이 약해짐)
묶음: ONT = dorado 판본별 (q01 공변량) · PacBio = jasmine / fibertools0.4.2 따로 (m6A 없는 2개체는 CpG 는 있음 → 'PB_noM6A' 따로)
EM: 성분 L(저메틸)·H(고메틸) 각 BetaBinomial(n, a, b), M 단계는 가중 로그우도 수치 최대화.
판정: P(L) >= 0.9 → L · P(L) <= 0.1 → H · 사이 → U(보류, 분석에서 제외)
출력 a13_unit_state_gmm.tsv.gz · a13_mixture_params.tsv · a13_beta_mixture.png
"""
import io, os, sys, gzip, collections
import numpy as np
from scipy.stats import betabinom
from scipy.optimize import minimize
from scipy.special import logsumexp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# Paths to files inside the repository are resolved relative to `__file__`
#  (`a_lib` 이 이 폴더에 같이 예치돼 있다) · 자료 뿌리는 `a_lib` 이 환경변수로 받는다
#  변수 표는 `README.md` §2a 에 있다.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import a_lib as L

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False
#  여기에 쓴다 (`RDNA_ACCESS_DIR`) — 발행 자료 폴더를 가리키면
#  `a13_unit_state_gmm.tsv.gz` 와 `a13_mixture_params.tsv` 를 덮는다.
D = L.D
COV = {}
for l in io.open(os.path.join(D, "q01_sample_qc.tsv"), encoding="utf8").readlines()[1:]:
    x = l.rstrip("\n").split("\t")
    COV[(x[0], x[1])] = x[3]

U = []
bad = 0
for src, meth in L.METH.items():
    for node in ("w3", "iB", "iA"):
        p = os.path.join(L.R, "08_methylation_s5", node, meth, "METH_READUNIT_%s.tsv.gz" % L.HOST[node])
        with gzip.open(p, "rt", encoding="utf8") as f:
            ix = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
            for l in f:
                x = l.split("\t", 12)
                if x[ix["data_source"]] != src or x[ix["in_final"]] != "Y":
                    continue
                n = int(x[ix["n_cpg"]]); m = int(x[ix["n_meth"]]); b = float(x[ix["beta"]])
                if n < 10:
                    continue
                if abs(m / n - b) > 0.006:
                    bad += 1
                smp = x[ix["sample"]]
                if src == "UCSC_ONT_FIBERSEQ":
                    grp = "ONT " + COV.get((src, smp), "dorado ?")
                else:
                    grp = ("PB fibertools0.4.2" if smp in L.PB_EXCL else
                           "PB noM6A" if smp in ("HG03672", "HG02624") else "PB jasmine2.7")
                U.append((src, smp, x[ix["read_id"]], x[ix["slot"]], x[ix["array_pos"]], n, m, b, grp))
print("유닛 %d · β≠n_meth/n_cpg %d" % (len(U), bad))


def fit(n, k, iters=60):
    pars = [(2.0, 10.0), (10.0, 2.0)]
    w = 0.3
    prev = -np.inf
    for it in range(iters):
        l0 = np.log(w) + betabinom.logpmf(k, n, *pars[0])
        l1 = np.log(1 - w) + betabinom.logpmf(k, n, *pars[1])
        tot = np.logaddexp(l0, l1)
        ll = tot.sum()
        r = np.exp(l0 - tot)
        w = min(max(r.mean(), 1e-4), 1 - 1e-4)
        newp = []
        for c, rr in ((0, r), (1, 1 - r)):
            def nll(th):
                a, b = np.exp(th)
                return -(rr * betabinom.logpmf(k, n, a, b)).sum()
            res = minimize(nll, np.log(pars[c]), method="L-BFGS-B", bounds=[(-5, 8), (-5, 8)])
            newp.append(tuple(np.exp(res.x)))
        pars = newp
        if abs(ll - prev) < 1e-3:
            break
        prev = ll
    l0 = np.log(w) + betabinom.logpmf(k, n, *pars[0])
    l1 = np.log(1 - w) + betabinom.logpmf(k, n, *pars[1])
    pL = np.exp(l0 - np.logaddexp(l0, l1))
    return w, pars, pL, ll, it + 1


groups = sorted(set(u[8] for u in U))
params = []
PL = {}
fig, axes = plt.subplots(int(np.ceil(len(groups) / 3)), 3, figsize=(15, 3.6 * int(np.ceil(len(groups) / 3))))
axes = axes.ravel()
for gi, grp in enumerate(groups):
    idx = [i for i, u in enumerate(U) if u[8] == grp]
    n = np.array([U[i][5] for i in idx]); k = np.array([U[i][6] for i in idx]); b = np.array([U[i][7] for i in idx])
    w, pars, pL, ll, nit = fit(n, k)
    for i, p in zip(idx, pL):
        PL[i] = p
    mL = pars[0][0] / sum(pars[0]); mH = pars[1][0] / sum(pars[1])
    nmed = int(np.median(n))
    grid = np.arange(nmed + 1)
    l0 = np.log(w) + betabinom.logpmf(grid, nmed, *pars[0]); l1 = np.log(1 - w) + betabinom.logpmf(grid, nmed, *pars[1])
    pg = np.exp(l0 - np.logaddexp(l0, l1))
    cut50 = grid[np.argmax(pg < 0.5)] / nmed
    lo90 = grid[np.argmax(pg < 0.9)] / nmed; hi10 = grid[np.argmax(pg <= 0.1)] / nmed
    st = np.where(pL >= 0.9, "L", np.where(pL <= 0.1, "H", "U"))
    old = np.where(b < 0.5, "L", "H")
    agree = float(np.mean(old[st != "U"] == st[st != "U"])) if (st != "U").any() else float("nan")
    params.append((grp, len(idx), len({U[i][1] for i in idx}), w, mL, mH, pars[0][0], pars[0][1], pars[1][0], pars[1][1],
                   nmed, cut50, lo90, hi10, float((st == "L").mean()), float((st == "U").mean()), float((b < 0.5).mean()), agree, nit))
    ax = axes[gi]
    ax.hist(b[st == "L"], bins=50, range=(0, 1), color="#c0392b", alpha=0.7, label="L (P≥0.9)")
    ax.hist(b[st == "H"], bins=50, range=(0, 1), color="#2c3e50", alpha=0.7, label="H (P≤0.1)")
    ax.hist(b[st == "U"], bins=50, range=(0, 1), color="#f39c12", alpha=0.9, label="보류 U")
    ax.axvline(0.5, color="grey", ls=":", lw=1)
    ax.set_title("%s · 유닛 %d · L평균 %.2f / H평균 %.2f\n보류 %.1f%% · n_cpg 중앙 %d 에서 경계 β %.2f"
                 % (grp, len(idx), mL, mH, 100 * (st == "U").mean(), nmed, cut50), fontsize=9)
    ax.legend(fontsize=7)
for ax in axes[len(groups):]:
    ax.axis("off")
fig.suptitle("5S 유닛 CpG β — 두 성분 beta-binomial 혼합 판정 (점선 = 옛 기준 0.5)", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(D, "a13_beta_mixture.png"), dpi=130)

with io.open(os.path.join(D, "a13_mixture_params.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("group\tn_unit\tn_sample\tpi_L\tmean_L\tmean_H\taL\tbL\taH\tbH\tn_cpg_median\tbeta_cut_P50\tbeta_P90_L_upper\tbeta_P10_H_lower\t"
            "frac_L\tfrac_U\tfrac_L_old05\tagree_old_on_called\titers\n")
    for r in params:
        g.write("\t".join(str(v) if isinstance(v, (str, int)) else "%.4f" % v for v in r) + "\n")
for r in params:
    print("%-22s 유닛 %6d 개체 %2d · πL %.3f · 평균 L %.3f H %.3f · n중앙 %d: 경계 β %.2f (L확실 <%.2f, H확실 >=%.2f) · L %.3f · 보류 %.3f · 옛 L %.3f · 옛기준 일치 %.4f · 반복 %d" %
          (r[0], r[1], r[2], r[3], r[4], r[5], r[10], r[11], r[12], r[13], r[14], r[15], r[16], r[17], r[18]))

with gzip.open(os.path.join(D, "a13_unit_state_gmm.tsv.gz"), "wt", encoding="utf8") as g:
    g.write("src\tsample\tread_id\tslot\tarray_pos\tn_cpg\tn_meth\tbeta\tgroup\tpL\tstate\n")
    for i, u in enumerate(U):
        p = PL[i]
        st = "L" if p >= 0.9 else "H" if p <= 0.1 else "U"
        g.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%s\t%.4f\t%s\n" % (u[:8] + (u[8], p, st)))
print("기록 a13_unit_state_gmm.tsv.gz")
