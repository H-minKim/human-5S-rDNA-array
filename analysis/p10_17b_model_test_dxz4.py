# -*- coding: utf-8 -*-
# A test of the claim:
#   5S makes its transition by changing the proportion of bimodally separated weakly and highly methylated units,
#    whereas in DXZ4 the methylation level itself rises in steps of about 0.3 at the level of the unit.
#
# The test that separates these already exists here: M against S in p10_17_model_test.py.
#     M: two components with mu0 and mu1 fixed, and a mixing fraction w_b per band = the proportion changes.
#     S: one component with a mean mu_b per band = the level shifts.
#     Both are beta-binomial with a shared overdispersion rho. delta BIC = BIC_S - BIC_M, and a positive value favours M.
#   That code reads the per-base P10_CPG_ALL.tsv.gz, which holds only 5S and D4Z4,
#     so it cannot test DXZ4, RNU2 and SST1.
#   fit() uses only (k, n) per unit, and those are n_meth and n_cpg in P10_FRAME_ALL.tsv,
#     so the same function can be imported unchanged and run on all six arrays.
#
# Two questions are asked separately:
#   [A] within an array, with the position band as the band: how the change proceeds along the array.
#   [B] between chromosomes, with the chromosome as the band: the form of the difference in level between chromosomes.
#       The claim that DXZ4's overall methylation rises in steps belongs to [B].
#   Pooling a female's two X chromosomes would itself create a spurious mixture in DXZ4, so the state is split using the sex table and run separately as well.
#
# The canonical tables are not touched; output goes only to 02_output_tables/claim_test_MvsS/ .
import io, os, math, json, collections, statistics
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, expit, logit

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


H = _need("RDNA_OUT_ROOT")
SRC = os.path.join(H, "02_output_tables", "P10")
SEXP = os.path.join(_need("RDNA_WORK_ROOT"), "P10", "out", "P10_SAMPLES.tsv")
OUT = os.path.join(H, "02_output_tables", "claim_test_MvsS")
os.makedirs(OUT, exist_ok=True)
FAMS = ("5S", "D4Z4_4q", "D4Z4_10q", "DXZ4", "RNU2", "SST1")
DIRN = {"5S": "distal", "D4Z4_4q": "proximal", "D4Z4_10q": "proximal",
        "DXZ4": "proximal", "RNU2": "proximal", "SST1": "proximal"}
MINCPG = 8   # not used; the filter is _complete()


# ---------------- taken from the original code unchanged
def lbb(k, n, mu, rho):
    a = mu * (1 - rho) / rho
    b = (1 - mu) * (1 - rho) / rho
    return (gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
            + gammaln(k + a) + gammaln(n - k + b) - gammaln(n + a + b)
            - gammaln(a) - gammaln(b) + gammaln(a + b))


def fit(k, n, bi, B):
    k, n, bi = np.array(k, float), np.array(n, float), np.array(bi, int)
    N = len(k)

    def nll_S(t):
        mu = expit(t[:B])[bi]
        rho = expit(t[B])
        return -np.sum(lbb(k, n, np.clip(mu, 1e-4, 1 - 1e-4), rho))

    def nll_M(t):
        m0, m1 = expit(t[0]), expit(t[1])
        rho = expit(t[2])
        w = expit(t[3:3 + B])[bi]
        l0 = lbb(k, n, np.clip(m0, 1e-4, 1 - 1e-4), rho) + np.log(np.clip(1 - w, 1e-12, 1))
        l1 = lbb(k, n, np.clip(m1, 1e-4, 1 - 1e-4), rho) + np.log(np.clip(w, 1e-12, 1))
        m = np.maximum(l0, l1)
        return -np.sum(m + np.log(np.exp(l0 - m) + np.exp(l1 - m)))

    mb = []
    for b in range(B):
        s = bi == b
        mb.append(k[s].sum() / max(1, n[s].sum()) if s.any() else 0.5)
    rS = minimize(nll_S, np.r_[logit(np.clip(mb, 0.02, 0.98)), logit(0.05)], method="L-BFGS-B")
    best = None
    for m0, m1 in ((0.08, 0.8), (0.3, 0.8), (0.45, 0.8), (0.1, 0.6)):
        w0 = [logit(min(0.98, max(0.02, (x - m0) / (m1 - m0)))) for x in mb]
        rM = minimize(nll_M, np.r_[logit(m0), logit(m1), logit(0.03), w0], method="L-BFGS-B")
        if best is None or rM.fun < best.fun:
            best = rM
    rM = best
    bicS = 2 * rS.fun + (B + 1) * math.log(N)
    bicM = 2 * rM.fun + (B + 3) * math.log(N)
    m0, m1 = sorted([expit(rM.x[0]), expit(rM.x[1])])
    return (bicS - bicM, m0, m1, expit(rM.x[2]), expit(rS.x[B]),
            [expit(x) for x in rS.x[:B]],
            sorted([expit(x) for x in rM.x[3:3 + B]], key=lambda q: 0) if False else
            [expit(x) for x in rM.x[3:3 + B]], (expit(rM.x[0]) > expit(rM.x[1])))


MINC, COVMIN = 20, 0.9          # complete units, by the same rule as fig39, fig43 and fig47
_TOTB = collections.defaultdict(collections.Counter)
with io.open(os.path.join(SRC, "P10_FRAME_ALL.tsv"), encoding="utf8") as _f:
    _h = _f.readline().rstrip(chr(10)).split(chr(9))
    _j = {k: _h.index(k) for k in _h}
    for _L in _f:
        _g = _L.rstrip(chr(10)).split(chr(9))
        if _g[_j['family']] not in FAMS:
            continue
        for _t in _g[_j['bins']].split(';'):
            if _t:
                _TOTB[_g[_j['family']]][int(_t.split(':')[0])] += 1
CPGB = {fam: {b for b, v in c.items() if v >= 0.10 * max(c.values())}
        for fam, c in _TOTB.items()}
print('  흔한 칸 : %s' % ' '.join('%s %d' % (k, len(v)) for k, v in sorted(CPGB.items())))


def _complete(g, ix):
    fam = g[ix['family']]
    if int(g[ix['n_cpg']]) < MINC:
        return False
    bs = {int(t.split(':')[0]) for t in g[ix['bins']].split(';') if t}
    cp = CPGB.get(fam) or set()
    return bool(cp) and len(bs & cp) / float(len(cp)) >= COVMIN


# ---------------- data: (k, n, position, chromosome) per unit
SEX = {}
for i, L in enumerate(io.open(SEXP, encoding="utf8")):
    if i:
        f = L.rstrip("\n").split("\t"); SEX[f[0]] = f[2]
U = collections.defaultdict(list)          # fam -> [(k, n, p_prox, p_dist, sample, hap)]
with io.open(os.path.join(SRC, "P10_FRAME_ALL.tsv"), encoding="utf8") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ix = {k: hdr.index(k) for k in hdr}
    for L in f:
        g = L.rstrip("\n").split("\t")
        fam = g[ix["family"]]
        if fam not in FAMS:
            continue
        u, nu = g[ix["unit_from_distal"]], g[ix["n_units"]]
        if not u.lstrip("-").isdigit() or not nu.isdigit():
            continue
        n = int(g[ix["n_cpg"]])
        if not _complete(g, ix):      # complete units only
            continue
        d, N = int(u), int(nu)
        if d > N:
            continue
        U[fam].append((int(g[ix["n_meth"]]), n, N + 1 - d, d, g[ix["sample"]], g[ix["hap"]]))
for fam in FAMS:
    print("  %-9s 유닛 %s" % (fam, format(len(U[fam]), ",")))

# ---------------- the DXZ4 chromosome state, taken from the sex table and not from a threshold on the data
mn = collections.defaultdict(lambda: [0, 0])
for k, n, pp, pd, s, h in U["DXZ4"]:
    a = mn[(s, h)]; a[0] += k; a[1] += n
bys = collections.defaultdict(list)
for (s, h), a in mn.items():
    bys[s].append((a[0] / float(a[1]), h))
DXG = {}
for s, v in bys.items():
    v.sort(reverse=True)
    if SEX.get(s) == "XY":
        for m, h in v:
            DXG[(s, h)] = "Xhigh"
    else:
        DXG[(s, v[0][1])] = "Xhigh"
        for m, h in v[1:]:
            DXG[(s, h)] = "Xlow"

BANDS_POS = {"5S": [(1, 3), (4, 6), (7, 10), (11, 15), (16, 20), (21, 25), (26, 30)],
             "D4Z4_4q": [(1, 2), (3, 4), (5, 6), (7, 8), (9, 11), (12, 15)],
             "D4Z4_10q": [(1, 2), (3, 4), (5, 6), (7, 8), (9, 11), (12, 15)],
             "DXZ4": [(1, 2), (3, 4), (5, 6), (7, 8), (9, 11), (12, 15)],
             "RNU2": [(1, 2), (3, 4), (5, 6), (7, 8), (9, 11), (12, 15)],
             "SST1": [(1, 2), (3, 4), (5, 6), (7, 8), (9, 11), (12, 15)]}
LINES = ["test\tfamily\tscope\tn_unit\tn_band\tdBIC_S_minus_M\twinner\tM_mu0\tM_mu1\tM_gap"
         "\tM_rho\tS_rho\tS_mu_min\tS_mu_max\tS_range\tband_detail"]


def run(test, fam, rows, bands, names, scope):
    K, N, BI = [], [], []
    for r in rows:
        b = bands(r)
        if b is not None:
            K.append(r[0]); N.append(r[1]); BI.append(b)
    B = len(names)
    if len(K) < 200 or len(set(BI)) < 3:
        print("    %-9s %-22s (표본 부족 : 유닛 %d · 구간 %d)" % (fam, scope, len(K), len(set(BI))))
        return
    d, m0, m1, rhoM, rhoS, muS, wM, flip = fit(K, N, BI, B)
    w = [(1 - x) if flip else x for x in wM]
    win = "M (mixture ratio)" if d > 0 else "S (level shift)"
    det = ";".join("%s:w=%.3f:mu=%.3f" % (names[i], w[i], muS[i]) for i in range(B))
    print("    %-9s %-22s 유닛 %6d  dBIC %+9.1f  -> %-18s  M %.3f/%.3f (gap %.3f)  S %.3f-%.3f"
          % (fam, scope, len(K), d, win, m0, m1, m1 - m0, min(muS), max(muS)))
    LINES.append("%s\t%s\t%s\t%d\t%d\t%.1f\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s"
                 % (test, fam, scope, len(K), B, d, win, m0, m1, m1 - m0, rhoM, rhoS,
                    min(muS), max(muS), max(muS) - min(muS), det))


print()
print("=" * 118)
print("[A] array 안 — 위치 구간이 band   (ΔBIC > 0 = M 섞임비율 · < 0 = S 수준이동)")
print("=" * 118)
for fam in FAMS:
    bs = BANDS_POS[fam]
    names = ["%d-%d" % t for t in bs]
    key = 3 if DIRN[fam] == "distal" else 2

    def mk(r, bs=bs, key=key):
        p = r[key]
        for i, (a, z) in enumerate(bs):
            if a <= p <= z:
                return i
        return None
    run("A_position", fam, U[fam], mk, names, "all chromosomes")
    if fam == "DXZ4":
        for st in ("Xhigh", "Xlow"):
            run("A_position", fam, [r for r in U[fam] if DXG.get((r[4], r[5])) == st], mk, names, st)

print()
print("=" * 118)
print("[B] 염색체 사이 — 염색체가 band   (DXZ4 의 '전체 메틸이 계단처럼 오른다' 는 여기에 해당)")
print("=" * 118)
for fam in FAMS:
    chs = sorted({(r[4], r[5]) for r in U[fam]})
    cnt = collections.Counter((r[4], r[5]) for r in U[fam])
    chs = [c for c in chs if cnt[c] >= 60]
    idx = {c: i for i, c in enumerate(chs)}
    names = ["%s_%s" % c for c in chs]
    run("B_chromosome", fam, [r for r in U[fam] if (r[4], r[5]) in idx],
        lambda r, idx=idx: idx.get((r[4], r[5])), names, "%d chromosomes" % len(chs))

io.open(os.path.join(OUT, "MODEL_TEST_ALL6.tsv"), "w", encoding="utf8",
        newline="\n").write("\n".join(LINES) + "\n")
print()
print("쓴 표 : %s" % os.path.join(OUT, "MODEL_TEST_ALL6.tsv"))
