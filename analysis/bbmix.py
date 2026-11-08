# -*- coding: utf-8 -*-
# 베타-이항 1 성분 / 2 성분 혼합 — 가중 우도 판(같은 (k,n) 을 묶어 한 번만 계산)
#   xa1_mixture_lrt.py 의 적합과 수학적으로 같다. 유닛 수가 많아도 (k,n) 가짓수만큼만 계산하므로 빠르다
#   (예: 유닛 93,323 개라도 CpG 가 15~20 이면 (k,n) 가짓수는 수백 개).
#   대조는 xa4_bbmix_check.py 에서 xa1 과 같은 값이 나오는지 실측한다.
import numpy as np
from scipy import optimize
from scipy.special import betaln, gammaln


def pack(k, n):
    """(k,n) 을 묶어 (k_u, n_u, 개수) 로 — 가중 우도용"""
    k = np.asarray(k, dtype=np.int64); n = np.asarray(n, dtype=np.int64)
    key = np.stack([k, n], 1)
    uq, cnt = np.unique(key, axis=0, return_counts=True)
    return uq[:, 0].astype(float), uq[:, 1].astype(float), cnt.astype(float)


def bb_logpmf(k, n, a, b):
    return (gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
            + betaln(k + a, n - k + b) - betaln(a, b))


def fit1(k, n, c):
    def nll(t):
        a, b = np.exp(t)
        if not np.isfinite(a) or not np.isfinite(b) or a <= 0 or b <= 0:
            return 1e18
        return -(c * bb_logpmf(k, n, a, b)).sum()
    p = float(np.clip((c * k).sum() / (c * n).sum(), 1e-4, 1 - 1e-4))
    pb = k / n
    m = float((c * pb).sum() / c.sum())
    v = max(float((c * (pb - m) ** 2).sum() / max(c.sum() - 1, 1)), 1e-8)
    t0 = max(p * (1 - p) / v - 1, 0.5)
    x0 = np.log([max(p * t0, 1e-3), max((1 - p) * t0, 1e-3)])
    r = optimize.minimize(nll, x0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-5, fatol=1e-5))
    return -r.fun, np.exp(r.x)


def fit2(k, n, c, starts=3, seed=0):
    best = (-np.inf, None)
    pb = k / n
    order = np.argsort(pb)
    cw = np.cumsum(c[order]) / c.sum()

    def q(x):      # 가중 분위수
        return float(pb[order][np.searchsorted(cw, x, side="left").clip(0, len(pb) - 1)])
    for s in range(starts):
        q0, q1 = (0.25, 0.80) if s == 0 else (q(0.15 + 0.1 * s), q(0.80 + 0.05 * s))
        q0 = min(max(q0, 1e-3), 0.95); q1 = min(max(q1, q0 + 0.02), 1 - 1e-3)
        conc = 20.0 * (s + 1)
        t0 = np.array([0.0, np.log(q0 * conc), np.log((1 - q0) * conc),
                       np.log(q1 * conc), np.log((1 - q1) * conc)])

        def nll(t):
            w = 1.0 / (1.0 + np.exp(-t[0]))
            a0, b0, a1, b1 = np.exp(t[1:])
            if not np.all(np.isfinite([a0, b0, a1, b1])) or min(a0, b0, a1, b1) <= 0:
                return 1e18
            l0 = bb_logpmf(k, n, a0, b0) + np.log(max(w, 1e-12))
            l1 = bb_logpmf(k, n, a1, b1) + np.log(max(1 - w, 1e-12))
            m = np.maximum(l0, l1)
            return -(c * (m + np.log(np.exp(l0 - m) + np.exp(l1 - m)))).sum()
        r = optimize.minimize(nll, t0, method="Nelder-Mead", options=dict(maxiter=8000, xatol=1e-5, fatol=1e-5))
        if -r.fun > best[0]:
            w = 1.0 / (1.0 + np.exp(-r.x[0]))
            a0, b0, a1, b1 = np.exp(r.x[1:])
            mu0, mu1 = a0 / (a0 + b0), a1 / (a1 + b1)
            if mu0 > mu1:
                w, a0, b0, a1, b1, mu0, mu1 = 1 - w, a1, b1, a0, b0, mu1, mu0
            best = (-r.fun, dict(w_low=w, a0=a0, b0=b0, a1=a1, b1=b1, mu_low=mu0, mu_high=mu1))
    return best


def lr_stat(k, n, c=None, starts=3, seed=0):
    if c is None:
        k, n, c = pack(k, n)
    l1, ab = fit1(k, n, c)
    l2, pr = fit2(k, n, c, starts=starts, seed=seed)
    return 2.0 * (l2 - l1), l1, l2, ab, pr


def boot_null(n_units, a, b, rng):
    """적합된 1 성분에서 같은 깊이로 한 벌 모의 → (k,n,개수)"""
    ks = rng.binomial(n_units.astype(int), rng.beta(a, b, size=len(n_units)))
    return pack(ks, n_units)


def mixture_peaks(pr, grid=None):
    """적합된 혼합이 만드는 분자 수준 p 분포 의 봉우리 수·자리·골깊이비"""
    from scipy import stats as _st
    if grid is None:
        grid = np.linspace(1e-4, 1 - 1e-4, 2001)
    d = (pr["w_low"] * _st.beta.pdf(grid, pr["a0"], pr["b0"])
         + (1 - pr["w_low"]) * _st.beta.pdf(grid, pr["a1"], pr["b1"]))
    up = np.r_[True, d[1:] > d[:-1]]
    dn = np.r_[d[:-1] > d[1:], True]
    pk = np.flatnonzero(up & dn)
    if len(pk) < 2:
        return len(pk), np.nan, (float(grid[pk[0]]) if len(pk) else np.nan), np.nan
    i, j = pk[0], pk[-1]
    ratio = float(d[i:j + 1].min() / min(d[i], d[j]))
    return len(pk), ratio, float(grid[i]), float(grid[j])
