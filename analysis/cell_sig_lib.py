# -*- coding: utf-8 -*-
"""세포주 / HPRC 공용 — 블라인드 격자 스캔 + 호출기 없는 확인검정 (2026-10-12)

통계는 HPRC 본표(zzt_gridwide_rank · zzx_site_confirm_callerfree)와 같은 식이다.
  ① 블라인드(FiberHMM) : 5 bp 크기 칸 x 자리 · 치환 한 번에 이동량 한 벌을 모든 칸에 같이 쓴다 ·
       p_grid = 귀무에서 격자 전체 최대 >= 관측 · p_win = 귀무의 그 창 안(칸 x 자리) 최대 >= 관측
  ② 확인(호출기 없음) : 분자마다 (A, m6A) 자리표를 통째로 원형 이동 · 창 안 자리 x 폭 최대
다른 점(적어 둔다)
  - 칸 사용 문턱 : HPRC 150 사건 -> 여기는 minev (분자 수가 적어서). 사건 수는 이동에 불변이라 편향 없음.
  - 귀무 평균 하한 : 사건 1 개에 해당하는 밀도(1000/(nmol*SM)). 기대 사건 <1 인 곳의 배수는 뜻이 없다.
"""
import numpy as np

PER, HALF = 2240, 1120
LO, HI = -150, 118
NP = HI - LO + 1
I0, I1 = LO + HALF, HI + HALF + 1
SM = 11
SBIN, SMIN, SMAX = 5, 20, 164
SITES = {"TFIIIB": (-28, -10), "ICR": (57, 68)}
SITE_BINS = {"TFIIIB": (30, 39), "ICR": (40, 49)}
WIDTHS = list(range(9, 72, 2))
GAP, FLANKW = 8, 17
REL = np.arange(LO, HI + 1)


def smooth_rows(X):
    h = SM // 2
    P = np.pad(X, ((0, 0), (h, h)))
    cs = np.concatenate([np.zeros((X.shape[0], 1)), np.cumsum(P, 1)], 1)
    return (cs[:, SM:] - cs[:, :-SM]) / SM


def hist_rows(rb, ce, nk, nmol):
    H = np.bincount(rb * PER + ce, minlength=nk * PER).reshape(nk, PER)
    return H[:, I0:I1] * 1000.0 / nmol


def blind(mi, ce, sz, nmol, NB, rng, minev=30, sbin=5):
    """블라인드 격자 스캔."""
    ok = (sz >= SMIN) & (sz <= SMAX) & (ce >= 0) & (ce < PER)
    mi, ce, sz = mi[ok], ce[ok], sz[ok]
    bid = (sz - SMIN) // sbin
    nb_all = (SMAX - SMIN) // sbin + 1
    cnt = np.bincount(bid, minlength=nb_all)
    keep = [b for b in range(nb_all) if cnt[b] >= minev]
    if not keep:
        return None
    remap = -np.ones(nb_all, int)
    remap[keep] = np.arange(len(keep))
    sel = remap[bid] >= 0
    mi, ce, bid = mi[sel], ce[sel], bid[sel]
    rb = remap[bid]
    nk = len(keep)
    OBS = smooth_rows(hist_rows(rb, ce, nk, nmol))
    PROF = np.empty((NB, nk, NP))
    for t in range(NB):
        d = rng.integers(0, PER, nmol)
        PROF[t] = smooth_rows(hist_rows(rb, (ce + d[mi]) % PER, nk, nmol))
    BASE = PROF.mean(0)
    FLOOR = 1000.0 / (nmol * SM)
    BF = np.maximum(BASE, FLOOR)
    ENR = OBS / BF
    NR = PROF / BF
    bins = [(SMIN + b * sbin, SMIN + b * sbin + sbin - 1) for b in keep]
    gmax = NR.reshape(NB, -1).max(1)
    thr = float(np.percentile(gmax, 95))
    out = dict(bins=bins, ENR=ENR, OBS=OBS, BASE=BASE, thr=thr, gmax=gmax, NB=NB, nk=nk)
    res = {}
    for nm, (a, b) in SITES.items():
        m = (REL >= a) & (REL <= b)
        obs = float(ENR[:, m].max())
        j = np.unravel_index(np.argmax(np.where(m[None, :], ENR, -1e9)), ENR.shape)
        nullw = NR[:, :, m].reshape(NB, -1).max(1)
        res[nm] = dict(
            obs=obs, bin=bins[j[0]], pos=int(REL[j[1]]),
            p_win=(float((nullw >= obs).sum()) + 1) / (NB + 1),
            p_grid=(float((gmax >= obs).sum()) + 1) / (NB + 1))
        # 사전 지정 크기칸(HPRC 가 정한 칸) 안만
        lo_s, hi_s = SITE_BINS[nm]
        kb = [i for i, bb in enumerate(bins) if bb[0] >= lo_s and bb[1] <= hi_s]
        if kb:
            obsb = float(ENR[np.ix_(kb, np.where(m)[0])].max())
            nullb = NR[:, kb, :][:, :, m].reshape(NB, -1).max(1)
            res[nm]["obs_bin"] = obsb
            res[nm]["p_win_bin"] = (float((nullb >= obsb).sum()) + 1) / (NB + 1)
    out["sites"] = res
    # 격자 전체 순위 : 봉우리(15 bp 이상 떨어진) 목록
    flat = []
    for bi in range(nk):
        e = ENR[bi]
        used = []
        for j in np.argsort(e)[::-1]:
            if e[j] <= thr:
                break
            if any(abs(REL[j] - q) < 15 for q in used):
                continue
            used.append(REL[j])
            flat.append((float(e[j]), int(REL[j]), bins[bi]))
    flat.sort(reverse=True)
    out["peaks"] = flat
    return out


def dens_fold(mi, ce, nmol, NB, rng, nmol_total=None):
    """한 사건 집합의 평활 밀도(1000 분자당)와 귀무 평균. mi 는 0..nmol-1."""
    ce = np.asarray(ce)
    mi = np.asarray(mi)
    ok = (ce >= 0) & (ce < PER)
    mi, ce = mi[ok], ce[ok]
    obs = smooth_rows(hist_rows(np.zeros(len(ce), int), ce, 1, nmol))[0]
    base = np.zeros(NP)
    for t in range(NB):
        d = rng.integers(0, PER, nmol)
        base += smooth_rows(hist_rows(np.zeros(len(ce), int), (ce + d[mi]) % PER, 1, nmol))[0]
    base /= NB
    return obs, base


# ── 호출기 없는 방법
def W(cs, c0, c1):
    o = np.full(NP, np.nan)
    i = np.arange(NP)
    x0, x1 = i + c0, i + c1 + 1
    ok = (x0 >= 0) & (x1 <= NP)
    o[ok] = cs[x1[ok]] - cs[x0[ok]]
    return o


def dprof(sa, sm, w):
    h = w // 2
    ca = np.concatenate([[0.0], np.cumsum(sa)])
    cm = np.concatenate([[0.0], np.cumsum(sm)])
    cA, cM = W(ca, -h, h), W(cm, -h, h)
    fA = W(ca, -h - GAP - FLANKW, -h - GAP - 1) + W(ca, h + GAP + 1, h + GAP + FLANKW)
    fM = W(cm, -h - GAP - FLANKW, -h - GAP - 1) + W(cm, h + GAP + 1, h + GAP + FLANKW)
    s = np.full(NP, np.nan)
    ok = np.isfinite(cA) & np.isfinite(fA) & (cA >= 200) & (fA >= 200)
    s[ok] = (fM[ok] / fA[ok]) / ((cM[ok] + 0.5) / (cA[ok] + 1.0))
    return s


def dub_paper(sa, sm, c):
    """그림 Fig5F 의 기하 : core c-7..c+1 · flank c-35..c-27, c+22..c+27."""
    def seg(a, c0, c1):
        i0, i1 = c0 - LO, c1 - LO + 1
        if i0 < 0 or i1 > NP:
            return np.nan
        return a[i0:i1].sum()
    if c + 27 > HI or c + 1 > HI:
        return np.nan
    cA, cM = seg(sa, c - 7, c + 1), seg(sm, c - 7, c + 1)
    fA = np.nansum([seg(sa, c - 35, c - 27), seg(sa, c + 22, c + 27)])
    fM = np.nansum([seg(sm, c - 35, c - 27), seg(sm, c + 22, c + 27)])
    if cA >= 200 and fA >= 200:
        return (fM / fA) / ((cM + 0.5) / (cA + 1.0))
    return float("nan")


def confirm(A, M, NB, rng):
    """창 안 자리 x 폭 최대의 원형 이동 귀무 (분자마다 통째로)."""
    N = len(A)
    sA, sM = A.sum(0).astype(float), M.sum(0).astype(float)
    masks = {nm: (REL >= a) & (REL <= b) for nm, (a, b) in SITES.items()}
    obs = {}
    for nm, m in masks.items():
        best = (-9.0, None, None)
        for w in WIDTHS:
            o = dprof(sA, sM, w)
            k = m & np.isfinite(o)
            if not k.any():
                continue
            j = int(np.nanargmax(np.where(k, o, np.nan)))
            if o[j] > best[0]:
                best = (float(o[j]), int(REL[j]), w)
        obs[nm] = best
    nullmax = {nm: np.full(NB, -9.0) for nm in masks}
    ar = np.arange(NP)[None, :]
    for t in range(NB):
        d = rng.integers(0, NP, N)
        idx = (ar - d[:, None]) % NP
        ra = np.take_along_axis(A, idx, 1).sum(0).astype(float)
        rm = np.take_along_axis(M, idx, 1).sum(0).astype(float)
        for w in WIDTHS:
            p = dprof(ra, rm, w)
            for nm, m in masks.items():
                k = m & np.isfinite(p)
                if k.any():
                    v = float(np.nanmax(np.where(k, p, np.nan)))
                    if v > nullmax[nm][t]:
                        nullmax[nm][t] = v
    res = {}
    for nm in masks:
        o = obs[nm]
        nv = nullmax[nm]
        res[nm] = dict(obs=o[0], pos=o[1], width=o[2],
                       null_med=float(np.median(nv)), null95=float(np.percentile(nv, 95)),
                       p=(float((nv >= o[0]).sum()) + 1) / (NB + 1))
    return res
