# -*- coding: utf-8 -*-
# The shared smoothing for the within-unit methylation curves: FigG13, FigH8 and FigI2 all use this same function.
#   Why: a curve that looks different from figure to figure for the same array misleads the reader
#        (a window of 25 CpG sites and a window of +-bp diverge to a rho of 0.51 in the gorilla).
#   The window width is +-half bp, with half = min(50, max(5, unit length / 20)).
#   Where the calls within the window fall short of minc (40), the window is widened at that position alone until minc is met;
#     the window width changes, not the data.
#     Where widening still cannot meet it, nothing is drawn: no data is invented.
import numpy as np


def half_bp(L):
    return int(min(50, max(5, L / 20.0)))


def smooth_bp(pos, n, k, L, half=None, minc=40, maxmul=16):
    half = half_bp(L) if half is None else int(half)
    N = np.zeros(L + 2)
    K = np.zeros(L + 2)
    N[np.clip(pos, 0, L + 1)] = n
    K[np.clip(pos, 0, L + 1)] = k
    cn = np.concatenate([[0.0], np.cumsum(N)])
    ck = np.concatenate([[0.0], np.cumsum(K)])
    xs, ys, ws = [], [], []
    for p in range(L + 2):
        h = half
        while True:
            a, b = max(0, p - h), min(L + 1, p + h)
            tot = cn[b + 1] - cn[a]
            if tot >= minc or h >= half * maxmul:
                break
            h *= 2
        if tot >= minc:
            xs.append(p)
            ys.append((ck[b + 1] - ck[a]) / tot)
            ws.append(h)
    return np.array(xs), np.array(ys), np.array(ws)
