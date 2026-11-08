# -*- coding: utf-8 -*-
# 유닛 안 메틸 곡선의 공용 평활 — FigG13 · FigH8 · FigI2 가 같은 함수 를 쓴다(2026-10-06)
#   왜 : 같은 array 인데 그림마다 곡선이 달라 보이면 독자가 오해한다(사용자 지적).
#        앞서 G13 은 'CpG 자리 25 개' 창, H8 은 '±bp' 창이라 고릴라에서 rho 0.51 까지 벌어졌다.
#   창 폭 : ±half bp, half = min(50, max(5, 유닛길이/20))
#   끝이 끊기던 문제 : 창 안 호출수가 minc(40) 에 못 미치면 그 자리를 안 그렸다 →
#     그 자리에서만 창을 넓혀 minc 를 채운다(자료를 채우는 게 아니라 창 폭만 바꾼다).
#     넓혀도 못 채우면 그리지 않는다(없는 자료를 만들지 않는다).
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
