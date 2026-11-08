# -*- coding: utf-8 -*-
# β 를 칸에 담을 때 쓰는 공용 자(2026-10-04 · 사용자 지적)
#   β = k / n_cpg 라서 값이 n+1 가지뿐 이다. 0~1 을 고정 폭(예 0.02)으로 자르면
#     · 값들이 칸에 고르게 안 떨어져 모양이 왜곡되고
#     · 칸 가운데 를 이어 그리면 오른쪽 끝이 1.0 에 못 미쳐 끝난다(CpG 8 → 마지막 가운데 0.95).
#   그래서 칸 가운데가 0, 1/n, …, 1 이 되도록 경계를 잡는다. n 이 크면(>CAP) 고정 칸으로 돌아간다.
import numpy as np


def value_bins(n_cpg, cap=60):
    """반환 : (edges, centers) — 칸 가운데가 β 가 가질 수 있는 값이 되게"""
    n = int(np.clip(round(float(n_cpg)), 4, cap))
    edges = np.linspace(-0.5 / n, 1 + 0.5 / n, n + 2)
    centers = np.arange(n + 1) / float(n)
    return edges, centers


def density(b, n_cpg, cap=60):
    """밀도(비율 ÷ 칸폭) 와 칸 가운데 — 0 과 1 이 양 끝 가운데가 된다"""
    edges, centers = value_bins(n_cpg, cap)
    h, e = np.histogram(np.asarray(b, dtype=float), bins=edges)
    w = np.diff(e)
    tot = h.sum()
    return (h / tot) / w if tot else h * 0.0, centers


def fraction(b, n_cpg, cap=60):
    """몫(그 칸에 든 read x 유닛 / 전체) 과 칸 가운데 — 2026-10-20

    `density()` 와 무엇이 다른가 : `density()` 는 몫 ÷ 칸폭 이라 면적이 1 인 밀도다.
    사용자가 축 이름을 `Fraction of read-units` 로 정했으므로 그려지는 양도 몫이어야 한다 —
    이름만 바꾸면 틀린 이름 이 된다(§2-12 : 값·자·창·집계단위를 함께 적는다).
    모양은 바뀌지 않는다 : 한 칸 안에서 칸폭 `w = 1/n` 이 상수 이고, 칸마다 y 를 그 칸의 최대에 맞추며
    y 눈금을 적지 않으므로 보이는 것은 그대로다(실측으로 확인한다).
    `density()` 는 그대로 둔다 — 다른 그림(FigA6 · FigJ13 등)이 쓰고 있다.
    """
    edges, centers = value_bins(n_cpg, cap)
    h, _e = np.histogram(np.asarray(b, dtype=float), bins=edges)
    tot = h.sum()
    return (h / tot if tot else h * 0.0), centers


def value_index(b, n_cpg, cap=50):
    """열지도 줄 번호 : β 를 가장 가까운 가능한 값 의 칸으로 (n 이 크면 cap 칸 고정)"""
    n = int(round(float(n_cpg)))
    b = np.asarray(b, dtype=float)
    if n <= cap:
        return np.clip(np.round(b * n).astype(int), 0, n), n + 1
    return np.clip((b * cap).astype(int), 0, cap - 1), cap
