# -*- coding: utf-8 -*-
# 틀(축 사각형) 기준 자리잡기 -, 배치 규칙
#
# 왜 : 조립기가 칸(그림 전체 상자)으로 붙이면 y 축 이름 길이가 다른 순간
# 그래프 선이 어긋난다. ⇒ 틀로 맞춘다.
#
# 이 파일이 하는 일 ⒜ 지금 패널들에서 틀과 네 여백을 잰다(여백은 글자가 정하므로
# 배율과 무관하다) ⒝ 줄마다 틀 폭,높이를 하나로 정한다 ⒞ SPEC_FRAME json 을 낸다.
# 그 json 을 넣고 다시 그리면 틀이 그 크기로 못박히고, 조립기는 틀로 붙인다.
#
# 기준 값 :
# 맨 왼쪽 틀 x 41.90, 오른쪽 끝 408.17, 패널 글자 19.42(= 41.90 - 22.48), 틀 위 -1.53
import io
import json
import os
import re
import sys

SP = os.path.dirname(os.path.abspath(__file__))
NL = chr(10)
NUM = r"[-+]?[\d.]+"
X_L, X_R = 41.90, 408.17        # 틀 왼쪽 선, 오른쪽 선 (모든 줄 공통)
GAP_X, GAP_Y = 22.0, 8.0    # 줄 사이 13 -> 8 # 2026-10-19 : 줄 사이를 26 -> 13 # 2026-10-19 : 칸이 붙어 보여 12 -> 22 # 칸 사이 가로 틈(틀 사이가 아니라 여백 사이), 줄 사이
INPAD = 6.0
# 칸마다 가로 몫을 손으로 기울인다 - {파일: 가중치}
WEIGHT = {}                     # 한 칸 안의 작은 틈


def frame_of(path):
    """패널 svg → (폭, 높이, 틀 x0, 틀 y0, 틀 폭, 틀 높이) pt"""
    s = io.open(path, encoding="utf8").read()
    W = float(re.search(r'width="([\d.]+)pt"', s[:400]).group(1))
    H = float(re.search(r'height="([\d.]+)pt"', s[:400]).group(1))
    i = s.find('<g id="axes_1">')
    if i < 0:
        return None
    m = re.search(r'<path d="M\s*(%s)\s+(%s)\s*(.*?)"' % (NUM, NUM), s[i:i + 4000], re.S)
    if not m:
        return None
    xs, ys = [float(m.group(1))], [float(m.group(2))]
    for a, b in re.findall(r"L\s*(%s)\s+(%s)" % (NUM, NUM), m.group(3)):
        xs.append(float(a))
        ys.append(float(b))
    # figure_1 에 통째 옮김(translate)이 걸린 패널이 있다
    # (S5 의 A - 왼쪽 빈 여백을 없애려고 밀었다).
    # 그 옮김을 안 더하면 틀 자리를 옛 자리로 읽어 칸이 쪽 밖으로 끌려 나간다.
    tm = re.search(r'<g id="figure_1" transform="translate\(([-\d.]+),\s*([-\d.]+)\)"', s)
    dx, dy = (float(tm.group(1)), float(tm.group(2))) if tm else (0.0, 0.0)
    return W, H, min(xs) + dx, min(ys) + dy, max(xs) - min(xs), max(ys) - min(ys)


def margins(path):
    """(왼, 오른, 위, 아래, 틀폭, 틀높이)"""
    f = frame_of(path)
    if f is None:
        return None
    W, H, fx, fy, fw, fh = f
    return fx, W - (fx + fw), fy, H - (fy + fh), fw, fh


def plan(rows, panels_dir, name):
    """줄마다 틀 폭·높이를 정해 {stem__label: [fw, fh]} 로 돌려준다"""
    out, report = {}, []
    for ri, row in enumerate(rows):
        cells = []                     # [(파일들 [[...]], letter)]
        for f, lab in row:
            cells.append((f, lab))
        # ── 한 줄의 칸들을 평평하게 펴서 (COL 은 세로로 쌓이므로 폭은 한 칸)
        flat = []                      # [(파일, 칸번호, 줄안의 줄번호)]
        for ci, (f, lab) in enumerate(cells):
            if f == "COL":
                for sub_f, _l in lab:
                    flat.append((sub_f, ci))
            elif isinstance(f, list):          # 격자
                for rr in f:
                    for sub in (rr if isinstance(rr, list) else [rr]):
                        flat.append((sub, ci))
            else:
                flat.append((f, ci))
        # ── 칸마다 가장 넓은 왼/오른 여백
        ncell = len(cells)
        mg = {}
        for f, ci in flat:
            p = os.path.join(panels_dir, f + ".svg")
            if not os.path.exists(p):
                continue
            m = margins(p)
            if m is None:
                continue
            mg[f] = m
        if not mg:
            continue
        # ── 세로칸(column) 을 제대로 센다 : COL 은 한 칸, 격자는 그 안에서 나뉜다
        colspecs = []                  # [(파일들, 왼여백, 오른여백, 틀폭)]
        for ci, (f, lab) in enumerate(cells):
            if f == "COL":
                fs = [sf for sf, _l in lab if sf in mg]
            elif isinstance(f, list):
                grid = [rr if isinstance(rr, list) else [rr] for rr in f]
                ncol = max(len(rr) for rr in grid)
                for j in range(ncol):
                    fs = [rr[j] for rr in grid if j < len(rr) and rr[j] in mg]
                    if fs:
                        colspecs.append((fs, max(mg[x][0] for x in fs),
                                         max(mg[x][1] for x in fs),
                                         max(mg[x][4] for x in fs)))
                continue
            else:
                fs = [f] if f in mg else []
            if fs:
                colspecs.append((fs, max(mg[x][0] for x in fs),
                                 max(mg[x][1] for x in fs),
                                 max(mg[x][4] for x in fs)))
        if not colspecs:
            continue
        ncols = len(colspecs)
        between = sum(colspecs[j - 1][2] + GAP_X + colspecs[j][1]
                      for j in range(1, ncols))
        avail = (X_R - X_L) - between
        def _w(fs, w):
            return w * max(WEIGHT.get(f, 1.0) for f in fs)
        tot = sum(_w(c[0], c[3]) for c in colspecs) or 1.0
        FWS = {}
        for fs, _L, _R, w in colspecs:
            for f in fs:
                FWS[f] = avail * _w(fs, w) / tot
        hs = sorted(mg[f][5] * (FWS[f] / mg[f][4]) for f in FWS)
        FH = hs[len(hs) // 2]
        FW = avail / max(1, ncols)
        for f in FWS:
            out[f] = [round(FWS[f], 3), round(FH, 3)]
        report.append("  줄 %d : 칸 %d · 가로칸 %d · 틀 %.1f x %.1f"
                      % (ri + 1, ncell, ncols, FW, FH))
    print("## %s" % name)
    print(NL.join(report))
    return out


if __name__ == "__main__":
    print("이 파일은 조립기가 불러 쓴다 (plan).")
