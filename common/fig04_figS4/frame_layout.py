# -*- coding: utf-8 -*-
# Placement from the frame, that is the axes rectangle, and the layout rules.
#
# Why: when the assembler joins panels on their cells, the whole figure box, a difference in the length of the y-axis label
# puts the plot lines out of line. The frame is used instead.
#
# What this file does: a) measure the frame and its four margins in the current panels (the margins are set by the type and so do not depend on scale);
# b) fix one frame width and height per row; c) write the SPEC_FRAME json.
# Feeding that json back in and redrawing fixes the frames at that size, and the assembler then joins them on the frame.
#
# The reference values:
# the leftmost frame at x 41.90, the right edge at 408.17, the panel letter at 19.42 (41.90 - 22.48), and the top of the frame at -1.53.
import io
import json
import os
import re
import sys

SP = os.path.dirname(os.path.abspath(__file__))
NL = chr(10)
NUM = r"[-+]?[\d.]+"
X_L, X_R = 41.90, 408.17        # the left and right frame lines, shared by every row
GAP_X, GAP_Y = 22.0, 8.0    # the horizontal gap between cells, measured between their margins rather than between their frames, and the gap between rows
INPAD = 6.0
# the horizontal share of each cell, tilted by hand: {file: weight}
WEIGHT = {}                     # the small gaps within one cell


def frame_of(path):
    """panel svg -> (width, height, frame x0, frame y0, frame width, frame height) pt"""
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
    # a panel in figure_1 carries a whole-panel translate
    # (A of S5, pushed to remove the empty margin on its left).
    # Without adding that translate in, the frame is read at its old position and the cell is dragged off the page.
    tm = re.search(r'<g id="figure_1" transform="translate\(([-\d.]+),\s*([-\d.]+)\)"', s)
    dx, dy = (float(tm.group(1)), float(tm.group(2))) if tm else (0.0, 0.0)
    return W, H, min(xs) + dx, min(ys) + dy, max(xs) - min(xs), max(ys) - min(ys)


def margins(path):
    """(left, right, top, bottom, frame width, frame height)"""
    f = frame_of(path)
    if f is None:
        return None
    W, H, fx, fy, fw, fh = f
    return fx, W - (fx + fw), fy, H - (fy + fh), fw, fh


def plan(rows, panels_dir, name):
    """set the frame width and height per row and return {stem__label: [fw, fh]}"""
    out, report = {}, []
    for ri, row in enumerate(rows):
        cells = []                     # [([[files, ...]], letter)]
        for f, lab in row:
            cells.append((f, lab))
        # -- flatten the cells of one row; a COL stacks vertically, so it counts as one cell of width
        flat = []                      # [(file, cell number, row number within the cell)]
        for ci, (f, lab) in enumerate(cells):
            if f == "COL":
                for sub_f, _l in lab:
                    flat.append((sub_f, ci))
            elif isinstance(f, list):          # a grid
                for rr in f:
                    for sub in (rr if isinstance(rr, list) else [rr]):
                        flat.append((sub, ci))
            else:
                flat.append((f, ci))
        # -- the widest left and right margin per cell
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
        # -- count the columns properly: a COL is one column, and a grid divides within it
        colspecs = []                  # [(files, left margin, right margin, frame width)]
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
