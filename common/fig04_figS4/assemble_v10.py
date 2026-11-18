# -*- coding: utf-8 -*-
# The assembler, version 10, working from the frame, that is the axes rectangle.
#
# Joining the panels on their cells, the whole figure box, puts the plot lines out of line wherever the y-axis labels differ in length.
# Version 10 joins them on the frame:
# , within a column the left and right x of the frame are equal (41.90 at the far left, 408.17 at the far right);
# , within a row the top y of the frame is equal and the frame heights are equal;
# , the panel letter sits at -22.48 from the left of that frame and -1.53 above its baseline, as an 8.2 pt capital.
# The frame sizes are fixed by paper_spec_hook, in SPEC_FRAME, and are neither grown nor shrunk here.
import collections
import io
import json
import os
import re
import shutil
import sys

SP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SP)
from frame_layout import X_L, X_R, GAP_X, GAP_Y, INPAD, margins   # noqa


def _need(k):
    """The deposited copy keeps no default: it stops rather than read the wrong place."""
    v = os.environ.get(k)
    if not v:
        raise SystemExit("environment variable %s is not set (the "
                         "deposited copy has no default; see README.md,"
                         " section 'Environment variables')" % k)
    return v


# The output root is taken from OUT_BASE, so that nothing is written outside the output folder.
BASE = (os.environ.get("OUT_BASE")
        or os.path.join(_need("RDNA_ARRAY_METH_DIR"), "03_figures"))
PREV = os.path.join(BASE, "00_final_set_2026-10-12")
OLDD = os.path.join(BASE, "Fig4_submission_set_2026-10-09")
OUTD = os.path.join(BASE, os.environ.get("OUT_SET", "00_시험_423pt_2026-10-18"))
PD = os.path.join(OUTD, "panels")
SRC = os.path.join(SP, os.environ.get("PANEL_SRC", "panelsplit_spec"), "_panels")
os.makedirs(PD, exist_ok=True)
NL = chr(10)
W_TOTAL = float(os.environ.get("SPEC_W", "423.5"))
LAB_DX, LAB_DY, FS_LAB = -22.48, -1.53, 8.2
# where a cell exceeds the page width of 423.5, the right frame line is pulled in by that much
XR_ADJ = [0.0]
# a figure which cannot be redrawn is scaled down whole
# only the cell, the plot, is scaled; the panel letter stays at 8.2 pt
# the 4.5 pt floor on type size is then broken, and that fact is printed and recorded
FIG_SCALE = {"S5": None}        # None means it is fitted automatically to the page width
# shifting individual cells further left, where a long y-axis label runs into the letter: (figure, letter) -> the amount to add
LAB_DX_EXTRA = {("Fig4_main", "N"): -9.0,
                # where the cell letter would touch the y-axis label of the cell to its left
                # (E against "Low-methylation blocks", F against "High-...", J against beta)
                ("S4", "E"): -5.0, ("S4", "F"): -5.0, ("S4", "J"): -2.5,
                # where a cell letter would touch the legend or the y-axis label above its own frame
                ("Fig4_main", "F"): -3.0,
                ("S6", "C"): -3.0, ("S6", "G"): -6.5, ("S6", "H"): -3.0,
                ("S6", "R"): -3.0, ("S6", "S"): -3.0}
# the amount by which a whole cell is pushed right: {(figure, letter): pt}
# the amount by which each row is pulled up: {(figure, row index from 0): pt}
ROW_DY = {("Fig4_main", 2): -10.0, ("Fig4_main", 3): -10.0,
          ("S4", 1): -10.0, ("S4", 2): -10.0,
          ("S6", 2): -10.0, ("S6", 5): -10.0}
CELL_DX = {("S6", "A"): 10.0, ("S6", "F"): 10.0, ("S6", "L"): 10.0,
           ("S6", "Q"): 10.0, ("Fig4_main", "G"): 10.0,
           ("Fig4_main", "H"): 6.0}

# -- the layout table is read from version 9 unchanged, so that only one copy exists
_src = io.open(os.path.join(SP, "assemble_v9_2026-10-18.py"), encoding="utf8").read()
_ns = {"__file__": os.path.join(SP, "assemble_v9_2026-10-18.py")}
exec(compile(_src[:_src.index("def fetch(")], "head", "exec"), _ns)
FIG4, S41, S42, S44 = _ns["FIG4"], _ns["S41"], _ns["S42"], _ns["S44"]


def fetch(b):
    got = []
    for ext in ("svg", "pdf"):
        s = os.path.join(SRC, b + "." + ext)
        for alt in (os.path.join(PREV, "panels", b + "." + ext),
                    os.path.join(OLDD, "panels_2026-10-11", b + "." + ext),
                    os.path.join(OLDD, "panels", b + "." + ext)):
            if not os.path.exists(s):
                s = alt
        if os.path.exists(s):
            shutil.copy2(s, os.path.join(PD, b + "." + ext))
            got.append(ext)
    if "svg" not in got:
        raise SystemExit("no panel: " + b)
    return "pdf" in got


def size_pt(p):
    s = io.open(p, encoding="utf8").read(400)
    return (float(re.search(r'width="([\d.]+)pt"', s).group(1)),
            float(re.search(r'height="([\d.]+)pt"', s).group(1)))


def inner(p, pref):
    s = io.open(p, encoding="utf8").read()
    i = s.find("<svg")
    j = s.find(">", i)
    body = s[j + 1:s.rfind("</svg>")]
    for k in sorted(set(re.findall(r'id="([^"]+)"', body)), key=len, reverse=True):
        nk = "%s_%s" % (pref, k)
        body = body.replace('id="%s"' % k, 'id="%s"' % nk)
        body = body.replace('#%s"' % k, '#%s"' % nk)
        body = body.replace("url(#%s)" % k, "url(#%s)" % nk)
    return body


def columns_of(row):
    """row -> [(letter, [vertically stacked (file, column number) ...])] in column order"""
    cols = []                      # [(letter, [[file, ...] per column], the share of the row width)]
    for f, lab in row:
        if f == "COL":
            sub = [(s_f, s_l) for s_f, s_l in lab]
            cols.append(("COL", sub))
        elif isinstance(f, list):
            grid = [rr if isinstance(rr, list) else [rr] for rr in f]
            cols.append(("GRID", (grid, lab)))
        else:
            cols.append(("P", (f, lab)))
    return cols


def _cell_dx(fig, lab):
    for (f, l), v in CELL_DX.items():
        if f in fig and l == lab:
            return v
    return 0.0


def _extra(fig, lab):
    for (f, l), v in LAB_DX_EXTRA.items():
        if f in fig and l == lab:
            return v
    return 0.0


RETRY = [False]
SCALE = [1.0]


def build(name, rows):
    need = []
    for r in rows:
        for kind, v in columns_of(r):
            if kind == "COL":
                need += [f for f, _l in v]
            elif kind == "GRID":
                need += [x for rr in v[0] for x in rr]
            else:
                need.append(v[0])
    haspdf = {b: fetch(b) for b in need}
    MG = {b: margins(os.path.join(PD, b + ".svg")) for b in need}

    place, labs, y = [], [], 0.0
    for _ri, row in enumerate(rows):
        n0 = len(place)
        for (_f, _r), _v in ROW_DY.items():
            if _f in name and _r == _ri:
                y += _v          # pulled up from that row onwards
        cols = columns_of(row)
        # -- the list of cells in this row; a grid divides further within them
        flat = []            # [(file, cell index)]
        for ci, (kind, v) in enumerate(cols):
            if kind == "COL":
                flat += [(f, ci) for f, _l in v]
            elif kind == "GRID":
                flat += [(x, ci) for rr in v[0] for x in rr]
            else:
                flat.append((v[0], ci))
        FH = max(MG[f][5] for f, _ci in flat)
        topmax = max(MG[f][2] for f, _ci in flat)
        botmax = max(MG[f][3] for f, _ci in flat)
        fy = y + topmax                      # the top y of the frames in this row, which are all equal
        # The gap is not fixed but computed for each row, so that the rightmost frame meets X_R.
        # (With the frames fixed the margins differ slightly from the planning stage, so the planned gap does not fit.)
        seq = []                              # [(left margin, frame width, right margin)]
        for ci, (kind, v) in enumerate(cols):
            if kind == "P":
                f = v[0]
                seq.append((MG[f][0], MG[f][4], MG[f][1]))
            elif kind == "COL":
                fs = [f for f, _l in v]
                seq.append((max(MG[f][0] for f in fs), max(MG[f][4] for f in fs),
                            max(MG[f][1] for f in fs)))
            else:
                grid = v[0]
                ncol = max(len(rr) for rr in grid)
                for j in range(ncol):
                    fs = [rr[j] for rr in grid if j < len(rr)]
                    seq.append((max(MG[f][0] for f in fs), max(MG[f][4] for f in fs),
                                max(MG[f][1] for f in fs)))
        _XR = X_R - XR_ADJ[0]
        need_w = sum(w for _l, w, _r in seq) + sum(seq[j - 1][2] + seq[j][0]
                                                 for j in range(1, len(seq)))
        GAP = ((_XR - X_L) - need_w) / (len(seq) - 1) if len(seq) > 1 else 0.0
        # The left x of the frame is set in advance for each column, because the next cell's left margin must be added in
        # for the rightmost frame to meet X_R.
        XS, cur = [], X_L
        for j, (L, w, R) in enumerate(seq):
            XS.append(cur)
            if j + 1 < len(seq):
                cur = cur + w + R + GAP + seq[j + 1][0]
        # of the cells stacked vertically within one row, as COL or a grid,
        # the nth row must share the top y of its frame, so that the lines below align as well
        # (as for the bottoms of C and B in S4).
        nsub = 1
        for kind, v in cols:
            if kind == "COL":
                nsub = max(nsub, len(v))
            elif kind == "GRID":
                nsub = max(nsub, len(v[0]))
        SUBY, yy0 = [], fy
        for j in range(nsub):
            SUBY.append(yy0)
            hmax = 0.0
            for kind, v in cols:
                fs = []
                if kind == "COL" and j < len(v):
                    fs = [v[j][0]]
                elif kind == "GRID" and j < len(v[0]):
                    fs = list(v[0][j])
                elif kind == "P" and j == 0:
                    fs = [v[0]]
                for f in fs:
                    hmax = max(hmax, MG[f][5] + MG[f][3] + MG[f][2])
            yy0 += hmax + GAP_Y * 0.5
        jj = 0
        for ci, (kind, v) in enumerate(cols):
            if kind == "P":
                f, lab = v
                L, R, T, B, fw, fh = MG[f]
                x = XS[jj] + _cell_dx(name, lab); jj += 1
                place.append((f, x - L, fy - T))
                labs.append((lab, x + LAB_DX + _extra(name, lab), fy + LAB_DY))
            elif kind == "COL":
                x = XS[jj]; jj += 1
                for j, (f, lab) in enumerate(v):
                    L, R, T, B, fw, fh = MG[f]
                    place.append((f, x - L, SUBY[j] - T))
                    labs.append((lab, x + LAB_DX + _extra(name, lab), SUBY[j] + LAB_DY))
            else:
                grid, lab = v
                ncol = max(len(rr) for rr in grid)
                xs_g = XS[jj:jj + ncol]; jj += ncol
                first = True
                for j, rr in enumerate(grid):
                    for j2, f in enumerate(rr):
                        L, R, T, B, fw, fh = MG[f]
                        place.append((f, xs_g[j2] - L, SUBY[j] - T))
                        if first:
                            labs.append((lab, xs_g[0] + LAB_DX + _extra(name, lab),
                                         SUBY[j] + LAB_DY))
                            first = False
        # The y of the next row is taken from the actual bottom of the boxes placed in this row
        # rather than from an estimate, so that the next row cannot ride up over it.
        bot = 0.0
        for f, px, py in place[n0:]:
            _w, _h = size_pt(os.path.join(PD, f + ".svg"))
            bot = max(bot, py + _h)
        y = bot + GAP_Y
    H = y - GAP_Y

    # Using 423.5 as the page width unconditionally clips any cell wider than that
    # (S5 at 1,116 pt).
    # So where the content is wider, that width is used: better to record it as out of specification than to ship a clipped figure.
    WNEED = 0.0
    for f, px, py in place:
        WNEED = max(WNEED, px + size_pt(os.path.join(PD, f + ".svg"))[0])
    # Frame fitting also pushes a cell to the left, giving px below 0.
    # Scaling the whole figure by looking only at the right edge (WNEED) would leave whatever runs off the left clipped
    # (cell D1 of S5 sits at -58.9 pt, where its "spacer-1 + ..." label begins).
    # So whatever runs off the left and the top is included in the width and the scale, and the content shifted in by that much.
    XMIN = min([px for f, px, py in place] + [lx for l, lx, ly in labs] + [0.0])
    YMIN = min([py for f, px, py in place] + [0.0])
    # a figure flagged for whole-figure scaling is scaled here
    _sc = None
    for k, v in FIG_SCALE.items():
        if k in name:
            _sc = v if v else (W_TOTAL - 1.0) / (WNEED - XMIN)
    if _sc and _sc < 0.999:
        place = [(f, (px - XMIN) * _sc, (py - YMIN) * _sc) for f, px, py in place]
        labs = [(l, (lx - XMIN) * _sc, (ly - YMIN) * _sc) for l, lx, ly in labs]
        H = (H - YMIN) * _sc
        WNEED = (WNEED - XMIN) * _sc
        SCALE[0] = _sc
        if XMIN < -0.05 or YMIN < -0.05:
            print("     overhang of left %.1f, top %.1f pt pushed back inside"
                  % (-XMIN * _sc, -YMIN * _sc))
        # the cell letter at 8.2 pt is not scaled, so in a scaled-down figure it protrudes above the page
        # (A of S5 at y0 = -4.9), so the content is pushed down by the height of that ink.
        _top = min([ly - FS_LAB * 1.05 for l, lx, ly in labs]
                   + [py for f, px, py in place] + [0.0])
        if _top < -0.05:
            _d = -_top + 1.0        # +1 is the share of the title that extends beyond the panel's own box
            place = [(f, px, py + _d) for f, px, py in place]
            labs = [(l, lx, ly + _d) for l, lx, ly in labs]
            H += _d
            print("     cell letter is %.1f pt above the page top, lowered by %.1f pt" % (-_top, _d))
        print("     scaled down as a whole by a factor of %.4f (panel that cannot be redrawn): the font-size floor of 4.5 is broken"
              % _sc)
    PAGEW = W_TOTAL
    if WNEED > W_TOTAL + 0.5 and not RETRY[0]:
        return WNEED            # the caller pulls in the right line and calls again
    if WNEED > W_TOTAL + 0.5:
        PAGEW = round(WNEED + 1.0, 2)
        print("     cells wider than the page %.1f > %.1f -> page set to %.1f (outside the spec, to avoid clipping)"
              % (WNEED, W_TOTAL, PAGEW))
    out, idx = [], 0
    for f, px, py in place:
        idx += 1
        ow, oh = size_pt(os.path.join(PD, f + ".svg"))
        out.append('<g transform="translate(%.4f,%.4f)">' % (px, py))
        out.append('<svg width="%.4f" height="%.4f" viewBox="0 0 %.4f %.4f" '
                   'overflow="visible">' % (ow * SCALE[0], oh * SCALE[0], ow, oh))
        out.append(inner(os.path.join(PD, f + ".svg"), "r%d" % idx))
        out.append("</svg></g>")
    for lab, lx, ly in labs:
        out.append('<text x="%.4f" y="%.4f" font-family="Arial" font-size="%.1f" '
                   'font-weight="bold" fill="#000000">%s</text>'
                   % (max(0.5, lx), ly, FS_LAB, lab))
    txt = ('<?xml version="1.0" encoding="utf-8"?>' + NL
           + '<svg xmlns="http://www.w3.org/2000/svg" '
             'xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" '
             'width="%.2fpt" height="%.2fpt" viewBox="0 0 %.2f %.2f">'
             % (PAGEW, H, PAGEW, H)
           + NL + '<rect width="%.2f" height="%.2f" fill="#ffffff"/>' % (PAGEW, H)
           + NL + NL.join(out) + NL + "</svg>" + NL)
    io.open(os.path.join(OUTD, name + ".svg"), "w", encoding="utf8",
            newline=NL).write(txt)
    # A single mismatched tag in a panel's svg
    # breaks the assembled XML, and the renderer stops at the first </svg>, which makes the cell letters disappear.
    # The png can still look intact, so the file is always parsed and the result printed.
    import xml.etree.ElementTree as _ET
    try:
        _ET.fromstring(txt.encode("utf8"))
    except Exception as _e:
        print("     assembly svg is broken as XML: %s -> check the <g> pairs of the panel svg" % _e)
    print("%-46s %7.2f x %7.2f pt, cells %2d" % (name, PAGEW, H, len(labs)))
    # -- measurements for the report
    if not all(haspdf[b] for b in need):
        print("     panel without pdf: png and pdf skipped")
        return
    import pymupdf
    doc = pymupdf.open()
    page = doc.new_page(width=PAGEW, height=H)
    page.draw_rect(pymupdf.Rect(0, 0, PAGEW, H), color=None, fill=(1, 1, 1))
    for f, px, py in place:
        ow, oh = size_pt(os.path.join(PD, f + ".svg"))
        src = pymupdf.open(os.path.join(PD, f + ".pdf"))
        page.show_pdf_page(pymupdf.Rect(px, py, px + ow * SCALE[0],
                                        py + oh * SCALE[0]), src, 0)
        src.close()
    for lab, lx, ly in labs:
        page.insert_text(pymupdf.Point(max(0.5, lx), ly), lab,
                         fontname="hebo", fontsize=FS_LAB, color=(0, 0, 0))
    doc.save(os.path.join(OUTD, name + ".pdf"))
    pg = pymupdf.open(os.path.join(OUTD, name + ".pdf"))[0]
    # the png is 300 dpi, the NAR raster floor
    pg.get_pixmap(matrix=pymupdf.Matrix(300 / 72.0, 300 / 72.0)).save(
        os.path.join(OUTD, name + ".png"))


ONLY = [x for x in os.environ.get("ONLY_FIG", "").split(",") if x]
for nm, R in (("Fig4_main_2026-10-16", FIG4),
              ("SuppFig_S4_blocks_and_molecules_2026-10-16", S41),
              ("SuppFig_S5_dispersed_5S_2026-10-16", S42),
              ("SuppFig_S6_other_macrosatellites_2026-10-16", S44)):
    if ONLY and not any(k in nm for k in ONLY):
        continue
    XR_ADJ[0] = 0.0
    RETRY[0] = False
    SCALE[0] = 1.0
    for _try in range(4):
        _over = build(nm, R)
        if not _over:
            break
        _d = _over - W_TOTAL + 0.5
        XR_ADJ[0] += _d
        print("     right frame line pulled in by %.1f pt (pass %d)" % (XR_ADJ[0], _try + 1))
        if _try == 2:
            RETRY[0] = True      # as a last resort the page is widened so that nothing is clipped
