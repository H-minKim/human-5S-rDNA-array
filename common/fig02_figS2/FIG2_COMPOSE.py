# -*- coding: utf-8 -*-
r"""Assemble the combined Fig 2 exactly as the layout table (canonical) and the 9 rows set by the user.
   - Cut only the panel regions out of the published SVGs and insert them as nested <svg viewBox> (no redrawing)
   - Width fixed at 414.4 pt, height free, panel heights matched per row with generous gaps between cells
   - Cell letters (a to v) are stamped fresh for this layout
   Output goes to the C: scratch folder. Installation is separate."""
import io, os, sys, re, math
import xml.etree.ElementTree as ET
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import SVGPANEL as S

BORROW = {"ko": (__import__("os").environ["RDNA_DATA_ROOT"] + r"/v11_readvar/03_figures/00_figure_set_to_table_mapping_2026-10-08/figures_ko"),
          "en": (__import__("os").environ["RDNA_DATA_ROOT"] + r"/v11_readvar/03_figures/00_figure_set_to_table_mapping_2026-10-08/figures_en")}


def resolve(figdir, name):
    """If it is not in the publish folder, look in the borrowed read folder (the English edition pairs with _EN)."""
    p = os.path.join(figdir, name)
    if os.path.exists(p): return p
    en = figdir.endswith("영어판")
    base, ext = os.path.splitext(name)
    cand = (base + "_EN" + ext) if en else name
    q = os.path.join(BORROW["en" if en else "ko"], cand)
    return q if os.path.exists(q) else p

SVG = S.SVG
W = 414.4                 # pt, typeset coordinate width
# spec figure Untitled-1.svg (= this composite x 1.022), published width 423.5 pt, all figures share one horizontal width
#   Only the outer <svg> width and height are enlarged and the viewBox is kept, so text, line widths, markers and margins all grow by the same factor (lines never fail to grow)
OUT_W = 423.5
SCALE = OUT_W / W
MARGIN, GAP_X, GAP_Y = 6.0, 10.0, 5.0     # reduce the vertical margin between rows (was 16)
LAB = 6.0                 # Space for cell letter height (2026-10-18 night: 9 -> 6 when letters were moved down to the graph frame)
# panel letters are anchored at the top-left corner of the graph frame (top end of the y axis), not the top-left of the cell
#   Letter baseline = top line of that row's graph frame - LET_DY, letter left edge = left line of the leftmost frame in that cell - LET_DX (clamped inside the cell left edge)
#   => the leftmost letter of each row sits on the same vertical line (since the frame left lines coincide), and letters in a row share the same height
LET_DX, LET_DY = 22.0, 1.5

GAP_COL = 12.0            # Between the two columns

# Main figure, two columns at the top and one full-width row at the bottom
# zero gap between cells, the gap is provided by the cell figure margins (cellstyle.ROWS), graph frames align at the top, bottom, left and right ends of each row
MAIN_BANDS = {"bands": [{"full": [["a", ("b", 1.0, 0.0)]]},
                        {"full": [["c", ("d", 1.0, 0.0), ("e", 1.0, 0.0)]]},
                        {"full": [["f", ("g", 1.0, 0.0)]]},          # old g and h merged into one cell g, later letters shift up by one (old i->h ... s->r)
                        {"full": [["h", ("i", 1.0, 0.0), ("j", 1.0, 0.0), ("k", 1.0, 0.0)]]},
                        {"full": [["l", ("m", 1.0, 0.0), ("n", 1.0, 0.0)]], "dy": -11.0},    # 'move l, m, n up a little' (-6), 2026-10-28 user 'move l,m,n up a bit': empty band above the row is 20.2 (EN) and 18.7 pt (KO) (other row gaps 13 to 17) -> 5 more (-11)
                        {"full": [["o", ("p", 1.0, 0.0), ("q", 1.0, 0.0), ("r", 1.0, 0.0)]], "dy": -5.0}]}   # 'move o, p, q, r up a little too' (5 on top of the 6 for the l row)

# Supplement, 2026-10-16: c, d, f, h and old e go to the supplement, e = fig_coshare, remaining cell letters unchanged (empty slots)
# the supplement also uses cell figures (cellS2_*, MK_SUPP_SPEC), zero gap between cells, graph frames align at the top, bottom, left and right ends of each row
# supplement cell letters run A to P in reading order; rows and positions unchanged
# (7 rows -> 6): 'narrow a,b and put c on the top row, d,e on the same row, f,g on the same row, narrow i so it shares a row with h', cell letters unchanged
#   Old rows (earlier that day): a b / c d (dy -9) / e f / g h / i (dy -7) / j k l / m n o
SUPP_BANDS = {"bands": [{"full": [["a", ("b", 1.0, 0.0), ("c", 1.0, 0.0)]]},
                        {"full": [["d", ("e", 1.0, 0.0)]], "dy": -9.0},     # Top margin 28.5 (slanted header of d), same pull as the old c, d row
                        {"full": [["f", ("g", 1.0, 0.0)]]},
                        {"full": [["h", ("i", 1.0, 0.0)]]},
                        {"full": [["j", ("k", 1.0, 0.0), ("l", 1.0, 0.0)]]},
                        {"full": [["m", ("n", 1.0, 0.0), ("o", 1.0, 0.0)]]}]}

# 'move opqr of main above jkl of supple and split supple in half at h (it is one figure but too long, so split it in two)'
#   + user answers 'second page from H' and 'continue the lettering' -> last main row (o p q r) moves to the supplement, supplement = one figure on two pages
#   Main: MAIN_BANDS above minus its last row (A to N)
#   Supplement page 1: A B C / D E / F G, page 2: H I / J K L M (old main O P Q R) / N O P (old supplement J K L) / Q R S (old supplement M N O)
#   The mapping from cell letter to cell figure is set by the composition table (FIG2 composition table in the tables folder), supplement j..o -> n..s, main o..r -> supplement j..m
#   Same night, user 'sorry, h,i should stay in one figure, split from j' -> page 1 A to I (four rows), page 2 J to S (three rows)
#   Same night, user 'reduce the vertical spacing between panels in supple' -> shrink the empty bands between rows (measured on the composite png, EN/KO)
#     Page 1: C->D 7.7/10.8, E->F 19.4/19.0, G->H 9.8/10.1 pt, page 2: M->N 19.2/10.1, P->Q 8.9/8.6 pt (before pulling)
MAIN_BANDS_1027S = {"bands": MAIN_BANDS["bands"][:-1]}
SUPP1_BANDS = {"bands": [SUPP_BANDS["bands"][0],
                         {"full": [["d", ("e", 1.0, 0.0)]], "dy": -30.0},   # 3 more on top of -9, same night user 'pull d,e up a bit more': after removing the A n text, empty band above the D header is 27.4 (EN) and 34.1 (KO), 34.8 pt below the C legend -> 18 more (-30)
                         {"full": [["f", ("g", 1.0, 0.0)]], "dy": -10.0},
                         {"full": [["h", ("i", 1.0, 0.0)]], "dy": -4.0}]}
SUPP2_BANDS = {"bands": [{"full": [["j", ("k", 1.0, 0.0), ("l", 1.0, 0.0), ("m", 1.0, 0.0)]]},
                         {"full": [["n", ("o", 1.0, 0.0), ("p", 1.0, 0.0)]], "dy": -22.0},   # Old bottom margin of the main row + top margin 11 of this row -> empty band 28.3 pt (EN), 5 more on top of -9, M ticks went from three lines to one, so 8 more
                         {"full": [["q", ("r", 1.0, 0.0), ("s", 1.0, 0.0)]], "dy": -3.0}]}

# user layout instruction: main 'remove graph e and put f in its place, the two plots of g on the left, h, i on the right, bottom row j left, k right, f of FigS2a in the middle'
# cell E is not included; cell letters run in reading order and FigS2b runs I to R
#   Main: A B / C D E(old F) / F(old G) G(old H) H(old I) / I(old J) J(old K) K(old supplement F) / L M N  (user 'order is i, k, then j'), supplement page 1: A B C / D E / F(old G) G(old H) H(old I), page 2 I to R (old J to S)
#   The mapping from cell letter to cell figure is set by the composition table (copy of MK_TAB_1028L), cell sizes by cellstyle.ROWS and MK_SUPP_SPEC (middle of main row 4 = cellS2_k)
MAIN_BANDS_1028L = {"bands": [MAIN_BANDS["bands"][0],
                              {"full": [["c", ("d", 1.0, 0.0), ("e", 1.0, 0.0)]]},
                              {"full": [["f", ("g", 1.0, 0.0), ("h", 1.0, 0.0)]]},
                              {"full": [["i", ("j", 1.0, 0.0), ("k", 1.0, 0.0)]], "dy": -14.0},   # Measured on preview: empty band below H ticks about 37 pt
                              MAIN_BANDS["bands"][4]]}
# Same day, user 2nd instruction (supplement) 'narrow d,e, put f above (order f,d,e), and g,h back to the original figures' -> reading order: D(old F split gain) E(old D rank table) F(old E) / G H (as published, pull -4 of old row h i)
#   Same day, user 3rd instruction 'move f before d' -> second row D(published E, length proportion) E(published G, split gain) F(published D, rank table)
#   Same day, user 4th instruction 'widen e, d to the left of g, d,g,h on one row' -> D(published G, split gain) E(published D, rank table) / F(published E, length proportion) G(published H) H(published I)
#   Same day, user 5th instruction 'just drop d, gh back to the original' -> length proportion (published E) removed, page 1 A B C / D E / F(published H) G(published I), page 2 H to Q
SUPP1_BANDS_1028L = {"bands": [SUPP1_BANDS["bands"][0],
                               {"full": [["d", ("e", 1.0, 0.0)]], "dy": -12.0},   # 4th: D(split gain) E(rank table), so the slanted header of E does not touch the legend below C (at -30 it touched)
                               {"full": [["f", ("g", 1.0, 0.0)]], "dy": -14.0}]}   # 8th, user '(S2a) move f,g up a bit': -4 -> -14   # 5th 'drop d, gh as original': F(published H boundary class) G(published I fingerprint), published row, pull -4 also as published
# 5th: supplement page 1 has 7 cells (A to G) -> page 2 H to Q (lettering continued, rows and pulls unchanged)
SUPP2_BANDS_1028L = {"bands": [{"full": [["h", ("i", 1.0, 0.0), ("j", 1.0, 0.0), ("k", 1.0, 0.0)]]},
                               {"full": [["l", ("m", 1.0, 0.0), ("n", 1.0, 0.0)]], "dy": -24.0},   # Measured on preview: empty band below the J axis label about 36 pt
                               {"full": [["o", ("p", 1.0, 0.0), ("q", 1.0, 0.0)]], "dy": -16.0}]}   # Measured on preview: empty band below M ticks about 24 pt
# 6th (FigS2b): H(terminal types by depth) I(diagnostic variants) J(gap distribution) K(gap fold) / L(exchange) M(neighbour difference) N(Jaccard) / O(carriage by population) P(depth sets) Q(cell lines), rows and letter positions unchanged, only the cell figures change (composition table)

SPEC = {"Fig2": MAIN_BANDS, "SFig2": SUPP_BANDS,
        "Fig2-1027s": MAIN_BANDS_1027S, "SFig2-1": SUPP1_BANDS, "SFig2-2": SUPP2_BANDS,
        "Fig2-1028L": MAIN_BANDS_1028L, "SFig2-1L": SUPP1_BANDS_1028L, "SFig2-2L": SUPP2_BANDS_1028L}


def spec_keys(spec):
    """Cell letters placed in that layout (band and row order)"""
    return [x[0] if isinstance(x, tuple) else x for band in spec["bands"] for r in band["full"] for x in r]


def _k(x): return x[0] if isinstance(x, tuple) else x
def _w(r): return max([x[1] for x in r if isinstance(x, tuple)] + [1.0])
def _g(x): return x[2] if isinstance(x, tuple) and len(x) > 2 else None    # separate leading gap (pt) per cell (GAP_X if absent)


def layout(spec, keys, aspect, y0=None):
    """Compute cell positions only (no drawing), shared by the SVG composer and the PNG builder.
       spec : list = one column, 2-tuple = two columns (stacked separately), {"bands": [(left, right), ...]} = heights matched per band
       returns : ([(k, x, y, w, h)], bottom y)"""
    y0 = MARGIN if y0 is None else y0
    keys = set(keys)
    cells = []

    def rows_of(sp):
        rr = [([_k(x) for x in r if _k(x) in keys], _w(r), [_g(x) for x in r if _k(x) in keys]) for r in sp]
        return [(r, w, g) for r, w, g in rr if r]

    def lay(rows, x0, availw, ytop):
        y = ytop
        for r, rw, gg in rows:
            asp = [aspect(k) for k in r]
            gaps = [0.0] + [GAP_X if g is None else g for g in gg[1:]]
            h = (availw * rw - sum(gaps)) / sum(asp)
            x = x0
            for k, a, gp in zip(r, asp, gaps):
                x += gp
                cells.append((k, x, y + LAB, h * a, h))
                x += h * a
            y += LAB + h + GAP_Y
        return y

    cw = (W - 2 * MARGIN - GAP_COL) / 2.0
    xr = MARGIN + cw + GAP_COL
    if isinstance(spec, dict):
        y = y0
        for band in spec["bands"]:
            if isinstance(band, dict):                       # One full-width band
                y += band.get("dy", 0.0)                     # reduce the gap before the band to move it up (cell figures unchanged, the bands below move up too)
                y = lay(rows_of(band["full"]), MARGIN, W - 2 * MARGIN, y)
            else:
                L, R = band
                y = max(lay(rows_of(L), MARGIN, cw, y), lay(rows_of(R), xr, cw, y))
    elif isinstance(spec, tuple):
        y = max(lay(rows_of(spec[0]), MARGIN, cw, y0), lay(rows_of(spec[1]), xr, cw, y0))
    else:
        y = lay(rows_of(spec), MARGIN, W - 2 * MARGIN, y0)
    got = [c[0] for c in cells]
    assert len(got) == len(set(got)), "한 칸이 두 번 놓였다 : %s" % got
    miss = sorted(k for k in keys if k not in set(got))
    assert not miss, "줄 묶음에서 빠진 칸 : %s" % miss
    return cells, y

ROWS = [["a", "b"],                       # Variants by position, variant load by array position
        ["c", "d", "e"],                  # Homogenization: range (small), distance dependence, exchange unit
        ["f", "g", "h"],                  # Composition by origin
        ["i"],                            # Distal CpG skew deficit (restored 2026-10-11), a wide cell so it takes a whole row
        ["j", "k", "l", "m"],             # Region, length
        [("n", 0.60)],                    # Segmentation, width 60%
        ["o", "p", "q", "r"],             # Terminal types, carriage by population
        ["s", "t", "u"]]                  # Array types, population-specific tokens

# The supplement follows the same flow as the main figure (2026-10-10)
ROWS_SA = [["a", "b"], ["c"], ["d"], ["e", "f", "g"], ["h", "i"], ["j", "k"], ["l", "m"]]
ROWS_SB = [["a"], ["b"], ["c", "d"], ["e", "f"], ["g"], ["h"], ["i", "j"], ["k", "l"], ["m"]]   # b = Fig9 (B) added 2026-10-11, later letters shift by one
ROWS_SC = [["a"], ["b", "c"], ["d", "e"], ["f", "g"], ["h", "i"], ["j", "k"],
           ["l", "m"], ["n", "o"], ["p", "q"]]


def load_plan(tab, which="Fig2"):
    P, ci = {}, None
    for l in io.open(tab, encoding="utf8"):
        if l.startswith("#") or "\t" not in l: continue
        f = l.rstrip("\n").split("\t")
        if f[0] == "그림":
            ci = f.index("본문칸그림") if "본문칸그림" in f else None     # cell figure column for the main figure
            continue
        if f[0] != which: continue
        if ci is not None and len(f) > ci and f[ci].strip():
            P[f[1]] = dict(파일=f[ci].strip(), 패널="(칸)", 제목=f[4])   # Use the whole cell figure (source figure and panel stay in the source file and source panel columns)
        else:
            P[f[1]] = dict(파일=f[2], 패널=f[3], 제목=f[4])
    return P


def cell_svg(figdir, plan, key):
    """(nested svg Element, aspect ratio), or (None, 1.4) if absent"""
    p = resolve(figdir, plan[key]["파일"])
    if not os.path.exists(p): return None, 1.4
    if plan[key]["패널"] == "(칸)":                 # main cell figure, figure size = cell size, used whole without cutting anything
        root = S.parse(p)
        box = S.viewbox(root)
        g = next((e for e in root.iter("{%s}g" % SVG) if (e.get("id") or "").startswith("figure_")), None)
        assert g is not None, p
        # make the full white background of the cell figure (first patch of the figure) transparent, when a row is pulled up (dy) and the cell box overlaps the bottom margin of the cell above
        #   the white background covers the text of the upper cell (raising the D and E row of supplement a made the legend below C disappear, a text content check does not catch it because the text is still in the svg)
        #   The white background inside axes (axes patch) is kept
        _bg = next((c for c in g if (c.get("id") or "").startswith("patch_")), None)
        if _bg is not None:
            for _pth in _bg.iter("{%s}path" % SVG): _pth.set("style", (_pth.get("style") or "").replace("fill: #ffffff", "fill: none"))
        sub = ET.Element("{%s}svg" % SVG)
        sub.set("viewBox", "%.4f %.4f %.4f %.4f" % box)
        dfs = S.needed_defs(root, g)          # Deferred fix: the order of needed_defs varies between runs (see the SVGPANEL.needed_defs comment), so a differing md5 does not prove a different figure
        if dfs:
            d = ET.SubElement(sub, "{%s}defs" % SVG)
            for e in dfs: d.append(e)
        sub.append(g)
        S.prefix_ids(sub, "c%s" % key)
        return sub, box[2] / box[3]
    root, (nm, g, box, ax, tx) = S.find_panel(p, plan[key]["제목"], plan[key]["패널"])
    if nm == "figure":                              # Cells that use the whole figure
        box = S.clean_figure_cell(root, g, box)
    elif S.DROP_TITLE and S.REGION in ("ink", "mid+ink"):
        S._drop_title_texts(g, ax)                  # Remove the original panel title (its letter), the composite adds its own letters
        S.drop_label_title(g, plan[key]["패널"])     # Table figure with axes off, title starting with the original letter
    else:
        S.strip_panel_letter(g, ax, plan[key]["패널"])
    pre = "c%s" % key
    sub = ET.Element("{%s}svg" % SVG)
    sub.set("viewBox", "%.4f %.4f %.4f %.4f" % box)
    dfs = S.needed_defs(root, g)          # Bring along the clipPaths and markers that the panel uses
    if dfs:
        d = ET.SubElement(sub, "{%s}defs" % SVG)
        for e in dfs: d.append(e)
    sub.append(g)
    S.prefix_ids(sub, pre)                # defs and references together, avoids id collisions between cells
    return sub, (box[2] / box[3] if box[3] else 1.4)


_PATCH = re.compile(r"M\s*([-\d.]+)\s+([-\d.]+)\s*L\s*([-\d.]+)\s+([-\d.]+)\s*L\s*([-\d.]+)\s+([-\d.]+)\s*L\s*([-\d.]+)\s+([-\d.]+)\s*z")


def frame_tl(sub, x, w):
    """(left x, top y) of the leftmost frame among the axes background rectangles (axes_N > patch_M) in the cell svg, in composite pt relative to the cell top-left, or None"""
    if sub is None: return None
    vb = [float(v) for v in sub.get("viewBox").split()]
    s = w / vb[2]
    best = None
    for g in sub.iter("{%s}g" % SVG):
        if not re.search(r"axes_\d+$", g.get("id") or ""): continue
        for c in g:
            if re.search(r"patch_\d+$", c.get("id") or ""):
                pth = c.find("{%s}path" % SVG)
                m = _PATCH.search(pth.get("d") or "") if pth is not None else None
                if m:
                    v = list(map(float, m.groups())); xs, ys = v[0::2], v[1::2]
                    cand = ((min(xs) - vb[0]) * s, (min(ys) - vb[1]) * s)
                    if best is None or cand[0] < best[0] - 0.01: best = cand
                break
    return best


def compose(figdir, tab, out, which="Fig2", title=None):
    """Typeset according to SPEC[which], positions are set by layout() alone (same positions as the PNG builder)."""
    plan = load_plan(tab, which.split("-")[0])
    if "-" in which:                                 # one page of a figure or part of the main figure, only the cells placed in that layout
        ks = spec_keys(SPEC[which])
        assert set(ks) <= set(plan), ("구성표에 없는 칸", sorted(set(ks) - set(plan)))
        plan = {k: plan[k] for k in ks}
    subs, asp, miss = {}, {}, []
    for k in plan:
        s, a = cell_svg(figdir, plan, k)
        subs[k], asp[k] = s, a
        if s is None: miss.append(k)
    y0 = MARGIN + (14.0 if title else 0.0)
    cells, y = layout(SPEC[which], list(plan), lambda k: asp[k], y0)
    doc = ET.Element("{%s}svg" % SVG, {"width": "%.2fpt" % W, "version": "1.1"})
    body = ET.SubElement(doc, "{%s}g" % SVG)
    if title:
        t = ET.SubElement(body, "{%s}text" % SVG,
                          {"x": "%.2f" % MARGIN, "y": "%.2f" % (MARGIN + 9),
                           "style": "font-family: Helvetica, Arial; font-size: 10px"})
        t.text = title
    built, LABS = [], []
    FR = {k: frame_tl(subs[k], x, w) for k, x, yy, w, h in cells}      # (frame left x, frame top y) cell coordinates -> composite coordinates
    ROWTOP = {}
    for k, x, yy, w, h in cells:
        if FR[k] is not None:
            ROWTOP[round(yy, 3)] = min(ROWTOP.get(round(yy, 3), 1e9), yy + FR[k][1])
    for k, x, yy, w, h in cells:
        fx = (x + FR[k][0]) if FR[k] is not None else None
        lx = max(x + 0.5, fx - LET_DX) if fx is not None else x + 0.5
        if x <= MARGIN + 0.01: lx = min(lx, MARGIN + 35.0 - LET_DX)   # leftmost cell letters of a row share one vertical line (based on cell left margin LM 35), when the supplement D table became the first cell of a row its letter sat 28 pt inside, measured from the table cell left edge
        ly = ROWTOP.get(round(yy, 3), yy + 7.0) - LET_DY
        lab = ET.Element("{%s}text" % SVG,             # Attach after all cells are placed (so the cell white background does not cover the letters)
                         {"x": "%.2f" % lx, "y": "%.2f" % ly,
                             "style": "font-family: Helvetica, Arial; font-size: 8px; font-weight: bold"})   # Helvetica, no bold -> same night 'make the panel numbers (A,B,C...) bold helvetica' (only panel letters bold, other text stays regular)
        lab.text = k.upper()                        # panel letters are capitals A B C ...
        LABS.append(lab)
        s = subs[k]
        if s is None:
            ET.SubElement(body, "{%s}rect" % SVG,
                          {"x": "%.2f" % x, "y": "%.2f" % yy, "width": "%.2f" % w,
                           "height": "%.2f" % h,
                           "style": "fill:none; stroke:#c0392b; stroke-dasharray:4,3; stroke-width:0.8"})
            tt = ET.SubElement(body, "{%s}text" % SVG,
                               {"x": "%.2f" % (x + w / 2), "y": "%.2f" % (yy + h / 2),
                                "text-anchor": "middle",
                                "style": "font-family: Helvetica, Arial, Malgun Gothic; font-size: 7px; fill:#c0392b"})
            tt.text = "%s — 빌려온 그림 (파일 없음)" % k
        else:
            # width and height with four decimals like the viewBox (was %.2f), a cell whose width was truncated downward (117.3247 -> 117.32) gets a meet scale below 1
            #   and in the Chrome pdf that cell's letters printed at 4.793 pt (caught in the user's 'font floor 4.8' measurement, main C and J, supplement J and K)
            s.set("x", "%.4f" % x); s.set("y", "%.4f" % yy)
            s.set("width", "%.4f" % w); s.set("height", "%.4f" % h)
            s.set("preserveAspectRatio", "xMidYMid meet")
            body.append(s)
        built.append((k, x, yy, w, h))
    for lab in LABS: body.append(lab)
    H = y - GAP_Y + MARGIN
    doc.set("width", "%.2fpt" % (W * SCALE))
    # height rounded UP to 0.01 pt (was %.2f rounding), when it rounds down the meet scale is set by the height and is 1e-5 smaller than the width scale,
    #   and since Chrome printing truncates font size to 0.01 px, 6.40 -> 6.39 px (published 4.7925 pt, below the 4.8 floor) (measured on the FigS2a Korean edition, 445.5742 -> 445.57)
    doc.set("height", "%.2fpt" % (math.ceil(H * SCALE * 100.0 - 1e-9) / 100.0))
    doc.set("viewBox", "0 0 %.2f %.2f" % (W, H))
    ET.ElementTree(doc).write(out, encoding="utf-8", xml_declaration=True)
    return built, miss, W * SCALE, H * SCALE


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    P = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout")
    TAB = os.path.join(P, "tables", "FIG2_layout_2026-10-09.tsv")
    HERE = os.path.dirname(os.path.abspath(__file__))
    OUT = os.path.join(HERE, "OUT_FIG2")
    os.makedirs(OUT, exist_ok=True)
    JOBS = [("Fig2", "fig2_main"), ("SFig2", "fig2_supp")]
    for which, stem in JOBS:
      for lang, d in (("ko", "그림"), ("en", "그림_영어판")):
        o = os.path.join(OUT, "%s_%s.svg" % (stem, lang))
        built, miss, w, h = compose(os.path.join(P, d), TAB, o, which=which, title=None)
        print("■ %-7s %s · 칸 %d · 빠진 칸 %s · 폭 %.1f pt · 높이 %.1f pt"
              % (which, lang, len(built), miss or "없음", w, h))
        print("   %s · %d B" % (o, os.path.getsize(o)))
        # Overlap check
        bad = 0
        for i in range(len(built)):
            for j in range(i + 1, len(built)):
                _, x1, y1, w1, h1 = built[i]; _, x2, y2, w2, h2 = built[j]
                if x1 < x2 + w2 and x2 < x1 + w1 and y1 < y2 + h2 and y2 < y1 + h1: bad += 1
        print("   칸끼리 겹침 %d 쌍" % bad)
