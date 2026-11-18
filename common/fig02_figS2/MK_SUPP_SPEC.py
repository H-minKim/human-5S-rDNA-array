# -*- coding: utf-8 -*-
r"""Supplementary cell figure design table -> spec for _retarget (Korean and English) (same alignment as the main figure)
   (1) graph frames in the same row share top and bottom (frame height 64 pt) (2) the left line of the leftmost frame is the same in every row (3) the right line of the rightmost frame is the same
   Rightmost margin RMS = 36 pt (measured overflow of the b line label "5S rDNA gene" 6.6 pt -> 28+8), so that the line labels of supplementary b sit outside the right of the graph (same for all rows)
   Output: spec_supp_ko.json, spec_supp_en.json, supp_layout.json (cell box sizes, for layout checking)"""
import json, os, sys
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "_fix"))
W_USE, LM, RMS, FH = 414.4 - 12.0, 35.0, 36.0, 64.0
# cell: (name, code, saved original name, axes, left {ko, en}, right, width ratio, split, options)
# applied: axis titles in the English version on one line (shortened meanings go into the Legend text), in the Korean version on two lines if long. Legends without a box at upper left (upper right if there is no room, o at lower right)
#   j: red explanatory text removed. v y max 1.2, w y max 0.2 (max plotted value v 1.0000, w 0.1765, no data clipped, measured from the value table)
_ETK = ["101–250", "251–500", "501–1,000", "1,001–2,240"]
# rows a b e / g j / k l / m y / o q r / v w x (panel letters keep reading order)
# For the first four rows the width-ratio slot holds the graph frame width (pt) directly (per row, sum = usable width - margins).
#   a is wide enough that the ticks 'Transmitted, Non-transmitted, Unrelated' do not touch (the English version puts Non-transmitted on two lines)
#   e: two panels of 50.7 each, x range -0.3 to 3.3 (so that the ticks '1,000' and '2,240' do not touch)
ROWS = [
  (11.0, (41.5, 38.0), [("a", "FIG_TRIOPOS.py", "fig_trio_unitmatch", [0], (LM, LM), 3.0, 76.0, ("one",),
                        {"labels": {"0": {"y": ["자리까지 같은 유닛 비율", "Fraction identical"]}},   # 'Position-identical fraction' is longer than the frame (64 pt), so 'Fraction identical' (the position meaning goes to the Legend)
                         # only the transmitted value 0.982 is kept: 0.053, 0.000, '436 / 447' and the 'n = ' text below the ticks are removed (n goes to the Legend text, values are in the FIG_TRIOPOS value table)
                         "drop_text": {"0": ["436 / 447", r"re:0\.0\d\d", "n = "]},
                         "xtick_sub": {"0": [["Non-transmitted", "Non-\ntransmitted"]]},
                         # small end padding on the x range widens the tick spacing, so that the English ticks 'Transmitted' and 'Non-transmitted' do not overlap. In the English version the ticks take two lines, so n = is moved down one line
                         "post": "A[0].set_xlim(-0.38, 2.38)\n"
                                 "if EN:\n"
                                 "    for t_ in A[0].texts:\n"
                                 "        if t_.get_text().startswith('n = '): t_.set_y(t_.get_position()[1] - CS.FS_TICK * 1.25 / 64.0)\n"
                                 "else:\n"   # in the Korean version n = is moved down 1 pt to clear the ticks
                                 "    for t_ in A[0].texts:\n"
                                 "        if t_.get_text().startswith('n = '): t_.set_y(t_.get_position()[1] - 1.0 / 64.0)\n"}),
                       ("b", "FIG_G6e.py", "fig_G6e_region_neighbour", [0], (22.0, 22.0), RMS + 3.0, 57.0, ("one",),   # right +3 so that the Korean line label for '5S rDNA gene' stays inside the cell
                        {"repel_out": {"0": 6.0},   # 6.0 gives a line label clearance of 0.53 pt in Helvetica
                          # same quantity as in main C and D (position difference between two units in the same array), so the same axis name
                          "labels": {"0": {"y": ["영역 일치 확률", "Probability of identity"], "x": ["유닛 거리 (copy)", "Separation (repeat units)"]}},   # common term instead of a formula form
                         "grid_off": ["0"]}),   # remove the horizontal grid lines at 0.8, 0.6, 0.4 and 0.2
                       ("e", "FIG_COSHARE.py", "fig_coshare", [0, 1], (21.0, 21.0), RMS, 113.4, ("h2", 12.0),
                         # the cell is narrow, so the bottom legend is one column. The Korean axis titles take two lines, so the row bottom margin is 41.5 (English 38)
                         {"keep_title": True, "fig_legend": {"ncol": 1}, "hide_ylabel": [1],
                          # no lines joining same-colored points: only points and error bars (also no horizontal line in the legend handles)
                          "post": "from matplotlib.lines import Line2D as _L2\n"
                                  "_NOM = (None, '', 'None', ' ')\n"
                                  "for a_ in (A[0], A[1]):\n"
                                  "    for l_ in a_.get_lines():\n"
                                  "        if l_.get_marker() not in _NOM: l_.set_linestyle('None')\n"
                                  "for lg_ in fig.legends:\n"
                                  "    for l_ in lg_.findobj(_L2):\n"
                                  "        if l_.get_marker() in _NOM and l_.get_linestyle() not in ('None', '', ' '): l_.set_visible(False)\n"
                                  "for a_ in (A[0], A[1]): a_.set_xlim(-0.3, 3.3)\n",
                          "labels": {"0": {"y": ["다른 변이 공유 비율", "Frac. shared"],
                                           "x": ["A 에서 그 변이까지의\n거리 (bp)", "Distance from A (bp)"]},   # the panels are narrow, so the Korean version uses two lines
                                     "1": {"x": ["A 에서 그 변이까지의\n거리 (bp)", "Distance from A (bp)"]}},
                          # en dashes in the tick text become a hyphen '-', with the tick split over two lines
                          "xtick_sub": {"0": [[t, t.replace("–", "-\n")] for t in _ETK], "1": [[t, t.replace("–", "-\n")] for t in _ETK]},
                          # value text sits on the connecting line -> orange (upper line) slightly up, gray (lower line) slightly down (pt)
                          #   (gray in panel A sits at 0, moving it down would leave the axes -> up)
                          "shift_pt": {"0": [[r"0\.\d{3}", 1.0, 1.6, "#e67e22", "bottom"], [r"0\.\d{3}", 1.0, 1.6, "#7f8c8d", "bottom"]],
                                       # panel B: only orange values above 0.9 (first point) go below the line (above is the frame top 1.0 and the subtitle). Other orange values above, gray below
                                       "1": [[r"0\.9\d\d", 1.0, -1.6, "#e67e22", "top"], [r"0\.[0-8]\d\d", 1.0, 1.6, "#e67e22", "bottom"],
                                             [r"0\.\d{3}", 1.0, -1.6, "#7f8c8d", "top"]]},
                          # the '×13, ×7, ×6, ×4' text is removed, and the two panels (orange and gray points) share the y-axis range (0 to 1.0)
                          "drop_text": {"0": ["×", r"re:0\.\d{3}"], "1": [r"re:0\.\d{3}"]},   # value text above the points removed (values in the value table)
                          "ylim": {"0": [0.0, 1.0], "1": [0.0, 1.0]}, "yticks": {"0": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], "1": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]}})]),
  # top margin 31.5, so that the slanted headers of supplementary D (such as 'G>A (non-CpG)') stay inside the cell
  (31.5, (20.0, 20.0), [("g","FIG_VTORG2C.py", "fig_VTORG2C_rank_table", [0], (LM, LM), 8.0, 155.0, ("one",),
                         # right 8: room for the end of the slanted header 'A>G'
                         {"table_axes": [0], "retext": {"0": [["proximal 1-3", "prox. 1-3"], ["proximal 4-17", "prox. 4-17"], ["distal 1-10", "dist. 1-10"]]},
                          # the top of the vertical group label goes below the horizontal line between the two groups (position only)
                          "post": "for t_ in A[0].texts:\n"
                                  "    if t_.get_text().strip() in ('단일유닛형만', 'single-unit only'):\n"
                                  "        t_.set_va('top'); t_.set_y(-4 * 0.78 - 0.10)\n"}),   # original CH 0.78, horizontal line y = -4 CH
                       ("j", "FIG_G3b4.py", "fig_G3b_region_character", [5], (LM, LM), RMS, 133.4, ("one",),   # last cell of the row, so the right margin is RMS
                        # lines and markers thin and small like main K (line 0.45, marker 1.1 pt), legend kept off the graph, 'U84C' deleted
                        #   the 'proportional line, fixed line' text goes into a small lower-right legend like main K
                        {"drop_text": {"5": ["#7b241c", "비례선", "고정선", "proportional line", "fixed line"]},
                         "line_abs": {"5": [0.45, 1.1]},
                         "legend": {"5": {"loc": "upper left", "borderaxespad": 0.0, "labelspacing": 0.1,
                                          "text": [["proportional individual variant 2205TC (U84C)", "Proportional 2205TC"],
                                                   ["fixed individual variant 2235TG", "Fixed 2235TG"],
                                                   ["fixed individual variant D1096:1 (terminal-type diagnostic)", "Fixed D1096:1 (terminal diag.)"]],
                                          "text_re": [[r"\s*\(U84C\)", ""]]}},
                         "post": "import matplotlib.lines as _ml\nax_ = A[5]\n_lg = ax_.get_legend(); ax_.add_artist(_lg)\n"
                                 "_d = [l for l in ax_.get_lines() if l.get_linestyle() == ':'][0]\n_f = [l for l in ax_.get_lines() if l.get_linestyle() == '--'][0]\n"
                                 "ax_.legend([_ml.Line2D([], [], color=_d.get_color(), ls=':', lw=0.45), _ml.Line2D([], [], color=_f.get_color(), ls='--', lw=0.45)],\n"
                                 "           ['Proportional line', 'Fixed line'] if EN else ['비례선', '고정선'], loc='lower right', frameon=False, fontsize=CS.FS_LEG,\n"
                                 "           handlelength=1.4, labelspacing=0.1, borderaxespad=0.0, handletextpad=0.3)\n",
                         "labels": {"5": {"y": ["유닛 수 (평균 = 1)", "Units of type (mean = 1)"],
                                          "x": ["array 길이 (평균 = 1)", "Array length (mean = 1)"]}}})]),
  (9.0, (30.0, 30.0), [("k", "FIG_G2b.py", "fig_G2b_zoning_map", [0, 1], (LM, LM), 4.0, 185.0, ("below", 10.0, 18.0),
                        {"labels": {"0": {"y": ["하플로타입", "Haplotype"]}},
                         # the 7-row legend is taller than the axes and covers data -> split the same handles and text in two, into empty space (4 boundary markers = lower right of the top axes, 3 lines = right of the zone strip)
                         "post": "a_, kw_ = A[0].get_legend()._rt_call\nhs_ = kw_['handles']\n"
                                 "for h_ in hs_:\n    h_.set_markersize(h_.get_markersize() * 0.55); h_.set_markeredgewidth(h_.get_markeredgewidth() * 0.55); h_.set_linewidth(h_.get_linewidth() * 0.55)\n"
                                 # with the zone strip at the full cell width, a 'cohort zone boundary' legend on its right would overlap '1-10' and '11-40' -> that row goes at the bottom of the boundary legend of the top axes
                                 "import matplotlib.lines as _ml\n"
                                 "_hb = _ml.Line2D([], [], color='#000000', lw=0.4, label='Cohort zone boundary' if EN else '코호트 구역 경계')\n"
                                 "l1 = A[0].legend(handles=hs_[:4] + [_hb], loc='lower right', frameon=False, fontsize=CS.FS_LEG, labelspacing=0.1, borderaxespad=0.0, handletextpad=0.5)\n"
                                 "relabel(l1, [['First boundary from the proximal end', 'First boundary, proximal'], ['First boundary from the distal end', 'First boundary, distal'],\n"
                                 "             ['Single boundary (counts for both ends)', 'Single boundary (both ends)']])\n"
                                 # cohort zone boundaries (proximal 3|4 and 17|18 vertical lines, distal 10|11 curve) as thin black lines. Drop the gray 'array end (distal junction)' line
                                 "import matplotlib.lines as _ml\n"
                                 "for l_ in list(A[0].get_lines()):\n"
                                 "    c_ = str(l_.get_color()).lower()\n"
                                 "    if c_ in ('#1a5276', '#943126'): l_.set_color('#000000'); l_.set_linewidth(0.4)\n"
                                 "    elif c_ == '#7f8c8d': l_.remove()\n"
                                 # the 'proximal side, distal side' text overlaps the strip -> just below the strip (position only)
                                 "for t_ in A[1].texts:\n"
                                 "    if t_.get_text().strip() in ('근위 쪽', '원위 쪽', 'proximal side', 'distal side'):\n"
                                 "        t_.set_y(-0.40); t_.set_va('top')\n"
                                 # stretch the zone strip axes to the strip ends across the full cell width (the strip is a schematic of a mean-length array), all strip text black
                                 "_x0 = min(p_.get_x() for p_ in A[1].patches); _x1 = max(p_.get_x() + p_.get_width() for p_ in A[1].patches)\n"
                                 "A[1].set_xlim(_x0, _x1)\n"
                                 "for t_ in A[1].texts: t_.set_color('#000000')\n"
                                 # gray background: the drawing script makes 542 horizontal lines, one per haplotype from 0 to array length N (their width comes from lw)
                                 #   -> one filled area (step polygon, no edge) with the same rows (y +/- 0.5), same lengths and same color, so that it does not depend on the global line width
                                 "import matplotlib.patches as _mp\n"
                                 "_bg = sorted([(int(round(l_.get_ydata()[0])), max(l_.get_xdata()), l_) for l_ in A[0].get_lines() if str(l_.get_color()).lower() == '#dfe6e9'], key=lambda r_: r_[0])\n"
                                 "assert _bg and [r_[0] for r_ in _bg] == list(range(_bg[0][0], _bg[0][0] + len(_bg))), len(_bg)\n"
                                 "_pts = [(0, _bg[0][0] - 0.5)] + [p_ for y_, n_, _ in _bg for p_ in ((n_, y_ - 0.5), (n_, y_ + 0.5))] + [(0, _bg[-1][0] + 0.5)]\n"
                                 "_z = _bg[0][2].get_zorder()\n"
                                 "for _, _, l_ in _bg: l_.remove()\n"
                                 "A[0].add_patch(_mp.Polygon(_pts, closed=True, facecolor='#dfe6e9', edgecolor='none', lw=0, zorder=_z))\n"
                                 "print('   k 회색 바탕 : 가로선 %d 개 → 채운 면 1 개' % len(_bg))\n"}),
                       ("l", "FIG_G2c.py", "fig_G2c_split_gain", [0, 1], (30.0, 30.0), RMS, 112.4, ("h2", 6.0),
                        {"keep_title": True,
                         # remove the values in the annotations (269, 164, 133, 209) (values go to the value table and the Legend)
                         "retext_re": {k_: [[r",\s*\d+", ""], [r"\(2 단계에서 \d+ 로 채택\)", "(2 단계에서 채택)"], [r"\(taken at step 2 with \d+\)", "(taken at step 2)"]] for k_ in ("0", "1")},
                         # threshold text to the left end of the axes (just above the line, position only)
                         "post": "for a_ in (A[0], A[1]):\n"
                                 "    for t_ in a_.texts:\n"
                                 "        if t_.get_text().startswith(('벌점 문턱', 'penalty threshold')):\n"
                                 "            t_.set_x(a_.get_xlim()[0] + 0.6); t_.set_ha('left')\n"
                                 # the panel is narrow, so the end of the threshold text above the line ('73.8') would overlap the curve -> just below the line (the left side is empty because the curve is above the threshold)
                                 "for a_ in (A[0], A[1]):\n"
                                 "    _th = [l_ for l_ in a_.get_lines() if len(set(l_.get_ydata())) == 1 and str(l_.get_color()).lower() == '#2c3e50']\n"
                                 "    assert len(_th) == 1, len(_th)\n"
                                 "    for t_ in a_.texts:\n"
                                 "        if t_.get_text().startswith(('벌점 문턱', 'penalty threshold')):\n"
                                 # 1.8 pt below (annotate ignores an added transform, so use data coordinates, frame height 64 pt)
                                 "            _y0, _y1 = a_.get_ylim(); t_.set_y(_th[0].get_ydata()[0] - (_y1 - _y0) * 1.8 / 64.0); t_.set_va('top')\n"
                                 "            if EN: t_.set_text(__import__('re').sub(r'\\s+(\\S+)$', r'\\n\\1', t_.get_text()))\n"   # the English text is long, so its end value ('73.8') reaches where the curve drops below the threshold (depth 32 and beyond) -> value alone on the next line
                                 "    a_.set_xticks([10, 20, 30, 40])\n"   # in the narrow panel the automatic ticks leave only 20 and 40 -> pin the ticks
                                 # the cut-position text becomes a single 'depth N' (and its Korean equivalent) (same form for distal 10), position unchanged
                                 "import re as _re\n"
                                 "for a_ in (A[0], A[1]):\n"
                                 "    for t_ in a_.texts:\n"
                                 "        s_ = t_.get_text()\n"
                                 "        if s_.startswith(('벌점 문턱', 'penalty threshold')): continue\n"
                                 "        m_ = _re.search(r'(?:깊이|depth)\\s*(\\d+)', s_)\n"
                                 "        if m_: t_.set_text(('depth %s' if EN else '깊이 %s') % m_.group(1))\n"
                                 # the left end of 'depth 3' touches the depth 4 point (curve) -> move up and right (the curve drops after depth 5)
                                 "        if m_ and m_.group(1) == '3':\n"
                                 "            x_, y_ = t_.get_position(); t_.set_position((x_ + 1.2, y_ + 10))\n",
                         "labels": {"0": {"y": ["분할 이득 g", "Split gain g"],
                                          "x": ["자르는 깊이", "Cut position (depth)"]},
                                    "1": {"x": ["자르는 깊이", "Cut position (depth)"]}}})]),
  # m 97 so that the Korean ticks do not overlap, y left 39 for the 'cluster 1 (n 51)' labels, row bottom 27 for the Korean 'variant token' label
  (3.0, (27.0, 27.0), [("m", "FIG_O2b.py", "fig_O2b_boundary_class3", [0], (LM, LM), 4.0, 97.0, ("one",),
                        {"legend": {"0": {"loc": "upper right", "ncol": 1, "bbox": [1.0, 1.0], "borderaxespad": 0.0, "labelspacing": 0.1}}, "ylim_top": {"0": 10.5},
                         "drop_text": {"0": [r"re:\d+\.\d+"]},   # the numbers above the graph are removed (values in the value table)
                         "labels": {"0": {"y": ["경계당 변이 수", "Variants per boundary"]}},
                         # position slots are the 4 cohort zones (proximal 1-3, 4-17, middle, distal 1-10). The junction and interior rows apply when drawing from a table that has those slots
                         "xtick_sub": {"0": [["근위 접합 1-6", "근위 접합\n1-6"], ["원위 접합 4-10", "원위 접합\n4-10"], ["내부 (양쪽 ≥15)", "내부\n(양쪽 ≥15)"],
                                             ["proximal junction 1-6", "Prox.\n1-6"], ["distal junction 4-10", "Dist.\n4-10"],
                                             ["interior (both ≥15)", "Interior\n(both ≥15)"], ["all", "All"],
                                             ["근위 1-3", "근위\n1-3"], ["근위 4-17", "근위\n4-17"], ["원위 1-10", "원위\n1-10"],
                                             ["proximal 1-3", "Prox.\n1-3"], ["proximal 4-17", "Prox.\n4-17"], ["distal 1-10", "Dist.\n1-10"],
                                             ["middle", "Middle"]]},
                         # hatching (random positions) appears in all three classes, purple, red and yellow -> one common legend row 'random positions (hatched)'
                         "post": "import matplotlib.patches as _mp\nax_ = A[0]\n_h, _l = ax_.get_legend_handles_labels()\n"
                                 "_k = [i for i, t in enumerate(_l) if not (('무작위' in t) or ('random' in t))]\n"
                                 "_hz = _mp.Patch(fc='white', ec='#333333', hatch='////', lw=0.4)\n"
                                 "ax_.legend([_h[i] for i in _k] + [_hz], [_l[i] for i in _k] + (['random positions (hatched)'] if EN else ['무작위 자리 (빗금)']),\n"
                                 "           loc='upper right', bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=CS.FS_LEG, labelspacing=0.1, borderaxespad=0.0)\n"}),
  # y (fingerprint, 16 columns): right margin RMS 36 (colorbar right end = frame right line of the other rows, colorbar title in the outer margin)
                       ("y", "FIG_X12f.py", "fig_X12f_type_fingerprint", [0, 1], (39.0, 39.0), RMS, 191.4, ("cbar", 3.0, 3.0), {"table_axes": [0],   # white cell separator 0.30, number ink vertically centered in the cell
                        "post": "from matplotlib.transforms import offset_copy as _oc\nfrom matplotlib.textpath import TextPath as _TP\n"   # move the cell number ink to the vertical center of the cell (va=center centers the text box, so the number ink sat 0.22 pt too high)
                                "for _t in A[0].texts:\n"
                                "    _e = _TP((0, 0), _t.get_text(), size=_t.get_fontsize(), prop=_t.get_fontproperties()).get_extents()\n"
                                "    _t.set_va('baseline'); _t.set_transform(_oc(_t.get_transform(), fig=fig, y=-(_e.y0 + _e.y1) / 2, units='points'))\n"
                                "for _g in A[0].xaxis.get_minor_ticks() + A[0].yaxis.get_minor_ticks():\n"   # original ax.grid(which='minor', color='white', lw=1.0): get_xgridlines in clamp_marks only sees the major grid
                                "    _g.gridline.set_linewidth(0.3)\n"})]),
  (17.0, (18.0, 18.0), [("o", "P_make_paper_figures.py", "fig_P4_terminal_character", [0], (LM, LM), 3.0, 0.75, ("one",),   # row top 17: room for the two-row legend above o (ncol 1)
                        {"legend": {"0": {"loc": "lower left", "bbox": [0.0, 1.01], "ncol": 1, "borderaxespad": 0.0, "columnspacing": 1.0, "labelspacing": 0.1}},
                         # pin the ticks at 0.2 spacing (the automatic ticks come at 0.25)
                         "xticks": {"0": [[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"]]}}),   # the legend (two symbol descriptions) sits above the graph so that it does not overlap the plot
                       ("q", "FIG_T2.py", "fig_T2_arrangement", [3], (20.0, 20.0), 3.0, 0.95, ("one",),
                        {"legend": {"3": {"loc": "upper left", "borderaxespad": 0.0, "labelspacing": 0.1,
                                          "text": [["distal terminal-type (groups 169)", "Distal terminal (groups 169)"],
                                                   ["other distal minor types (groups 1477)", "Other distal minor (groups 1477)"],
                                                   ["proximal minor type (groups 1368)", "Proximal minor (groups 1368)"],
                                                   # under the terminal-type definition (10 or more of 12 diagnostics), group counts 167 and 1479
                                                   ["distal terminal-type (groups 167)", "Distal terminal (groups 167)"],
                                                   ["other distal minor types (groups 1479)", "Other distal minor (groups 1479)"],
                                                   # group counts on 541: 163, 1442, 1331
                                                   ["distal terminal-type (groups 163)", "Distal terminal (groups 163)"],
                                                   ["other distal minor types (groups 1442)", "Other distal minor (groups 1442)"],
                                                   ["proximal minor type (groups 1331)", "Proximal minor (groups 1331)"]]}},
                         "ylim_top": {"3": 4.0},
                         # the 'null = 1' text goes to the Legend. The numbers above the bars (1.36, 1.31 ...) overlap each other, so they are removed
                         "drop_text": {"3": ["귀무 =", "null =", r"re:\d+\.\d+"]}}),
                       ("r", "FIG_ZBDIV.py", "fig_term_neighbour_div", [0], (28.0, 28.0), RMS, 1.3, ("one",),
                        {"labels": {"0": {"y": ["이웃 유닛 서열 차이 (%)", "Neighbor diff. (%)"]}},
                         # the post below puts all value text at one row height
                         "ylim_top": {"0": 2.6}, "yticks": {"0": [0.0, 0.5, 1.0, 1.5, 2.0]},   # ticks only up to 2.0 so that the 2.5 grid line does not pass through the value text row (2.19-2.55)
                         "xtick_sub": {"0": [["Terminal\n↔\nneighbour\nordinary", "Term.\n↔ nbr.\nord."], ["Ordinary\n↔\nordinary\n(distal)", "Ord.\n↔ ord.\n(distal)"],
                                             ["Ordinary\n↔\nordinary\n(middle)", "Ord.\n↔ ord.\n(middle)"], ["Ordinary\n↔\nordinary\n(proximal)", "Ord.\n↔ ord.\n(prox.)"],
                                             ["Terminal\n↔\nterminal", "Term.\n↔\nterm."]]},
                         # the drawing script draws the quartiles (q25-q75) as lines of width 4.2, but
                         #   the global line width 0.45 (clamp_marks) also applies to them and the bars vanish -> replace them with filled bars (no edge) at the same position, color and alpha (values unchanged)
                         "post": "import matplotlib.patches as _mp\n"
                                 "for l_ in list(A[0].get_lines()):\n"
                                 "    if l_.get_linewidth() > 2.0:\n"
                                 "        _x = l_.get_xdata()[0]; _y = sorted(l_.get_ydata())\n"
                                 # quartile bar width 0.5 (at 0.18 the quartile 2.01-2.10 of the first branch looks smaller than the point)
                                 "        A[0].add_patch(_mp.Rectangle((_x - 0.25, _y[0]), 0.5, _y[1] - _y[0], facecolor=l_.get_color(), alpha=l_.get_alpha(), lw=0, zorder=2))\n"
                                 "        l_.remove()\n"
                                 #   value text (median % and n): all at the same height just above the branches (above the highest p90)
                                 "_T = [t_ for t_ in A[0].texts if '%' in t_.get_text()]\n"
                                 "_yl = max(t_.get_position()[1] for t_ in _T)\n"
                                 "for t_ in _T: t_.set_position((round(t_.get_position()[0]), _yl)); t_.set_ha('center'); t_.set_va('bottom')\n"})]),
  (3.0, (21.0, 21.0), [("v", "FIG_T1.py", "fig_T1_terminal_lineage", [2], (LM, LM), 3.0, 1.15, ("one",),
                        {"ylim_top": {"2": 1.2}, "legend": {"2": {"loc": "upper left", "ncol": 2, "bbox": [0.0, 1.0], "borderaxespad": 0.0, "columnspacing": 1.0,
                                                                "text": [["within the same person", "Same person"], ["between different people", "Different people"]]}},
                         "xtick_sub": {"2": [["terminal-type\ndiagnostic included", "Term.\n+diag."], ["terminal-type\ndiagnostic excluded", "Term.\n−diag."],
                                             ["ordinary\ndiagnostic included", "Ord.\n+diag."], ["ordinary\ndiagnostic excluded", "Ord.\n−diag."]]},
                         "labels": {"2": {"y": ["자카드 닮음 (평균)", "Jaccard similarity (mean)"]}},
                         # remove the numbers above the bars
                         "drop_text": {"2": [r"re:\d\.\d\d"]}}),
                       ("w", "FIG_F9.py", "fig_F9_population", [3], (28.0, 28.0), 3.0, 1.0, ("one",),
                        {"ylim_top": {"3": 0.2}, "legend": {"3": {"loc": "upper left", "borderaxespad": 0.0, "labelspacing": 0.1,
                                          # legend text: population name only (the counts go to the Legend text)
                                          "text_re": [[r"^(\S+)\s+(?:전체|all)\s+(\d+)\s*[,·]\s*(?:말단형 보유|terminal-type)\s+(\d+)$", r"\1"]]}},
                         "labels": {"3": {"y": ["집단 내 하플로타입 비율", "Frac. of haplotypes"],
                                          "x": ["말단형 깊이 집합 크기", "Size of terminal-depth set (positions)"]}}}),
                       ("x", "P_make_paper_figures.py", "fig_P3_stability", [1], (40.0, 40.0), RMS, 0.85, ("one",),
                        {"xticks": {"1": [[1, 5, 10, 15, 20], ["1", "5", "10", "15", "20"], ["1", "5", "10", "15", "20"]]},
                         "labels": {"1": {"y": ["하플로타입", "Haplotype"]}},   # ticks 'cell line A/B', axis title haplotype
                         # remove the HepG2 rows (the remaining rows are pulled together with no gaps)
                         "post": "ax_ = A[1]\nlabs_ = [t.get_text() for t in ax_.get_yticklabels()]\n"
                                 "drop_ = [i for i, s in enumerate(labs_) if s.startswith('HepG2')]\nkeep_ = [i for i in range(len(labs_)) if i not in drop_]\n"
                                 "ny_ = {o: n for n, o in enumerate(keep_)}\n"
                                 "for c_ in list(ax_.collections):\n"
                                 "    off_ = c_.get_offsets()\n"
                                 "    if len(off_) == 0: continue\n"
                                 "    yy_ = int(round(float(off_[0][1])))\n"
                                 "    if yy_ in drop_: c_.remove()\n"
                                 "    else: off_[:, 1] = ny_[yy_]; c_.set_offsets(off_)\n"
                                 "for l_ in list(ax_.get_lines()):\n"
                                 "    yd_ = l_.get_ydata(); yy_ = int(round(float(yd_[0])))\n"
                                 "    if yy_ in drop_: l_.remove()\n"
                                 "    else: l_.set_ydata([ny_[yy_]] * len(yd_))\n"
                                 # haplotype names: slot 1 = A, slot 2 = B (hap1 -> A, hap2 -> B)
                                 # ticks combine cell line and A/B, like 'K562 A'
                                 "_HM = {'GM12878_hap1': 'GM12878 A', 'GM12878_hap2': 'GM12878 B', 'GM12878_hapA': 'GM12878 A', 'GM12878_hapB': 'GM12878 B',\n"
                                 "       'K562_hap1': 'K562 A', 'K562_hap2': 'K562 B', 'K562_hapA': 'K562 A', 'K562_hapB': 'K562 B'}\n"
                                 "assert all(labs_[i] in _HM for i in keep_), [labs_[i] for i in keep_]\n"
                                 "ax_.set_yticks(range(len(keep_))); ax_.set_yticklabels([_HM.get(labs_[i], labs_[i]) for i in keep_]); ax_.set_ylim(-0.6, len(keep_) - 0.4)\n"
                                 "print('   x 줄 뺌', [labs_[i] for i in drop_], '남김', [labs_[i] for i in keep_])\n"})]),
]


# The ROWS redefined below is the arrangement used.
#   k (segmentation table): middle cell of the fourth main row. Size and margins from zz in cellstyle.ROWS (frame top and bottom match the main row), content settings unchanged
#   second row = l (split gain), g (rank table). l width 165 (one panel about 77 pt), g 158.4, top 31.5 (g slanted header), bottom 30
#     g left margin 4 (table row names are drawn inside the axes)
#   _JOPT (j: y-axis top 2.1 into the empty band above the data (max 1.60), ticks 0.6 to 1.6, the D1096:1 legend row split into two lines before the parenthesis)
#     and _YOPT (y: cell numbers removed, YlGnBu gradient, vertical column names) are not used by ROWS
_C ={c[0]: c for _t, _b, cells in ROWS for c in cells}
def _mod(name, **kw):
    c = list(_C[name]); ix = {"axes": 3, "l": 4, "r": 5, "wt": 6, "sub": 7, "opt": 8}
    for k_, v_ in kw.items(): c[ix[k_]] = v_
    return tuple(c)
_JOPT = dict(_C["j"][8]); _JOPT["ylim_top"] = {"5": 2.1}; _JOPT["yticks"] = {"5": [0.6, 0.8, 1.0, 1.2, 1.4, 1.6]}
_JOPT["post"] = _C["j"][8]["post"] + r"""for t_ in _lg.get_texts():
    s_ = t_.get_text()
    if 'D1096' in s_: t_.set_text(__import__('re').sub(r'\s*\(', '\n(', s_, count=1))
"""
_YOPT = {"table_axes": [0], "drop_text": {"0": [r"re:\d\.\d\d"]}, "xtick_rot": {"0": 90},
         "post": "import numpy as _np\nfrom matplotlib.colors import ListedColormap as _LC\nimport matplotlib as _mpl\n"
                 "_im = A[0].images[0]\n_a = _np.ma.masked_less_equal(_np.asarray(_im.get_array(), dtype=float), 0.0)\n"
                 "_cm = _LC(_mpl.colormaps['YlGnBu'](_np.linspace(0.25, 1.0, 256))); _cm.set_bad('#ffffff')\n"
                 "_im.set_data(_a); _im.set_cmap(_cm); _im.set_clim(0.0, 1.0)\n"
                 "for _g in A[0].xaxis.get_minor_ticks() + A[0].yaxis.get_minor_ticks():\n"
                 "    _g.gridline.set_linewidth(0.3)\n"}
# q legend (gap fold, width 64): placed in the empty space above the graph (within the row top margin 19.5) so that it does not cross the right frame line,
#   and '(groups N)' and its Korean equivalent are removed from its text (the counts go to the Legend text)
#   v (Jaccard, width 89): the Korean ticks 'diagnostic included' and 'diagnostic excluded' touch -> Korean '+diag.' and '-diag.' forms (same form as the English '+diag.' and '-diag.', meaning in the Legend text)
#   r: alpha of the quartile bars (the filled bars made in post) is 0.95 (color unchanged)
_ZOPT = {"q": dict(_C["q"][8], legend={"3": dict(_C["q"][8]["legend"]["3"], loc="lower left", bbox=[0.0, 1.01],
                                                text_re=[[r"\s*\((?:groups|그룹)\s*\d+\)", ""]])}),
         "r": dict(_C["r"][8], post=_C["r"][8]["post"] + r"""import matplotlib.patches as _mp2
for p_ in A[0].patches:
    if isinstance(p_, _mp2.Rectangle) and p_.get_zorder() == 2: p_.set_alpha(0.95)
"""),
         "v": dict(_C["v"][8], post=_C["v"][8].get("post", "") + r"""if not EN:
    A[2].set_xticks(A[2].get_xticks()); A[2].set_xticklabels([t_.get_text().replace('진단 포함', '+진단').replace('진단 제외', '-진단') for t_ in A[2].get_xticklabels()])
""")}
_ZZ = ("MAIN", "zz", _C["k"])
ROWS = [ROWS[0],
        (31.5, (30.0, 30.0), [_mod("l", l=(LM, LM), r=4.0, wt=165.0), _mod("g", l=(4.0, 4.0), r=RMS, wt=158.4)]),
        _ZZ,
        # m y row: the cells of y carry their numbers, so y has one axes and no colorbar (the colorbar axes A[1] is hidden when dropped from the list, _retarget).
        #   The fingerprint frame reaches the end of the cell width = frame right line of the other rows (right margin RMS unchanged)
        (lambda r_: (r_[0], r_[1], [c_ if c_[0] != "y" else _mod("y", axes=[0], sub=("one",)) for c_ in r_[2]]))(
            [r_ for r_ in ROWS if [c[0] for c in r_[2]] == ["m", "y"]][0])] + \
       [("MAIN", "z" + n_, _mod(n_, opt=_ZOPT.get(n_, _C[n_][8]))) for n_ in ("o", "q", "r", "v", "w", "x")]   # supplementary o q r v w x: sizes from the z* slots of cellstyle.ROWS
assert [c[0] for r_ in ROWS if r_[0] != "MAIN" for c in r_[2]] == ["a", "b", "e", "l", "g", "m", "y"], "re-layout of rows"


def build(en):
    k = 1 if en else 0
    spec, lay = {}, {}
    for t, b, cells in ROWS:
        if t == "MAIN":                                  # cell inside a main row. Size and margins from that slot in cellstyle.ROWS (b = slot name)
            import cellstyle as _CS
            _d = _CS.layout(en)[b]; cells = [cells]
        else:
            assert cells[0][4][k] == LM and cells[-1][5] == RMS, cells[0][0]
            h = t + FH + b[k]
            fixed = sum(c[4][k] + c[5] for c in cells)
            avail = W_USE - fixed
            ws = sum(c[6] for c in cells)
        for name, code, save, axes, l, r, wt, sub, opt in cells:
            if t == "MAIN":
                assert abs(_d["h"] - _d["t"] - _d["b"] - FH) < 1e-9, _d
                axw, w, h, x0, y0 = _d["axw"], _d["w"], _d["h"], _d["l"], _d["b"]
                l = (_d["l"], _d["l"]); r = _d["r"]; bk = (_d["b"], _d["b"])
            else:
                axw = avail * wt / ws
                w = l[k] + axw + r
                x0, y0 = l[k], b[k]; bk = b
            if sub[0] == "one":
                rects = [[x0, y0, axw, FH]]
            elif sub[0] == "h2":
                g = sub[1]; aw = (axw - g) / 2.0
                rects = [[x0, y0, aw, FH], [x0 + aw + g, y0, aw, FH]]
            elif sub[0] == "v":
                fr, g = sub[1], sub[2]; ht = (FH - g) * fr
                rects = [[x0, y0 + FH - ht, axw, ht], [x0, y0, axw, FH - g - ht]]
            elif sub[0] == "below":       # the main graph axes match the row graph frame (64 pt), the attached strip goes in the bottom margin under the frame
                g, h2 = sub[1], sub[2]
                assert bk[k] >= g + h2, (name, bk[k], g, h2)
                rects = [[x0, y0, axw, FH], [x0, y0 - g - h2, axw, h2]]
            elif sub[0] == "cbar":
                cw, g = sub[1], sub[2]
                rects = [[x0, y0, axw - g - cw, FH], [x0 + axw - cw, y0, cw, FH]]
            cell = {"name": "cellS2_%s" % name, "w": w, "h": h, "axes": axes, "rects": rects}
            cell.update({kk: vv for kk, vv in opt.items() if kk != "fig_legend"})
            if opt.get("fig_legend"):
                cell["fig_legend"] = {"x": x0 + axw / 2.0, "y": 1.0, "loc": "lower center"}
                if isinstance(opt["fig_legend"], dict): cell["fig_legend"].update(opt["fig_legend"])   # ncol and similar
            spec.setdefault(save, []).append(cell)
            lay[name] = {"code": code, "save": save, "w": w, "h": h, "frame": [x0, y0, axw, FH]}
    return spec, lay


for L, en in (("ko", False), ("en", True)):
    spec, lay = build(en)
    json.dump(spec, open(os.path.join(HERE, "spec_supp_%s.json" % L), "w", encoding="utf8", newline="\n"), ensure_ascii=False, indent=1)
    json.dump(lay, open(os.path.join(HERE, "supp_layout_%s.json" % L), "w", encoding="utf8", newline="\n"), ensure_ascii=False, indent=1)
    print(L, " ".join("%s:%.0fx%.0f" % (n, d["w"], d["h"]) for n, d in lay.items()))
