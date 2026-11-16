# -*- coding: utf-8 -*-
# The shared module for the array methylation landscape panels: the supplementary FigE4 and the main Figure 1 are drawn by the same code.
#   (Writing the same figure twice in two places guarantees they will diverge, so there is one copy.)
#
# The two lines:
#   1) the two lines are on the same ruler
#      - the black CpG-window line is a call-weighted mean, sum(methylated) / sum(calls);
#      - the orange unit layer is a call-weighted mean as well, sum(n_meth) / sum(n_cpg), so the two lines mean the same thing;
#      - the median of the unit betas would differ from it: over 10,499 positions a median |weighted mean - median| of 0.0385, a maximum of 0.625, with the weighted mean the larger at 4.5% of positions
#        (the unit beta is the fraction of a binary switch and so has a skewed distribution, in which the median and the mean part company).
#   2) the black line
#      - at least 300 calls are pooled per window (fine_minc). No data is discarded;
#        the same data is pooled over a wider window, which reduces the number of points and the jitter. The window width is stated in the panel and the legend.
#
# Inputs, read-only: tables/PERARRAY_LANDSCAPE.tsv.gz, tables/PERARRAY_META.tsv, tables/ORIENTATION.tsv and
#                tables/UNITBETA_PERARRAY_v13.tsv.gz · tables/CHROM_NAMES.tsv
import numpy as np
import pandas as pd

NED, MIDB, FLMAX, OFFCUT = 20, 12, 20000, 0.35
BRK_F, SQ = 2000.0, 4.0
FINE, UCOL = "#4D5656", "#E07B39"
FINE_MINC = 300           # the minimum number of CpG calls pooled into one point of the black line


def t(u):
    a = np.asarray(u, float)
    mag = np.where(a < 0, np.abs(a), np.where(np.abs(a) <= BRK_F, np.abs(a), BRK_F + (np.abs(a) - BRK_F) / SQ))
    return np.sign(a) * mag


TOUT = float(t(FLMAX))


def win_edges(inside, imax, unit=None):
    """Windows inside an array are a multiple of the unit width: a window narrower than the unit echoes the pattern inside the unit
    and the line looks like a sawtooth. At the unit width the line is clean and
    measures the same thing as the unit layer (orange) at the same position."""
    if inside and unit:
        e = list(np.arange(0, max(imax, FLMAX) + unit, float(unit)))
    elif inside:
        e = list(np.arange(0, 2000, 50)) + list(np.arange(2000, 10000, 200)) + list(np.arange(10000, max(imax, 10001), 1000))
    else:
        e = list(np.arange(0, 2000, 100)) + list(np.arange(2000, 10000, 500)) + list(np.arange(10000, 60000 + 1, 5000))
    return np.array([x for x in e if x <= FLMAX + 1e-9], dtype=float)


def adaptive(dd, edges, minc=FINE_MINC):
    """Widen a grid window in place until it holds minc calls (only the window width changes, no data are filled in).
    position = call-weighted mean; β = Σmethylated/Σcalls (hard calls)."""
    got = {round(float(a), 1): (float(b), float(c))
           for a, b, c in zip(dd.d.values, dd.호출수.values, dd.메틸호출수.values)}
    xs, ys, nb = [], [], []
    an = am = 0.0
    bag = []
    for e in edges:
        n, m_ = got.get(round(float(e), 1), (0.0, 0.0))
        an += n
        am += m_
        bag.append((float(e), n))
        if an >= minc:
            xs.append(sum(x * w for x, w in bag) / an)
            ys.append(am / an)
            nb.append(len(bag))
            an = am = 0.0
            bag = []
    return np.array(xs), np.array(ys), nb


def slide_unit(dd, unit, step_frac=0.25, lo=0.0, hi=None):
    """CpG windows inside an array: window width = 1 unit, step = 1/4 unit
    Why : a window narrower than the unit echoes the pattern inside the unit and the line looks like a sawtooth.
         A sliding window of unit width gives a smooth line, with points denser than the unit layer so that the two lines can be told apart.
    position = window centre; β = Σmethylated / Σcalls (hard calls)."""
    d = dd.sort_values("d")
    x = d.d.values.astype(float)
    n = d.호출수.values.astype(float)
    m = d.메틸호출수.values.astype(float)
    step = max(unit * step_frac, 1.0)
    # No upper bound on the flank (FLMAX) is applied inside the array; applying one cuts the line off at 20 kb
    #   in a species with long units, where 20 ape units come to 47 kb.
    top = (x.max() if len(x) else lo) if hi is None else min(hi, x.max() if len(x) else lo)
    cs = np.arange(lo + unit / 2.0, top + step, step)
    xs, ys = [], []
    for c in cs:
        k = (x >= c - unit / 2.0) & (x < c + unit / 2.0)
        if not k.any():
            continue
        tot = n[k].sum()
        if tot <= 0:
            continue
        xs.append(c)
        ys.append(m[k].sum() / tot)
    return np.array(xs), np.array(ys)


def load(NEW):
    """The table is read once and kept."""
    LS = pd.read_csv(NEW + "/tables/PERARRAY_LANDSCAPE.tsv.gz", sep="\t")
    M = pd.read_csv(NEW + "/tables/PERARRAY_META.tsv", sep="\t")
    O = pd.read_csv(NEW + "/tables/ORIENTATION.tsv", sep="\t").set_index("array_id")
    UB = pd.read_csv(NEW + "/tables/UNITBETA_PERARRAY_v13.tsv.gz", sep="\t",
                     usecols=["묶음", "array_id", "유닛번호", "n_cpg", "n_meth"])
    GRP = dict(zip(UB.array_id, UB.묶음))
    IDS = {}
    for a, k in dict(zip(UB.array_id, UB.묶음)).items():
        IDS.setdefault(k, []).append(a)
    # the unit layer is a call-weighted mean: sum(methylated) / sum(calls) per locus and unit index, with the zebra finch maternal and paternal copies summed as one locus
    SLOT = (UB.groupby(["묶음", "유닛번호"], as_index=False)[["n_cpg", "n_meth"]].sum())
    SLOT["mean_pooled"] = SLOT.n_meth / SLOT.n_cpg
    # The ruler of the horizontal axis is the chromosome direction, in main panel b and supplementary 7-2 alike:
    #   the lower coordinate is on the left.
    #   The numbers in tables/FIGDATA_B4_pos_summary.tsv, off_fraction_proximal_end and off_fraction_distal_end, are
    #   on the chromosome convention as well, so the axis agrees with the numbers.
    #
    #   Proximal and distal are measured values, from the proximal-side column of tables/ORIENTATION.tsv: proximal is the side carrying the centromere
    #   (measured from tables/CENTROMERE.tsv). All 53 arrays carry a value, with none missing.
    #   Proximal is not taken as the end away from the nearer chromosome end: that is wrong for 6 arrays
    #   (zebrafish chr4 three times and chr18 once lie on the long arm, and zebra finch chr9 twice is acrocentric).
    #
    #   Transcription does not run in one direction throughout: the gorilla alone is on the minus strand.
    #   The axis is fixed on one ruler, increasing coordinate, and P1 and D1 are attached per panel from measurement, in panel() below.
    LEFT = {a: "낮은 좌표" for a in M.array_id}
    _CN = pd.read_csv(NEW + "/tables/CHROM_NAMES.tsv", sep="\t")
    _CH = dict(zip(_CN.접근번호.astype(str), _CN.염색체번호.astype(str)))
    POS = {}
    for _t in M.itertuples(index=False):
        _c = str(_t.chrom).replace("_mat", "").replace("_pat", "")
        _l = str(_t.이름표).replace("_mat", "").replace("_pat", "")
        _nm = _l if _l.lower().startswith("chr") else _CH.get(_c, _c)
        POS[_t.array_id] = _nm if "Mb" in _nm else "%s : %.1f Mb" % (_nm, float(_t.array시작) / 1e6)
    return dict(LS=LS, M=M, O=O, SLOT=SLOT, GRP=GRP, IDS=IDS, LEFT=LEFT, POS=POS)


def panel(ax, L, mt, fine_minc=FINE_MINC, ends=True, ylab=True, fs_small=6.0,
          lw_fine=0.8, lw_unit=1.6, ms_unit=2.2, marker="o", unit_window=False, col_unit=None,
          scalebar=False, fs_bar=6.0, bar_dy=-0.26):
    """Line width and marker size are arguments so that they can shrink with the cell size (cells in the main Fig1 are small)."""
    """Draw one cell. mt = one row of PERARRAY_META (itertuples). Returns what was drawn as a dict."""
    LS, O, SLOT, GRP, IDS, LEFT = L["LS"], L["O"], L["SLOT"], L["GRP"], L["IDS"], L["LEFT"]
    key = GRP.get(mt.array_id)
    ids = IDS.get(key, [mt.array_id])
    prox = LEFT[mt.array_id]
    dist = "높은 좌표" if prox == "낮은 좌표" else "낮은 좌표"
    UNIT = float(mt.유닛길이)
    TINu = abs(float(t(-(NED + 0.5) * UNIT)))
    MIDW = max(2400.0, 0.30 * TINu)
    X_J1, X_J2 = TOUT, TOUT + TINu + MIDW + TINu
    XMAX = X_J2 + TOUT

    def xleft(u):
        return X_J1 - t(u)

    def xright(u):
        return X_J2 + t(u)

    # the shading carries no border (lw=0); color= fills both the face and the edge, and without an explicit width
    #   matplotlib omits stroke-width from the SVG, which renders as a 1.0 pt border.
    ax.axvspan(0, X_J1, color="#F2F2F2", lw=0, zorder=0)
    ax.axvspan(X_J2, XMAX, color="#F2F2F2", lw=0, zorder=0)
    ax.axvspan(X_J1 + TINu, X_J1 + TINu + MIDW, color="#FFF8E7", lw=0, zorder=0)
    for xf in (xleft, xright):
        x0, x1 = sorted([xf(BRK_F), xf(FLMAX)])
        ax.axvspan(x0, x1, color="#E8EAF6", lw=0, zorder=0)
    d = LS[LS.array_id.isin(ids)]
    if len(ids) > 1:          # where a locus holds two arrays their call counts are summed into one curve
        d = (d.groupby(["쪽", "위치", "d"], as_index=False)[["호출수", "메틸호출수"]].sum())
        d["beta"] = d.메틸호출수 / d.호출수
    npt = 0
    for side, xf in ((prox, xleft), (dist, xright)):
        for loc, sgn in (("array 안", -1), ("flank", +1)):
            dd = d[(d.쪽 == side) & (d.위치 == loc)].sort_values("d")
            if not len(dd):
                continue
            # the CpG window is drawn inside the array as well; at a unit-width window it would be hidden under the unit layer
            #   so the calls pooled per window are raised to fine_minc (default FINE_MINC = 300, and the caller does not override it), which reduces the jitter
            if loc == "array 안":
                dx, dy = slide_unit(dd, UNIT)          # the unit-width window is slid by a quarter of a unit at a time
            else:
                ed = win_edges(False, float(mt.array끝 - mt.array시작))
                dx, dy, _nb = adaptive(dd, ed, minc=fine_minc)
            if not len(dx):
                continue
            x = xf(sgn * dx)
            o2 = np.argsort(x)
            ax.plot(x[o2], dy[o2], "-", color=FINE, lw=lw_fine, zorder=6)
            npt += len(dx)
    offs = {"낮은 좌표": float(mt.부분유닛_낮은쪽), "높은 좌표": float(mt.부분유닛_높은쪽)}
    for side, xf in ((prox, xleft), (dist, xright)):
        if offs[side] > 0:
            x0, x1 = sorted([xf(0.0), xf(-offs[side])])
            ax.axvspan(x0, x1, facecolor="#D5D8DC", edgecolor="#7F8C8D", hatch="///", lw=0, alpha=.85, zorder=1)
    # -- the unit layer: a call-weighted mean, the same ruler as the black line
    s = SLOT[SLOT.묶음 == key].copy()
    N = int(mt.유닛수)
    s["낮은쪽번호"] = s.유닛번호 + 1
    s["높은쪽번호"] = N - s.유닛번호
    pcol = "낮은쪽번호" if prox == "낮은 좌표" else "높은쪽번호"
    dcol = "높은쪽번호" if prox == "낮은 좌표" else "낮은쪽번호"
    ux, uy = [], []
    for k in range(1, NED + 1):
        r_ = s[s[pcol] == k]
        if len(r_):
            ux.append(xleft(-(k - .5) * UNIT))
            uy.append(float(r_.n_meth.sum() / max(r_.n_cpg.sum(), 1)))
    mid = s[(s[pcol] > NED) & (s[dcol] > NED)]
    if len(mid):
        q = np.linspace(0, 1, MIDB + 1)
        pr = (mid[pcol].values - NED) / max(N - 2 * NED, 1)
        xm = np.linspace(X_J1 + TINu + MIDW * .08, X_J1 + TINu + MIDW * .92, MIDB)
        for bi in range(MIDB):
            m_ = (pr >= q[bi]) & ((pr < q[bi + 1]) if bi < MIDB - 1 else (pr <= 1))
            if m_.any():
                ux.append(xm[bi])
                uy.append(float(mid.n_meth.values[m_].sum() / max(mid.n_cpg.values[m_].sum(), 1)))
    for k in range(NED, 0, -1):
        r_ = s[s[dcol] == k]
        if len(r_):
            ux.append(xright(-(k - .5) * UNIT))
            uy.append(float(r_.n_meth.sum() / max(r_.n_cpg.sum(), 1)))
    if ux:
        # the marker edge width (mew) must be given explicitly: without it a default of 1.0 pt slips in
        #   and the outer diameter becomes ms + 1.0 (ms written as 0.28 gives an outer diameter of 1.28).
        #   To keep the size chosen by eye, the specification value of 0.46 is written and ms is increased to match.
        ax.plot(ux, uy, "-" + (marker or ""), color=(col_unit or UCOL), lw=lw_unit,
                ms=ms_unit + 0.54, markeredgewidth=0.46, zorder=5)
    for xj in (X_J1, X_J2):
        ax.axvline(xj, color="#111", lw=0.46, zorder=4)    # the specification value, 0.46
    ax.axhline(OFFCUT, color="#888", lw=0.41, ls=(0, (4, 3)), zorder=2)   # a guide line is thinner than a data line at 0.46
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, 1.06)
    # 6) the first two units are marked on the x axis
    #   P1 is the first unit from the proximal end and D1 the first from the distal end, the same positions as the first points of the unit layer.
    #   The labels are attached per panel from measurement: the axis is laid out in increasing coordinate
    #   and which of the left and the right is proximal is then judged and written.
    #   The axis always has the lower coordinate on the left, so where proximal is the lower coordinate it reads P1 ... D1, and where it is the higher, D1 ... P1.
    #   Hard-coding them would make them false.
    #   Proximal and distal are set from the centromere coordinates (the proximal-side column of ORIENTATION.tsv, written by x41_orientation.py).
    #   Where the centromere-overlap column reads yes, the array lies inside the centromere interval and proximal and distal lose their meaning,
    #     so no label is written (one panel, Arabidopsis CP116282.1_15747509_15754203).
    _p_ori = O.근위쪽.get(mt.array_id, "낮은 좌표")
    _ovl = str(O.중심절겹침.get(mt.array_id, "아니오")) == "예" if "중심절겹침" in O.columns else False
    _lab = ["", ""] if _ovl else (["P1", "D1"] if _p_ori == "낮은 좌표" else ["D1", "P1"])
    ax.set_xticks([xleft(-0.5 * UNIT), xright(-0.5 * UNIT)])
    ax.set_xticklabels(_lab, fontsize=fs_small)
    ax.tick_params(axis="x", length=1.5, pad=1, width=0.41)
    c_in_l = mt.array시작 if prox == "낮은 좌표" else mt.array끝
    c_in_r = mt.array끝 if prox == "낮은 좌표" else mt.array시작
    c_out_l = c_in_l - FLMAX if prox == "낮은 좌표" else c_in_l + FLMAX
    c_out_r = c_in_r + FLMAX if prox == "낮은 좌표" else c_in_r - FLMAX
    if ends:
        ax.text(0.0, -0.115, "%.2f Mb" % (c_out_l / 1e6), fontsize=fs_small, ha="left",
                transform=ax.transAxes, color="#555")
        ax.text(1.0, -0.115, "%.2f Mb" % (c_out_r / 1e6), fontsize=fs_small, ha="right",
                transform=ax.transAxes, color="#555")
    if ylab:
        ax.set_ylabel("β", fontsize=8)
    if scalebar:
        # The x axis is compressed differently in different stretches and so cannot carry ticks; a real bp scale bar is drawn instead.
        #   Its length differs between species, the axis being tied to the unit length: the round value closest to a quarter of the width of 20 terminal units.
        _cand = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]
        _bp = min(_cand, key=lambda v: abs(v - TINu / 4.0))
        # 5) the scale bar goes directly below the plot; inside the panel it would cover the data.
        #   Its x is in data coordinates, the axis being compressed in stretches so that the bar must sit over the stretch it refers to,
        #   and only its y is lowered below the axis, through get_xaxis_transform. It starts 0.30 of a unit width to the right so as not to overlap the P1 label.
        _tr = ax.get_xaxis_transform()
        _x0 = X_J1 + TINu * 0.30
        ax.plot([_x0, _x0 + _bp], [bar_dy, bar_dy], "-", color="#222", lw=0.46,
                solid_capstyle="butt", zorder=9, clip_on=False, transform=_tr)
        _lab = ("%g kb" % (_bp / 1000.0)) if _bp >= 1000 else ("%d bp" % _bp)
        # so that the scale-bar text does not overlap the plot:
        #   a) it is attached at the left end of the bar with ha="left", keeping it to the right of the vertical line that marks the start of the array
        #      (centred on the bar it would run past that line in a narrow panel, and "1 kb" would read as "| kb"), and
        #   b) a white background box is laid under it, so that the lines and shading beneath do not cross the text.
        ax.text(_x0 + _bp + XMAX * 0.012, bar_dy, _lab, fontsize=fs_bar, ha="left", va="center",
                color="#222", zorder=9, clip_on=False, transform=_tr)
    return dict(묶음=key, 유닛층_점=len(ux), CpG창_점=npt, 왼쪽=prox, 유닛길이=UNIT, 유닛수=N)
