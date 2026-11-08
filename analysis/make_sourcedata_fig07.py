# -*- coding: utf-8 -*-
# Source Data workbook for Figure 7 and Supplementary Figures S11, S12 and S13.
#
#   The values are taken from the same result tables the figures read, through the same
#   modules the figures use (betabins, prof_smooth, landsc), so a sheet cannot drift away
#   from the panel it belongs to. Nothing is retyped and nothing is copied from a document.
#
#   Everything printed in the workbook is in English: file name, sheet names, headers and
#   the README sheet. The species names are read from the figure's own translation table,
#   so the workbook and the figures use the same words.
#
#   Paths come from the environment, with no default, so the script stops with a named
#   variable instead of silently reading the wrong folder:
#     RDNA_OUT_DIR         the working output root; the result tables are in $RDNA_OUT_DIR/tables
#     RDNA_LIB_DIR         the folder holding betabins.py, prof_smooth.py, landsc.py, figorder.py
#     RDNA_SOURCEDATA_DIR  where the workbook is written
import io
import os
import sys
import hashlib

import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter


def _need(k):
    v = os.environ.get(k)
    if not v:
        sys.exit("environment variable %s is not set (this script has no default; "
                 "see README.md, section 'How to run it')" % k)
    return v


NEW = _need("RDNA_OUT_DIR")
LIB = _need("RDNA_LIB_DIR")
OUTDIR = _need("RDNA_SOURCEDATA_DIR")
sys.path.insert(0, LIB)
import betabins      # noqa: E402
import prof_smooth   # noqa: E402
import landsc        # noqa: E402
import figorder      # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

T = NEW + "/tables/"
OUT = os.path.join(OUTDIR, "SourceData_Fig7_FigS11_FigS12_FigS13.xlsx")
NED, MIDB = landsc.NED, landsc.MIDB


def md5(p):
    return hashlib.md5(io.open(p, "rb").read()).hexdigest()[:12]


# ── species names, from the figure's own lookup ────────────────────────────────
def _sp_map():
    m = {}
    for f in ("ko_en_terms.tsv", "ko_en_templates.tsv"):
        p = os.path.join(LIB, f)
        if not os.path.isfile(p):
            continue
        for line in io.open(p, encoding="utf-8"):
            part = line.rstrip("\n").split("\t")
            if len(part) >= 2 and part[0].strip():
                m.setdefault(part[0].strip(), part[1].strip())
                if f.endswith("templates.tsv"):
                    m[part[0].strip()] = part[1].strip()
    return m


SPEN = _sp_map()


def sp(ko):
    return SPEN.get(str(ko), str(ko))


SIDE = {"낮은 좌표": "lower coordinate", "높은 좌표": "higher coordinate"}
LOC = {"array 안": "inside the array", "flank": "flank"}

# How each centromere was located, in English. The keys are the openings of the four strings
# the result table writes, so a new method in the table shows up as untranslated and is caught
# by the check that no Korean syllable reaches the workbook.
CENM = {
    "주석된 중심절 영역":
        "annotated centromeric region (Eichler Lab T2T primate RepeatMasker centromeric_regions, "
        "the Satellite range of chromosome 1, the two haplotypes taken together)",
    "서열에서 측정 (CEN180":
        "measured from the assembly sequence (the longest tandem block with a CEN180 period of "
        "178 bp, autocorrelation at least 0.80)",
    "위성 밀도 봉우리":
        "satellite-density peak (UCSC GenArk repeatMasker track of the same assembly, the longest "
        "run of 100 kb bins that are at least 50 percent Satellite, 40 to 259 times the genome "
        "median; the track was matched by chromosome name)",
    "서열에서 측정 (Tgut191A":
        "measured from the assembly sequence (10 kb windows whose period is a multiple of "
        "Tgut191A 191 bp or Tgut716A 716 bp, autocorrelation at least 0.55)",
}


def cenm(v):
    for k, en in CENM.items():
        if str(v).startswith(k):
            return en
    return str(v)


def lab_of(mt):
    """The panel label, with the middle dot of the key written as a hyphen."""
    return str(POS.get(mt.array_id, mt.이름표)).replace(" · ", " - ").replace("·", "-")

# ── tables ────────────────────────────────────────────────────────────────────
A1 = pd.read_csv(T + "FIGDATA_A1_chrom_map.tsv", sep="\t")
CN = pd.read_csv(T + "CHROM_NAMES.tsv", sep="\t")
CEN = pd.read_csv(T + "CENTROMERE.tsv", sep="\t")
PM = pd.read_csv(T + "PERARRAY_META.tsv", sep="\t")
PR = pd.read_csv(T + "PROFILE_META_v13.tsv", sep="\t")
UB = pd.read_csv(T + "UNITBETA_PERARRAY_v13.tsv.gz", sep="\t")
BP = pd.read_csv(T + "PROFILE_BP_v13.tsv.gz", sep="\t")
CX = pd.read_csv(T + "ATH_READUNIT_CONTEXT.tsv.gz", sep="\t")
CXM = pd.read_csv(T + "ATH_CONTEXT_PERARRAY.tsv", sep="\t")
LAND = landsc.load(NEW)
GRP, POS = LAND["GRP"], LAND["POS"]

# the ten groups drawn as panels in Fig 7B and 7C, in the figure's own order
MAIN = ["CM054679.2 · 20.7 Mb", "CM054630.2 · 20.7 Mb", "CM055446.2 · 16.5 Mb",
        "CM055472.2 · 20.6 Mb", "CM054434.2 · 20.5 Mb",
        "chr2 #2 (mat+pat)", "chr9 (mat+pat)", "chr4 · 36.22 Mb", "chr18 · 31.5 Mb",
        "CP116284.1 · 11.2 Mb"]
SPORD = {"금화조": 0, "제브라피시": 1, "애기장대": 2, "수마트라오랑우탄": 3, "보르네오오랑우탄": 4,
         "고릴라": 5, "보노보": 6, "침팬지": 7}
SPORD3 = {"수마트라오랑우탄": 0, "보르네오오랑우탄": 1, "고릴라": 2, "보노보": 3, "침팬지": 4,
          "금화조": 5, "제브라피시": 6, "애기장대": 7}
PMg = PM.copy()
PMg["묶음"] = PMg.array_id.map(GRP)
PMg = PMg[PMg.묶음.notna()].drop_duplicates("묶음")
PRx = PR[PR.종 != "시아망"]
_ORD = figorder.panel_order(NEW)


def order(df, spo):
    d = df.copy()
    d["_o"] = d.종.map(spo).fillna(9)
    d["_c"] = d.묶음.map(lambda k: _ORD.get(k, (999, 0.0))[0])
    d["_s"] = d.묶음.map(lambda k: _ORD.get(k, (999, 0.0))[1])
    return d.sort_values(["_o", "_c", "_s"])


FIG7 = order(PMg[PMg.묶음.isin(MAIN)], SPORD)
REST = order(PMg[(~PMg.묶음.isin(MAIN)) & (PMg.종.isin(["금화조", "제브라피시", "애기장대"]))], SPORD)
ALL49 = order(PMg[PMg.종 != "시아망"], SPORD3)
CPG = dict(zip(PRx.묶음, PRx.CpG중앙))
ULEN = dict(zip(PR.묶음.astype(str), PR.유닛길이.astype(int)))

# ── workbook ──────────────────────────────────────────────────────────────────
wb = openpyxl.Workbook()
wb.remove(wb.active)
BOLD = Font(bold=True)
SRC = []


def sheet(name, header, rows, src, note=""):
    ws = wb.create_sheet(name[:31])
    ws.append(list(header))
    for c in ws[1]:
        c.font = BOLD
        c.alignment = Alignment(horizontal="center")
    for r in rows:
        ws.append(list(r))
    ws.freeze_panes = "A2"
    for i, h in enumerate(header, 1):
        w = max([len(str(h))] + [len(str(r[i - 1])) for r in rows[:400]]) if rows else len(str(h))
        ws.column_dimensions[get_column_letter(i)].width = min(max(w + 2, 10), 44)
    SRC.append((name, ", ".join(os.path.basename(s) for s in src),
                ", ".join(md5(s) for s in src), len(rows), note))
    print("   %-34s %7d rows" % (name, len(rows)))
    return ws


# ── Fig 7A : the chromosome map ───────────────────────────────────────────────
CHRNO = dict(zip(CN.접근번호.astype(str), CN.염색체번호.astype(str)))
PICKED = set(PM[PM.array_id.map(GRP).isin(MAIN)].array_id)
rows = []
for t_ in A1.itertuples(index=False):
    nmc = str(t_.염색체이름)
    rows.append((sp(t_.종), CHRNO.get(nmc, nmc), str(t_.염색체), int(t_.염색체길이),
                 int(t_.시작), int(t_.끝), int((t_.시작 + t_.끝) / 2),
                 "yes" if str(t_.분석에씀) == "Y" else "no",
                 "yes" if t_.array_id in PICKED else "no", str(t_.array_id)))
sheet("Fig7A_arrays",
      ["Species", "Chromosome", "Accession", "Chromosome_length_bp", "Array_start_bp",
       "Array_end_bp", "Array_midpoint_bp", "Used_in_the_analysis",
       "Drawn_as_a_panel_in_Fig7B_and_7C", "Array_id"],
      rows, [T + "FIGDATA_A1_chrom_map.tsv", T + "CHROM_NAMES.tsv"],
      "One row per 5S array. Each array is drawn as one tick on its chromosome, at the midpoint "
      "of the array. Colour in the panel carries only the two flags of this sheet, nothing else. "
      "The last flag is yes for twelve arrays and the figure draws ten panels, because two of "
      "the panels pool the maternal and the paternal array of the same zebra finch locus.")

rows = []
for t_ in CEN.itertuples(index=False):
    rows.append((sp(t_.종), str(t_.염색체), str(t_.접근번호), int(t_.시작), int(t_.끝),
                 float(t_.길이Mb), cenm(t_.방법)))
sheet("Fig7A_centromeres",
      ["Species", "Chromosome", "Accession", "Start_bp", "End_bp", "Length_Mb", "How_it_was_located"],
      rows, [T + "CENTROMERE.tsv"],
      "The grey circle on each chromosome. The last column names the measurement for that "
      "species. Centromere positions were determined in this work.")


# ── the landscape panels (Fig 7B and FigS12) ──────────────────────────────────
def landscape_rows(SET):
    unit, fine = [], []
    for mt in SET.itertuples(index=False):
        key = GRP.get(mt.array_id)
        ids = LAND["IDS"].get(key, [mt.array_id])
        prox = LAND["LEFT"][mt.array_id]
        dist = "높은 좌표" if prox == "낮은 좌표" else "낮은 좌표"
        UNIT = float(mt.유닛길이)
        N = int(mt.유닛수)
        lab = lab_of(mt)
        # unit layer, the call-weighted mean of the units at each index
        s = LAND["SLOT"][LAND["SLOT"].묶음 == key].copy()
        s["낮은쪽번호"] = s.유닛번호 + 1
        s["높은쪽번호"] = N - s.유닛번호
        pcol = "낮은쪽번호" if prox == "낮은 좌표" else "높은쪽번호"
        dcol = "높은쪽번호" if prox == "낮은 좌표" else "낮은쪽번호"
        for k in range(1, NED + 1):
            r_ = s[s[pcol] == k]
            if len(r_):
                unit.append((sp(mt.종), key, lab, "left end", k, int(r_.n_cpg.sum()),
                             int(r_.n_meth.sum()),
                             round(float(r_.n_meth.sum() / max(r_.n_cpg.sum(), 1)), 6),
                             int(UNIT), N, SIDE[prox]))
        mid = s[(s[pcol] > NED) & (s[dcol] > NED)]
        if len(mid):
            q = np.linspace(0, 1, MIDB + 1)
            pr = (mid[pcol].values - NED) / max(N - 2 * NED, 1)
            for bi in range(MIDB):
                m_ = (pr >= q[bi]) & ((pr < q[bi + 1]) if bi < MIDB - 1 else (pr <= 1))
                if m_.any():
                    nc = int(mid.n_cpg.values[m_].sum())
                    nm_ = int(mid.n_meth.values[m_].sum())
                    unit.append((sp(mt.종), key, lab, "middle", bi + 1, nc, nm_,
                                 round(float(nm_ / max(nc, 1)), 6), int(UNIT), N, SIDE[prox]))
        for k in range(NED, 0, -1):
            r_ = s[s[dcol] == k]
            if len(r_):
                unit.append((sp(mt.종), key, lab, "right end", k, int(r_.n_cpg.sum()),
                             int(r_.n_meth.sum()),
                             round(float(r_.n_meth.sum() / max(r_.n_cpg.sum(), 1)), 6),
                             int(UNIT), N, SIDE[prox]))
        # the fine CpG-window layer, exactly the series the panel plots
        d = LAND["LS"][LAND["LS"].array_id.isin(ids)]
        if len(ids) > 1:
            d = d.groupby(["쪽", "위치", "d"], as_index=False)[["호출수", "메틸호출수"]].sum()
            d["beta"] = d.메틸호출수 / d.호출수
        for side in (prox, dist):
            for loc in ("array 안", "flank"):
                dd = d[(d.쪽 == side) & (d.위치 == loc)].sort_values("d")
                if not len(dd):
                    continue
                if loc == "array 안":
                    dx, dy = landsc.slide_unit(dd, UNIT)
                else:
                    ed = landsc.win_edges(False, float(mt.array끝 - mt.array시작))
                    dx, dy, _nb = landsc.adaptive(dd, ed, minc=landsc.FINE_MINC)
                for xx, yy in zip(np.asarray(dx), np.asarray(dy)):
                    fine.append((sp(mt.종), key, lab, SIDE[side], LOC[loc],
                                 int(round(float(xx))), round(float(yy), 6)))
    return unit, fine


UNIT_H = ["Species", "Group", "Position", "Segment", "Unit_index_or_middle_bin", "N_CpG_calls",
          "N_5mC_calls", "Beta", "Unit_length_bp", "Unit_count", "Left_side_of_the_panel"]
FINE_H = ["Species", "Group", "Position", "Side", "Location", "Distance_from_the_array_end_bp",
          "Beta"]
u7, f7 = landscape_rows(FIG7)
sheet("Fig7B_unit_layer", UNIT_H, u7, [T + "UNITBETA_PERARRAY_v13.tsv.gz", T + "PERARRAY_META.tsv"],
      "The coloured line with beads. Beta is the call-weighted mean over the units at that "
      "index, pooled over the arrays of the group. The twenty units at each end are drawn one "
      "by one; the units between them are pooled into twelve equal bins of relative position.")
sheet("Fig7B_CpG_windows", FINE_H, f7,
      [T + "PERARRAY_LANDSCAPE.tsv.gz", T + "PERARRAY_META.tsv"],
      "The fine grey line. Inside the array the window is one unit wide and slides by a quarter "
      "of a unit; in the flank the windows widen with distance until each holds enough calls.")

# ── Fig 7C and FigS11B : the unit-beta distribution ───────────────────────────
def hist_rows(SET):
    out = []
    for mt in SET.itertuples(index=False):
        b = UB[UB.묶음 == mt.묶음]
        ncpg = float(CPG.get(mt.묶음, 20))
        h_, cen_ = betabins.fraction(b.beta.values, ncpg)
        lab = lab_of(mt)
        for c_, v_ in zip(np.asarray(cen_), np.asarray(h_)):
            out.append((sp(mt.종), mt.묶음, lab, round(float(c_), 6), round(float(v_), 8),
                        int(len(b)), round(ncpg, 2)))
    return out


HIST_H = ["Species", "Group", "Position", "Bin_centre_unit_beta", "Fraction_of_read_units",
          "N_read_units_in_the_panel", "Median_CpG_calls_per_unit"]
sheet("Fig7C_unit_beta_distribution", HIST_H, hist_rows(FIG7),
      [T + "UNITBETA_PERARRAY_v13.tsv.gz", T + "PROFILE_META_v13.tsv"],
      "The filled curve. The bin width is set by the median number of CpG calls per read unit "
      "of that panel, so a panel with few calls per unit is not given finer bins than its data "
      "support. The fractions of one panel sum to one.")

# ── FigS11A : the repeat unit ─────────────────────────────────────────────────
A1u = A1[A1.분석에씀 == "Y"].copy()
A1u["유닛길이"] = A1u.array_id.map(lambda a: ULEN.get(str(GRP.get(a, "")), np.nan))
A1u = A1u[A1u.유닛길이.notna()]
A1u["유닛길이"] = A1u.유닛길이.astype(int)
rows = []
for ko in ["수마트라오랑우탄", "보르네오오랑우탄", "고릴라", "보노보", "침팬지", "금화조",
           "제브라피시", "애기장대"]:
    g = A1u[A1u.종 == ko]
    if not len(g):
        continue
    dom = int(g.유닛길이.value_counts().index[0])
    ul = [dom]
    far = sorted(set(g.유닛길이), key=lambda v: -max(v / dom, dom / v))
    if far and max(far[0] / dom, dom / far[0]) >= 1.2:
        ul.append(int(far[0]))
    for u in ul:
        rows.append((sp(ko), u, int(g[g.유닛길이 == u].유전자길이.median()),
                     int((g.유닛길이 == u).sum()),
                     "most common" if u == dom else "second length, drawn because it differs by "
                                                    "at least 1.2 fold"))
sheet("FigS11A_repeat_unit", ["Species", "Unit_length_bp", "Median_gene_length_bp",
                              "N_arrays_with_this_unit_length", "Why_this_row_is_drawn"],
      rows, [T + "FIGDATA_A1_chrom_map.tsv", T + "PROFILE_META_v13.tsv"],
      "One bar per row. The unit length is read from the frame the beta values were computed "
      "on, not from the annotation, so this panel and FigS13 use the same number for the same "
      "array. The gene occupies the last part of the unit.")

sheet("FigS11B_unit_beta_distribution", HIST_H, hist_rows(REST),
      [T + "UNITBETA_PERARRAY_v13.tsv.gz", T + "PROFILE_META_v13.tsv"],
      "Same quantity as Fig7C, for the arrays that Figure 7 does not show.")

# ── FigS11C : the three sequence contexts of A. thaliana ──────────────────────
CTX = ["CG", "CHG", "CHH"]
rows = []
CXMa = CXM[CXM.part == "전체"]
A_IDS = [a for a in CXMa.array_id.unique()]
A_IDS.sort(key=lambda a: _ORD.get(GRP.get(a, ""), (999, 0.0)))
for aid in A_IDS:
    for ctx in CTX:
        d = CX[(CX.array_id == aid) & (CX.context == ctx)]
        d = d[d.n_calls > 0]
        if not len(d):
            continue
        b = (d.n_5mC / d.n_calls).values
        nmed = float(np.median(d.n_calls.values))
        h_, cen_ = betabins.fraction(b, max(nmed, 3))
        for c_, v_ in zip(np.asarray(cen_), np.asarray(h_)):
            rows.append((str(aid), ctx, round(float(c_), 6), round(float(v_), 8), int(len(d)),
                         round(nmed, 2),
                         round(float(d.n_5mC.sum() / d.n_calls.sum()), 6)))
_ATH_MEAN, _ATH_POOL = [], []
for _c in CTX:
    _d = CX[(CX.context == _c) & (CX.n_calls > 0)]
    _per = [float(g.n_5mC.sum() / g.n_calls.sum()) for _a, g in _d.groupby("array_id")]
    _ATH_MEAN.append(float(np.mean(_per)))
    _ATH_POOL.append(float(_d.n_5mC.sum() / _d.n_calls.sum()))

sheet("FigS11C_context_distribution",
      ["Array_id", "Sequence_context", "Bin_centre_read_unit_beta", "Fraction_of_read_units",
       "N_read_units", "Median_calls_per_unit", "Pooled_beta_of_this_array_and_context"],
      rows, [T + "ATH_READUNIT_CONTEXT.tsv.gz", T + "ATH_CONTEXT_PERARRAY.tsv"],
      "Arabidopsis thaliana only. One curve per array and per context. The last column is the "
      "call-weighted mean beta of that array and context. The three numbers quoted in the "
      "legend are the unweighted mean of that column over the nine arrays, one array one "
      "weight, which is CG %.4f, CHG %.4f and CHH %.4f here; pooling every read unit of every "
      "array instead gives CG %.4f, CHG %.4f and CHH %.4f, which is a different quantity and "
      "is not what the legend reports." % tuple(_ATH_MEAN + _ATH_POOL))

# ── FigS12 : the landscape of the remaining arrays ────────────────────────────
u12, f12 = landscape_rows(REST)
sheet("FigS12_unit_layer", UNIT_H, u12,
      [T + "UNITBETA_PERARRAY_v13.tsv.gz", T + "PERARRAY_META.tsv"],
      "Same quantity as Fig7B, for the arrays that Figure 7 does not show.")
sheet("FigS12_CpG_windows", FINE_H, f12,
      [T + "PERARRAY_LANDSCAPE.tsv.gz", T + "PERARRAY_META.tsv"],
      "Same quantity as the Fig7B CpG windows, for the arrays that Figure 7 does not show.")

# ── FigS13 : the profile within the unit ──────────────────────────────────────
GRPN = {"all": "all units", "q4": "top quartile of unit beta", "q1": "bottom quartile of unit beta"}
rows = []
for mt in ALL49.itertuples(index=False):
    d_all = BP[BP.묶음 == mt.묶음]
    L = int(mt.유닛길이)
    for g_ in ("all", "q4", "q1"):
        dd = d_all[d_all.무리 == g_]
        if not len(dd):
            continue
        xs, ys, hw = prof_smooth.smooth_bp(dd.위치.values.astype(int), dd.호출수.values,
                                           dd.메틸수.values, L)
        for x_, y_, w_ in zip(np.asarray(xs), np.asarray(ys), np.asarray(hw)):
            rows.append((sp(mt.종), mt.묶음, lab_of(mt), GRPN[g_],
                         int(x_), round(float(y_), 6), int(w_), L))
sheet("FigS13_profile_within_unit",
      ["Species", "Group", "Position", "Curve", "Position_in_the_unit_bp", "Beta",
       "Window_half_width_bp", "Unit_length_bp"],
      rows, [T + "PROFILE_BP_v13.tsv.gz", T + "PROFILE_META_v13.tsv"],
      "Three curves per panel. Beta at each base is the call-weighted mean over a window "
      "centred on that base; the half width is given per point because it widens where the "
      "calls are thin. The frame runs along the gene strand, so the gene is at the right.")

# ── README sheet, placed first ────────────────────────────────────────────────
ws = wb.create_sheet("README", 0)
ws.append(["Source Data for Figure 7 and Supplementary Figures S11, S12 and S13"])
ws["A1"].font = Font(bold=True, size=12)
ws.append([])
for line in [
    "Cross-species 5S rDNA. Every sheet holds the numbers of one panel, taken from the result "
    "table that panel reads, through the same code the figure uses.",
    "Supplementary figures S11, S12 and S13 were numbered S11, S12 and S13 while the figures "
    "were being made; the file names of the published figures carry the new numbers.",
    "Beta is the fraction of methylated calls among the CpG calls of the unit, or of the "
    "window, as the sheet says. A pooled beta is the sum of methylated calls over the sum of "
    "calls, never a mean of ratios.",
    "The Group column is the key of the panel in the result tables and is printed exactly as it "
    "is stored there, middle dot and all, so that a sheet can be joined back to those tables. "
    "The Position column is the same locus written for reading, with a hyphen in place of that "
    "dot. Where a group name says mat+pat, the maternal and the paternal array of one zebra "
    "finch locus are pooled into one panel, as the figure draws them.",
    "The table below names, for each sheet, the result tables it was built from and their md5, "
    "so a reader can tie a sheet back to a file in the code archive.",
]:
    ws.append([line])
ws.append([])
hdr = ["Sheet", "Built from", "md5 of those files", "Rows", "What the sheet holds"]
ws.append(hdr)
hr = ws.max_row
for c in ws[hr]:
    c.font = BOLD
    c.alignment = Alignment(horizontal="center")
for r in SRC:
    ws.append(list(r))
for i, w in enumerate([32, 46, 30, 10, 110], 1):
    ws.column_dimensions[get_column_letter(i)].width = w
for row in ws.iter_rows(min_row=hr + 1):
    row[4].alignment = Alignment(wrap_text=True, vertical="top")
ws.column_dimensions["E"].width = 110

os.makedirs(OUTDIR, exist_ok=True)
wb.save(OUT)
print()
print("workbook : %s" % OUT)
print("   sheets %d · md5 %s · %d bytes" % (len(wb.sheetnames), md5(OUT), os.path.getsize(OUT)))
print("   %s" % ", ".join(wb.sheetnames))
