# -*- coding: utf-8 -*-
"""Source Data for Figure 4 and Supplementary Figures S4, S5, S6 (new numbering of 2026-11-06; internal names
S4, S5, S6). One workbook per figure, English only. Built 2026-11-06.

Every value in a data sheet is read from the file named in the Panels sheet. Two kinds of file:
  result tables   the tables the figures read, each re-built on 2026-11-06 by its producer and found identical
  panel tables    tables the CURRENT drawing scripts write while they draw (the plotted values), taken from the
                  2026-11-06 runs of those scripts and copied to panel_tables/ with their provenance
Tables too large for a sheet (L1_READ_UNIT.parquet, P10_FRAME_ALL.tsv, P10_CPG_ALL.tsv.gz, P10_CPG_ADD3.tsv.gz)
are named with their md5 and where they are deposited, and are not copied into the workbook.
The script stops if any cell of any workbook still holds a Hangul syllable.

usage: python -B make_sourcedata_2026-11-06.py <out folder>"""
import csv, hashlib, io, os, re, shutil, sys
import numpy as np
import openpyxl
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


def _need(k):
    """The deposited copy keeps no default: it stops rather than read the wrong place."""
    v = os.environ.get(k)
    if not v:
        raise SystemExit("environment variable %s is not set (the "
                         "deposited copy has no default; see README.md,"
                         " section 'Environment variables')" % k)
    return v


OUT = os.path.abspath(sys.argv[1])
AU = os.path.dirname(os.path.abspath(__file__))
H = _need("RDNA_ARRAY_METH_DIR")
T2 = os.path.join(H, "02_output_tables")
T4 = os.path.join(T2, "Fig4_2026-10-09")
P10 = os.path.join(T2, "P10")
REPO = _need("RDNA_DEPOSIT_DIR")
DRAW = _need("RDNA_FIGURE_CODE_DIR")
SUB = _need("RDNA_SUBMISSION_FIG_DIR")
D14 = os.path.join(H, "03_figures", "14_edition17_S42_NA21133_set_2026-11-01", "data")
LAYER = os.path.join(_need("RDNA_METH_BLOCKS_DIR"), "layer_v11_fix4_HPRC", "out")
m12 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()[:12]
HANGUL = re.compile("[가-힣]")


def sandbox(key, name):
    """last write of <name> under key <key>"""
    hit = None
    for l in io.open(os.path.join(AU, "logs", key + ".reads.tsv"), encoding="utf-8"):
        r = l.rstrip("\n").split("\t")
        if r[0] == "W" and os.path.basename(r[1]) == name:
            hit = r[2] or r[1]
    if not hit or not os.path.exists(hit):
        sys.exit("STOP no sandbox copy of %s in run %s" % (name, key))
    return hit


# ---------------- panel tables : (name in panel_tables, source path, drawing script, expected body md5 or None)
PT = [
    ("Fig4A_density_matrix.npy", sandbox("f4vmax", "f4a_Mn.npy"), "fig04/f4a_vmax_2026-11-01.py"),
    ("Fig4A_colourbar_upper_limit.txt", sandbox("f4vmax", "f4a_vmax.txt"), "fig04/f4a_vmax_2026-11-01.py"),
    ("FIG4C_SIX_UNITS.tsv", sandbox("f4d", "FIG4C_SIX_UNITS.tsv"), "fig04/fig04_D_six_units_across_boundary.py"),
    ("FIG4C2_UNIT_STEPS.tsv", sandbox("f4d", "FIG4C2_UNIT_STEPS.tsv"), "fig04/fig04_D_six_units_across_boundary.py"),
    ("FIG4C3C4_BLOCK_ANATOMY.tsv", sandbox("f4ef", "FIG4C3C4_BLOCK_ANATOMY.tsv"), "fig04/fig04_E_figS4_E_F_block_length_and_anatomy.py"),
    ("FIG4C3C4_BOUNDARY_POSITION.tsv", sandbox("f4ef", "FIG4C3C4_BOUNDARY_POSITION.tsv"), "fig04/fig04_E_figS4_E_F_block_length_and_anatomy.py"),
    ("FIG4C2_BLOCK_SPACING.tsv", sandbox("f4ef", "FIG4C2_BLOCK_SPACING.tsv"), "fig04/fig04_E_figS4_E_F_block_length_and_anatomy.py"),
    ("FIG4_V3_NUMBERS.tsv", sandbox("f4abc", "FIG4_V3_NUMBERS.tsv"), "fig04/fig04_A_B_C_array_axis_density_and_molecule_map.py"),
    ("LAND_SEGMENTS.tsv", sandbox("perarray", "LAND_SEGMENTS.tsv"), "common/fig04_figS4/per_array_panels_fig04_G-L_figS4_G-J_figS6_A-U.py"),
    ("FIG4D_SINGLE_MOLECULES.tsv", sandbox("s41b", "FIG4D_SINGLE_MOLECULES.tsv"), "figS4/figS4_B_four_molecules_through_boundary.py"),
    ("FIG4D_SELECTED_MOLECULES.tsv", sandbox("s41b", "FIG4D_SINGLE_MOLECULES_선택.tsv"), "figS4/figS4_B_four_molecules_through_boundary.py"),
    ("F2_EXT_STATE_LEVELS.tsv", sandbox("mno", "F2_EXT_STATE_LEVELS.tsv"), "fig04/fig04_M_N_O_base_resolution_transition.py"),
    ("F2_EXT_SUMMARY.tsv", sandbox("mno", "F2_EXT_SUMMARY.tsv"), "fig04/fig04_M_N_O_base_resolution_transition.py"),
    ("FigS5_B_persite_tone_points.tsv", sandbox("s42b", "FigS5_B_persite_tone_5col_points.tsv"), "figS5/figS5_B_persite_tone_five_columns_dedup_2026-11-01.py"),
]
# ---------------- result tables : (name, path, producer, expected body md5 measured 2026-11-06)
RT = {
    "ARRAYPROF_noHG008T_2026-10-16.tsv": (os.path.join(T4, "ARRAYPROF_noHG008T_2026-10-16.tsv"), "analysis/arrayprof_noHG008T_2026-10-16.py", "a80b86bc5741"),
    "DISP_PROFILES_S42SET29_2026-10-31.tsv": (os.path.join(T4, "DISP_PROFILES_S42SET29_2026-10-31.tsv"), "analysis/f4g_profile_2026-10-31.py", "80fd540b6c55"),
    "TRANSCRIBED_FOR_PANEL_A.tsv": (os.path.join(T2, "transcribed_region_Atlas", "TRANSCRIBED_FOR_PANEL_A.tsv"), "analysis/atlas2_map_to_our_frame.py, then two recorded steps", "0eab96ae032d"),
    "P10_ARRAY_BOUNDS.tsv": (os.path.join(P10, "P10_ARRAY_BOUNDS.tsv"), "analysis/p10_20_array_bounds.py", "7c92f790f85c"),
    "P10_ARRAYS.tsv": (os.path.join(P10, "P10_ARRAYS.tsv"), "analysis/p10_10_arrays.py", "4299508a8ca8"),
    "FIGX19_E_HAP_cut0508.tsv": (os.path.join(REPO, "analysis", "producers_v11_2026-10-11", "inputs", "FIGX19_E_HAP_cut0508.tsv"), "analysis/producers_v11_2026-10-11/figx19_e_haplevel.py", "5841963d4c7a"),
    "FIGX19_F_SEG.tsv": (os.path.join(REPO, "analysis", "producers_v11_2026-10-11", "inputs", "FIGX19_F_SEG.tsv"), "analysis/producers_v11_2026-10-11/figx19_prep.py", "81ac4e84a389"),
    "FIGX25_BIMODAL.tsv": (os.path.join(REPO, "analysis", "producers_v11_2026-10-11", "inputs", "FIGX25_BIMODAL.tsv"), "analysis/producers_v11_2026-10-11/figx25_shape_prep.py", "b8299d203b35"),
    "FIGX25_DENSITY.tsv": (os.path.join(REPO, "analysis", "producers_v11_2026-10-11", "inputs", "FIGX25_DENSITY.tsv"), "analysis/producers_v11_2026-10-11/figx25_shape_prep.py", "92b829838fe4"),
    "SUBUNIT_BIMODAL2.tsv": (os.path.join(T2, "5S_v11_fix4_HPRC", "SUBUNIT_BIMODAL2.tsv"), "subunit_bimodality2.py run on the HPRC layer", "24c5cc0b9ea9"),
    "CONCORDANCE_SUMMARY.tsv": (os.path.join(LAYER, "CONCORDANCE_SUMMARY.tsv"), "unit_concordance_and_shape.py run on the HPRC layer", "df001fdd3ab0"),
    "STEPS_UNITLEVEL_ALL6.tsv": (os.path.join(T2, "array_steps_extended", "STEPS_UNITLEVEL_ALL6.tsv"), "analysis/p10_19c_steps_unitlevel.py", "cb9fd83238f6"),
    "TRANSITION_ZONE_UNITDIST.tsv": (os.path.join(T2, "array_steps_extended", "TRANSITION_ZONE_UNITDIST.tsv"), "analysis/p10_19c_steps_unitlevel.py", "096053d42c92"),
    "MODEL_TEST_ALL6.tsv": (os.path.join(T2, "claim_test_MvsS", "MODEL_TEST_ALL6.tsv"), "analysis/p10_17b_model_test_dxz4.py", "ae81acc0f287"),
    "P10_UNIT_STEPS_ABS.tsv": (os.path.join(P10, "P10_UNIT_STEPS_ABS.tsv"), "p10_19_unit_steps.py", "f95b9e89a86c"),
    "PROB_A4_SITE.tsv": (os.path.join(D14, "PROB_A4_SITE.tsv"), "analysis/s42_dedup_tables_2026-11-01.py", "926386189854"),
    "PROB_A4_SITE_TONE.tsv": (os.path.join(D14, "PROB_A4_SITE_TONE.tsv"), "analysis/s42_dedup_tables_2026-11-01.py", "0f6dc4055dbf"),
    "P10_DXZ4_POLARITY_FINAL.tsv": (os.path.join(P10, "P10_DXZ4_POLARITY_FINAL.tsv"), "no producer preserved", "a6a120f498ee"),
}
BIG = {  # not copied into a sheet
    "L1_READ_UNIT.parquet": ("983ea5625d32 (raw)", "upstream/v13_L1_read_unit.py", "Zenodo record of the code archive (tables/README.md)"),
    "P10_FRAME_ALL.tsv": ("fdd03601f3cb (raw)", "analysis/p10_08_extract_frame.py", "Zenodo record of the code archive (tables/README.md)"),
    "P10_CPG_ALL.tsv.gz": ("4f9eea48a4c5 (raw)", "analysis/p10_12_extract_cpg.py", "Zenodo record of the code archive (tables/README.md)"),
    "P10_CPG_ADD3.tsv.gz": ("c11682f9e965 (raw)", "analysis/p10_12b_extract_cpg_DXZ4_RNU2_SST1.py", "Zenodo record of the code archive (tables/README.md)"),
}
TR = {  # the only non-English strings in these tables
    "색조 낮음": "lowest flank-tone third", "색조 높음": "highest flank-tone third",
    "이봉": "bimodal", "단봉": "unimodal",
    "하드콜": "hard call", "확률평균": "probability mean",
    "전체": "all", "중간": "middle flank-tone third",
}


def body(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for line in f:
            if not line.startswith(b"#"):
                h.update(line.replace(b"\r\n", b"\n"))
    return h.hexdigest()[:12]


bad = 0
for n, (p, prod, exp) in RT.items():
    if not os.path.exists(p) or body(p) != exp:
        print("STOP result table", n, body(p) if os.path.exists(p) else "missing", "expected", exp); bad += 1
if bad:
    sys.exit(1)

os.makedirs(os.path.join(OUT, "panel_tables"), exist_ok=False)
prov = io.open(os.path.join(OUT, "panel_tables", "00_provenance.tsv"), "w", encoding="utf-8", newline="\n")
prov.write("file\tmd5\tdrawing script (not in this deposit)\tscript md5\tnote\n")
PTP = {}
for n, src, scr in PT:
    q = os.path.join(OUT, "panel_tables", n)
    shutil.copyfile(src, q)
    PTP[n] = (q, scr)
    prov.write("\t".join([n, m12(q), scr, m12(os.path.join(DRAW, scr.replace("/", os.sep))), "yes"]) + "\n")
prov.close()


def tsv_rows(p):
    with io.open(p, encoding="utf-8", newline="") as f:
        return [r for r in csv.reader((l for l in f if not l.startswith("#")), delimiter="\t") if r]


def cell(v, text=False):
    """text=True keeps the value as written (ids and md5 strings: 6963902e6546 would read as a number)"""
    v = TR.get(v, v)
    if isinstance(v, str):
        for k, e in TR.items():
            v = v.replace(k, e)
        if text:
            return v
        try:
            if re.fullmatch(r"-?\d+", v) and not (len(v) > 1 and v.lstrip("-").startswith("0")):
                return int(v)
            if re.fullmatch(r"-?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?", v):
                f = float(v)
                return f if np.isfinite(f) else v
        except ValueError:
            return v
    return v


def textcol(h):
    h = str(h).lower()
    return "md5" in h or h.endswith("_id") or h in ("id", "read_id", "sample", "site", "chr", "slot")


HEAD, FILL = Font(bold=True, color="FFFFFF"), PatternFill("solid", fgColor="2F5597")


def sheet(wb, name, header, rows, note=None):
    ws = wb.create_sheet(name[:31])
    r0 = 1
    if note:
        ws.cell(1, 1, note).font = Font(italic=True, color="555555"); r0 = 3
    for j, h in enumerate(header, 1):
        c = ws.cell(r0, j, cell(h)); c.font = HEAD; c.fill = FILL; c.alignment = Alignment(wrap_text=True, vertical="top")
    keep = [textcol(h) for h in header]
    for i, row in enumerate(rows, r0 + 1):
        for j, v in enumerate(row, 1):
            ws.cell(i, j, cell(v, keep[j - 1] if j - 1 < len(keep) else False))
    for j in range(1, len(header) + 1):
        ws.column_dimensions[get_column_letter(j)].width = 18
    ws.freeze_panes = ws.cell(r0 + 1, 1)


def table_sheet(wb, sname, fname, note):
    if fname in RT:
        p = RT[fname][0]
    else:
        p = PTP[fname][0]
    if fname.endswith(".npy"):
        M = np.load(p)
        cols = ["proximal unit %d" % i for i in range(1, 21)] + ["interior bin %d of 16" % i for i in range(1, 17)] + \
               ["distal unit %d" % i for i in range(20, 0, -1)]
        rows = [["%.2f-%.2f" % (k / 50.0, (k + 1) / 50.0)] + ["%.6f" % x for x in M[k]] for k in range(M.shape[0])]
        sheet(wb, sname, ["beta bin"] + cols, rows, note)
        return
    if fname.endswith(".txt"):
        sheet(wb, sname, ["quantity", "value"], [["colour bar upper limit (99th percentile of the non-zero cells of A)",
                                                    io.open(p, encoding="utf-8").read().strip()]], note)
        return
    t = tsv_rows(p)
    sheet(wb, sname, t[0], t[1:], note)


def src_of(fname):
    if fname in RT:
        p, prod, _ = RT[fname]
        return fname, "result table", prod, m12(p)
    if fname in PTP:
        p, scr = PTP[fname]
        return fname, "panel table (plotted values)", scr, m12(p)
    md, prod, where = BIG[fname]
    return fname, "large table, not copied (" + where + ")", prod, md


# ---------------- figures : panel -> (title from the legend, [(sheet, file, note)], other sources)
FIGS = {
 "Fig4": ("Figure 4", "Fig4.pdf", [
   ("A", "Read-unit methylation along the array", [("A_density", "Fig4A_density_matrix.npy", "Rows are beta bins of 0.02. Columns are the array axis of panel A. Each column is normalised to sum to 1 (share of that column's read-units)."), ("A_colour_bar", "Fig4A_colourbar_upper_limit.txt", None)], ["L1_READ_UNIT.parquet"]),
   ("B", "Distribution of unit methylation", [("B_numbers", "FIG4_V3_NUMBERS.tsv", "Peak positions and other numbers that panels A, B and E print.")], ["L1_READ_UNIT.parquet"]),
   ("C", "Unit states along single molecules", [], ["L1_READ_UNIT.parquet"]),
   ("D", "Mean methylation across a boundary", [("D_six_units", "FIG4C_SIX_UNITS.tsv", "Mean hard-call beta in 50 bp bins over the six units of each boundary window."), ("D_unit_steps", "FIG4C2_UNIT_STEPS.tsv", None)], ["BOUNDARY_WIDE_H2_508.tsv"]),
   ("E", "Ends of low-methylation blocks", [("EF_block_anatomy", "FIG4C3C4_BLOCK_ANATOMY.tsv", None), ("EF_boundary_position", "FIG4C3C4_BOUNDARY_POSITION.tsv", None)], ["L1_READ_UNIT.parquet"]),
   ("F", "Ends of high-methylation blocks", [], ["L1_READ_UNIT.parquet"]),
   ("G", "Within-unit methylation of array units and dispersed 5S loci", [("G_array", "ARRAYPROF_noHG008T_2026-10-16.tsv", None), ("G_dispersed", "DISP_PROFILES_S42SET29_2026-10-31.tsv", None)], []),
   ("H", "Transcribed parts of the 5S rDNA and D4Z4 4q repeat units", [("H_transcribed", "TRANSCRIBED_FOR_PANEL_A.tsv", "Rows are transcribed parts. The figure uses the rows of 5S rDNA and D4Z4 4q.")], []),
   ("I", "Methylation along the whole array in 5S rDNA and D4Z4 4q", [("I_segments", "LAND_SEGMENTS.tsv", None), ("I_array_bounds", "P10_ARRAY_BOUNDS.tsv", None)], ["P10_FRAME_ALL.tsv"]),
   ("J", "Unit methylation distribution in 5S rDNA and D4Z4 4q", [("J_arrays", "P10_ARRAYS.tsv", "Per-array summary. hist50 is the share of units in 50 beta bins of 0.02.")], []),
 ]),
 "FigS4": ("Supplementary Figure S4 (internal name S4)", "FigS4.pdf", [
   ("A", "Unit methylation distribution in three datasets", [("A_density", "FIGX25_DENSITY.tsv", None), ("A_bimodality", "FIGX25_BIMODAL.tsv", None), ("A_haplotype_level", "FIGX19_E_HAP_cut0508.tsv", None), ("A_gene_body_variant", "FIGX19_F_SEG.tsv", None)], []),
   ("B", "Single molecules through a boundary", [("B_molecules", "FIG4D_SINGLE_MOLECULES.tsv", "Sliding mean of nine CpG hard calls (roll9) along each drawn molecule."), ("B_selection", "FIG4D_SELECTED_MOLECULES.tsv", None)], ["BOUNDARY_WIDE_H2_508.tsv"]),
   ("C", "Fitted state levels by sub-unit element", [("C_fit", "SUBUNIT_BIMODAL2.tsv", "Rows with kind = fit are drawn (p_low and p_high).")], []),
   ("D", "Number of high-methylation elements per unit", [("D_concordance", "CONCORDANCE_SUMMARY.tsv", "Rows nh_frac give the share of units with 0 to 6 high-methylation elements.")], []),
   ("E", "Spacing between block starts of the same state", [("E_block_spacing", "FIG4C2_BLOCK_SPACING.tsv", None)], ["L1_READ_UNIT.parquet"]),
   ("F", "Example read-by-unit map of 5S rDNA", [], ["P10_FRAME_ALL.tsv"]),
   ("G", "Example read-by-unit map of D4Z4 4q", [], ["P10_FRAME_ALL.tsv"]),
   ("H", "Transition at the proximal end of 5S rDNA at base resolution", [("HIJ_state_levels", "F2_EXT_STATE_LEVELS.tsv", None), ("HIJ_summary", "F2_EXT_SUMMARY.tsv", None)], ["P10_CPG_ALL.tsv.gz", "P10_CPG_ADD3.tsv.gz"]),
   ("I", "Transition at the distal end of 5S rDNA at base resolution", [], ["P10_CPG_ALL.tsv.gz"]),
   ("J", "Transition at the proximal end of D4Z4 4q at base resolution", [], ["P10_CPG_ADD3.tsv.gz"]),
   ("K", "Within-unit methylation of the two unit classes in 5S rDNA", [], ["P10_CPG_ALL.tsv.gz"]),
   ("L", "Within-unit methylation of the two unit classes in D4Z4 4q", [], ["P10_CPG_ADD3.tsv.gz"]),
   ("M", "Decomposition of each molecule's total rise by step size", [("M_steps", "STEPS_UNITLEVEL_ALL6.tsv", None), ("M_transition_zone", "TRANSITION_ZONE_UNITDIST.tsv", None), ("M_model_test", "MODEL_TEST_ALL6.tsv", None), ("M_unit_steps", "P10_UNIT_STEPS_ABS.tsv", None)], ["P10_FRAME_ALL.tsv"]),
 ]),
 "FigS7": ("Supplementary Figure S7 (internal name S5)", "FigS7.pdf", [
   ("A", "Element content of the dispersed 5S copies", [("A_sites", "PROB_A4_SITE.tsv", "29 sites that carry at least three elements, NA21133 counted once.")], []),
   ("B", "Locus methylation at each dispersed site by flank tone", [("B_points", "FigS5_B_persite_tone_points.tsv", "Plotted values. loci_frac is the share of a site's loci in each beta bin, for the lowest and the highest flank-tone third."), ("B_site_tone", "PROB_A4_SITE_TONE.tsv", None)], []),
 ]),
 "FigS6": ("Supplementary Figure S6 (internal name S6)", "FigS6.pdf", [
   ("A F L Q", "Transcribed parts of the repeat unit (D4Z4 10q, DXZ4, RNU2, SST1)", [("AFLQ_transcribed", "TRANSCRIBED_FOR_PANEL_A.tsv", None)], []),
   ("B G M R", "Methylation along the array", [("BGMR_segments", "LAND_SEGMENTS.tsv", None), ("BGMR_array_bounds", "P10_ARRAY_BOUNDS.tsv", None), ("G_DXZ4_polarity", "P10_DXZ4_POLARITY_FINAL.tsv", None)], ["P10_FRAME_ALL.tsv"]),
   ("C H I N S", "Unit methylation distribution", [("CHINS_arrays", "P10_ARRAYS.tsv", "Per-array summary. hist50 is the share of units in 50 beta bins of 0.02.")], []),
   ("D J O T", "Within-unit methylation", [], ["P10_CPG_ALL.tsv.gz", "P10_CPG_ADD3.tsv.gz"]),
   ("E K P U", "Read-by-unit map", [], ["P10_FRAME_ALL.tsv"]),
 ]),
}
RULES = [
    ["Rule", "Every value in this workbook is read from the file named in the Panels sheet. md5 is the first 12 hex digits of that file."],
    ["Methylation beta", "Hard call. Each CpG counts 1 when the methylation likelihood is at least 128 of 255 and 0 otherwise."],
    ["State cut", "Unit methylation is split at beta 0.508 in the HPRC cohort (H2). The 1kGP set (cohort P10) uses 0.4611."],
    ["Large tables", "Tables listed as large are not copied into this workbook. They are deposited with the code archive."],
]
out_md5 = []
for key, (title, pdf, panels) in FIGS.items():
    wb = openpyxl.Workbook()
    pdfp = os.path.join(SUB, "02_보충그림" if key != "Fig4" else "01_본문그림", pdf)
    sheet(wb, "README", ["Item", "Value"], [["Figure", title], ["Figure file", pdf], ["Figure md5", m12(pdfp)],
                                           ["Built", "2026-11-06"]] + RULES)
    prow, done = [], {}
    for letter, ptitle, sheets, others in panels:
        for sname, fname, note in sheets:
            f, kind, prod, md = src_of(fname)
            prow.append([letter, ptitle, sname, f, kind, prod, md])
            if sname not in done:
                table_sheet(wb, sname, fname, note); done[sname] = 1
        for fname in others:
            f, kind, prod, md = src_of(fname) if fname in BIG or fname in RT else (fname, "result table, see the Figure 4 workbook", "-", "-")
            if fname == "BOUNDARY_WIDE_H2_508.tsv":
                f, kind, prod, md = fname, "result table, not copied (deposited in tables/)", "analysis/boundary_profile_H2_2026-10-11.py", "963331f95ce8"
            prow.append([letter, ptitle, "-", f, kind, prod, md])
    sheet(wb, "Panels", ["Panel", "What the panel shows", "Sheet", "Source file", "Kind", "Producer or drawing script", "md5"], prow,
          "One row per panel and source file. Panels without a sheet are drawn from a large table named in the row.")
    del wb["Sheet"]
    wb.move_sheet("Panels", offset=-(len(wb.sheetnames) - 2))
    hits = []
    for ws in wb.worksheets:
        if HANGUL.search(ws.title):
            hits.append(("sheet name", ws.title))
        for row in ws.iter_rows(values_only=True):
            for v in row:
                if isinstance(v, str) and HANGUL.search(v):
                    hits.append((ws.title, v[:60]))
    q = os.path.join(OUT, "SourceData_%s.xlsx" % key)
    if hits:
        sys.exit("STOP Hangul left in %s : %s" % (key, hits[:8]))
    wb.save(q)
    out_md5.append((os.path.basename(q), m12(q), len(wb.sheetnames)))
for n, md, k in out_md5:
    print("%-24s %s  sheets %d" % (n, md, k))
