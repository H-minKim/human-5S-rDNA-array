# -*- coding: utf-8 -*-
# Fig 4 세트의 결과를 엑셀(원자료 제출용)로 -
#
# 규칙 : 수는 표에서 읽는다(내 기억에서 적지 않는다), 행마다 출처 파일과 md5 를 함께 적는다.
# 말은 영어로(논문에 그대로 낼 수 있게).
import collections
import hashlib
import io
import os
import sys

import openpyxl
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

SP = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SP)
H = _need("RDNA_ARRAY_METH_DIR")
T4 = os.path.join(H, "02_output_tables", "Fig4_2026-10-09")
P10 = os.path.join(H, "02_output_tables", "P10")
ATL = os.path.join(H, "02_output_tables", "transcribed_region_Atlas")
FIGD = os.path.join(H, "03_figures", "00_frame_reference_2026-10-19")
OUT = os.path.join(FIGD, "Fig4_set_SourceData_2026-10-20.xlsx")
NL = chr(10)


def md5(p):
    return hashlib.md5(io.open(p, "rb").read()).hexdigest()[:12] if os.path.exists(p) else "file not found"


def tsv(p):
    with io.open(p, encoding="utf8") as f:
        rows = [L.rstrip(NL).split("\t") for L in f if L.strip()]
    return rows


# ── 배치표에서 칸 ↔ 패널 파일을 기계로 읽는다
_src = io.open(os.path.join(SP, "assemble_v9_2026-10-18.py"), encoding="utf8").read()
_ns = {"__file__": os.path.join(SP, "assemble_v9_2026-10-18.py")}
exec(compile(_src[:_src.index("def fetch(")], "head", "exec"), _ns)
LAY = [("Fig 4", _ns["FIG4"]), ("Suppl. Fig S4", _ns["S41"]),
       ("Suppl. Fig S5", _ns["S42"]), ("Suppl. Fig S6", _ns["S44"])]


def panels(rows):
    out = []
    for row in rows:
        for item, lab in row:
            if item == "COL":
                for f, l in lab:
                    out.append((l, f))
            elif isinstance(item, list):
                out.append((lab, item[0][0] if isinstance(item[0], list) else item[0]))
            else:
                out.append((lab, item))
    return sorted(out, key=lambda t: t[0])


# ── 칸 설명 (영어) - 수는 넣지 않는다. 수는 Key_numbers 에 출처와 함께 둔다.
WHAT = {
 ("Fig 4", "A"): "Read-unit methylation beta along the array axis: 20 proximal units, the interior, 20 distal units",
 ("Fig 4", "B"): "Distribution of unit methylation beta (same cohort as A), drawn vertically next to A",
 ("Fig 4", "C"): "Molecule by unit map for one haplotype; one row per molecule",
 ("Fig 4", "D"): "Mean within-unit profile of the six units that span a boundary",
 ("Fig 4", "E"): "Block length and the spacing between block starts, in units of array elements",
 ("Fig 4", "F"): "Within-unit profile, array units versus dispersed 5S sites, by methylation state",
 ("Fig 4", "G"): "What is transcribed in the 5S rDNA repeat unit",
 ("Fig 4", "H"): "Methylation landscape along the whole 5S rDNA array",
 ("Fig 4", "I"): "Distribution of read-unit methylation beta in the 5S rDNA array",
 ("Fig 4", "J"): "What is transcribed in the D4Z4 4q repeat unit",
 ("Fig 4", "K"): "Methylation landscape along the whole D4Z4 4q array",
 ("Fig 4", "L"): "Distribution of read-unit methylation beta in the D4Z4 4q array",
 ("Fig 4", "M"): "Base-resolution mean profile across a transition, 5S rDNA",
 ("Fig 4", "N"): "Base-resolution mean profile across a transition, D4Z4 4q",
 ("Fig 4", "O"): "Base-resolution mean profile across a transition, D4Z4 10q",
 ("Suppl. Fig S4", "A"): "Unit methylation beta from three data sources, with the valley between the two modes",
 ("Suppl. Fig S4", "B"): "Four single molecules crossing a boundary, base resolution",
 ("Suppl. Fig S4", "C"): "Fitted component means per array element",
 ("Suppl. Fig S4", "D"): "Number of high-methylation elements per unit",
 ("Suppl. Fig S4", "E"): "Anatomy of low-methylation blocks",
 ("Suppl. Fig S4", "F"): "Anatomy of high-methylation blocks",
 ("Suppl. Fig S4", "G"): "Read by unit map, 5S rDNA; one row per molecule",
 ("Suppl. Fig S4", "H"): "Read by unit map, D4Z4 4q; one row per molecule",
 ("Suppl. Fig S4", "I"): "Within-unit profile by state, 5S rDNA",
 ("Suppl. Fig S4", "J"): "Within-unit profile by state, D4Z4 4q",
 ("Suppl. Fig S4", "K"): "Share of one molecule's total rise carried by steps of each size",
 ("Suppl. Fig S5", "A"): "Which part of the array unit each dispersed 5S site carries",
 ("Suppl. Fig S5", "B"): "Per-site distribution of hard-call methylation beta at the dispersed 5S sites",
}
for _f, _p in (("Suppl. Fig S6", x) for x in "ABCDEFGHIJKLMNOPQRSTU"):
    WHAT.setdefault((_f, _p), "")
S44W = {"A": "D4Z4 10q: what is transcribed in the repeat unit",
        "B": "D4Z4 10q: methylation landscape along the whole array",
        "C": "D4Z4 10q: distribution of unit methylation beta",
        "D": "D4Z4 10q: within-unit profile by state",
        "E": "D4Z4 10q: read by unit map; one row per molecule",
        "F": "DXZ4: what is transcribed in the repeat unit",
        "G": "DXZ4: methylation landscape along the array, Xa and Xi with the male X",
        "H": "DXZ4: distribution of unit methylation beta, male versus female",
        "I": "DXZ4: the one woman whose two X's separate",
        "J": "DXZ4: within-unit profile by state",
        "K": "DXZ4: read by unit map; one row per molecule",
        "L": "RNU2: what is transcribed in the repeat unit",
        "M": "RNU2: methylation landscape along the whole array",
        "N": "RNU2: distribution of unit methylation beta",
        "O": "RNU2: within-unit profile by state",
        "P": "RNU2: read by unit map; one row per molecule",
        "Q": "SST1: what is transcribed in the repeat unit",
        "R": "SST1: methylation landscape along the whole array",
        "S": "SST1: distribution of unit methylation beta",
        "T": "SST1: within-unit profile by state",
        "U": "SST1: read by unit map; one row per molecule"}
for k, v in S44W.items():
    WHAT[("Suppl. Fig S6", k)] = v

COH = {"Fig 4": {"A": "H2", "B": "H2", "C": "H2 (one haplotype)", "D": "H2",
                 "E": "H2", "F": "array side + dispersed sites",
                 "G": "annotation", "J": "annotation"},
       "Suppl. Fig S4": {"A": "HPRC + 1kGP + GIAB", "B": "H2"},
       "Suppl. Fig S5": {"A": "dispersed 5S sites", "B": "dispersed 5S sites"}}

wb = openpyxl.Workbook()
HEAD = Font(bold=True, color="FFFFFF")
FILL = PatternFill("solid", fgColor="2F5597")


def sheet(name, header, rows, widths=None, note=None):
    ws = wb.create_sheet(name[:31])
    r0 = 1
    if note:
        ws.cell(1, 1, note).font = Font(italic=True, color="555555")
        r0 = 3
    for j, h in enumerate(header, 1):
        c = ws.cell(r0, j, h)
        c.font = HEAD
        c.fill = FILL
        c.alignment = Alignment(vertical="center", wrap_text=True)
    for i, row in enumerate(rows, r0 + 1):
        for j, v in enumerate(row, 1):
            ws.cell(i, j, v)
    for j, w in enumerate(widths or [], 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = ws.cell(r0 + 1, 1)
    return ws


# ── 1. README
figs = [("Fig 4", "Fig4_main_2026-10-16"),
        ("Suppl. Fig S4", "SuppFig_S4_blocks_and_molecules_2026-10-16"),
        ("Suppl. Fig S5", "SuppFig_S5_dispersed_5S_2026-10-16"),
        ("Suppl. Fig S6", "SuppFig_S6_other_macrosatellites_2026-10-16")]
rows = [["Source data for the 5S rDNA macrosatellite figures", "", "", ""],
        ["Built", "2026-10-20", "", ""],
        ["Figure set (local path; the folder names are Korean)", FIGD, "", ""],
        ["", "", "", ""],
        ["Figure", "file", "md5 (svg)", "md5 (pdf)"]]
for lab, stem in figs:
    rows.append([lab, stem, md5(os.path.join(FIGD, stem + ".svg")),
                 md5(os.path.join(FIGD, stem + ".pdf"))])
rows += [["", "", "", ""],
         ["Rule", "every number in this workbook is read from the file named in its row; "
          "the md5 column is the first 12 hex digits of that file", "", ""],
         ["Methylation beta", "hard call: each CpG counts 0 or 1 at ML >= 128; "
          "a fully unmethylated unit lands at 0", "", ""],
         ["Three beta definitions", "whole unit, window (319-1178 + 1552-2240), element-wise; "
          "all three are cut at 0.508", "", ""]]
ws = sheet("README", ["", "", "", ""], rows, [26, 78, 18, 18])

# ── 2. Panels
prows = []
for lab, rowsdef in LAY:
    for letter, fstem in panels(rowsdef):
        prows.append([lab, letter, fstem, WHAT.get((lab, letter), ""),
                      COH.get(lab, {}).get(letter, "")])
sheet("Panels", ["Figure", "Panel", "Panel file (stem)", "What the panel shows", "Cohort"],
      prows, [18, 8, 44, 86, 26],
      note="One row per panel of the four figures. Panel letters are those drawn in the figure.")

# ── 3. Key numbers (전부 표에서 읽는다)
krows = []
p = os.path.join(T4, "FIG4_V3_NUMBERS.tsv")
for r in tsv(p)[1:]:
    krows.append(["Fig 4", "A, B, E", r[0], r[1],
                  (r[2] + " " + r[3]) if len(r) > 3 else "", os.path.basename(p), md5(p)])
p = os.path.join(T4, "F2_EXT_STATE_LEVELS_2026-10-11.tsv")
hdr = tsv(p)[0]
for r in tsv(p)[1:]:
    for j in range(2, len(r)):
        krows.append(["Fig 4", "M, N, O", "%s %s %s" % (r[0], r[1], hdr[j]), r[j], "",
                      os.path.basename(p), md5(p)])
sheet("Key_numbers",
      ["Figure", "Panel(s)", "Quantity", "Value", "Note", "Source file", "md5"],
      krows, [14, 12, 46, 16, 16, 40, 14],
      note="Numbers that the figure or the Results text quotes. Read from the file in the "
           "last two columns; not retyped by hand.")

# ── 4. Arrays (P10) - G~L 과 S6 의 바탕 수
p = os.path.join(P10, "P10_ARRAYS.tsv")
KO2EN = {"이봉": "bimodal", "단봉": "unimodal"}
rows = [[KO2EN.get(x, x) for x in r] + [os.path.basename(p), md5(p)] for r in tsv(p)[1:]]
sheet("Arrays_P10", tsv(p)[0] + ["Source file", "md5"], rows, [16, 30, 26, 24, 14],
      note="Per-array summary used by Fig 4 G-L and Suppl. Fig S6 "
           "(unit counts, mixture components, Ashman D, cuts).")

# ── 5. Transcripts - 2026-10-20 에 고친 자리
p = os.path.join(ATL, "TRANSCRIBED_FOR_PANEL_A.tsv")
hd = tsv(p)[0]
keep = ["array_our", "feature", "grade", "cons_strand", "canon_segments", "canon_len", "wraps", "note"]
idx = [hd.index(k) for k in keep]
rows = [[r[i] for i in idx] + [os.path.basename(p), md5(p)] for r in tsv(p)[1:]]
sheet("Transcripts", keep + ["Source file", "md5"], rows,
      [14, 30, 8, 10, 26, 10, 8, 80, 34, 14],
      note="Transcribed parts drawn in the 'a' panels. NOTE: rows are PARTS, not transcripts. "
           "DXZ4 has TWO transcripts (one per strand); the plus-strand one is listed as part1 and "
           "part2_next_unit because it runs past the unit end into the next unit. "
           "The figure draws one lane per transcript.")

# ── 6. Atlas summary (전사체 수의 근거)
p = os.path.join(ATL, "ATLAS_TRANSCRIBED_IN_OUR_FRAME.tsv")
if os.path.exists(p):
    hd = tsv(p)[0]
    rows = [r + [os.path.basename(p), md5(p)] for r in tsv(p)[1:]]
    sheet("Transcripts_atlas", hd + ["Source file", "md5"], rows, [14, 28, 8] + [12] * 12,
          note="The Atlas reply counted per transcript, not per drawn part. "
               "This sheet is the evidence that DXZ4 has two transcripts.")

# ── 6b. Cohorts (P10) - ARRAY_COHORT_2026-10-18.txt 의 표를 그대로 읽는다
p = os.path.join(T4, "ARRAY_COHORT_2026-10-18.txt")
rows = []
if os.path.exists(p):
    for L in io.open(p, encoding="utf8"):
        if not L.startswith("| `"):
            continue
        cc = [x.strip() for x in L.strip().strip("|").split("|")]
        if len(cc) < 4:
            continue
        rows.append([cc[0].strip("`"), cc[1], cc[2].replace("", "").replace("**", ""),
                     cc[3].replace("", "").replace("**", ""),
                     os.path.basename(p), md5(p)])
sheet("Cohorts_P10",
      ["Array", "Complete units", "Individuals", "Haplotypes", "Source file", "md5"],
      rows, [14, 16, 14, 14, 34, 14],
      note="P10 cohort (1kGP ONT R10). Gate identical to the figures: n_cpg >= 20 and at least "
           "90% of the common CpG columns covered. The same nine individuals carry all six arrays.")

# ── 6c. Samples of the H2 cohort
p = os.path.join(T4, "SAMPLE_DATASOURCE_2026-10-16.tsv")
if os.path.exists(p):
    t = tsv(p)
    cnt = collections.Counter(r[1] for r in t[1:])
    rows = [[k, v, os.path.basename(p), md5(p)] for k, v in sorted(cnt.items())]
    rows.append(["TOTAL", sum(cnt.values()), os.path.basename(p), md5(p)])
    sheet("Samples_by_source", ["Data source", "Samples", "Source file", "md5"], rows,
          [34, 12, 40, 14],
          note="How many samples come from each data source in the figure set.")

# ── 6d. 뺀 칸 G 의 수 + F 의 네 곡선 - 칸은 뺐어도 수는 남긴다
import subprocess


def _need(k):
    """The deposited copy keeps no default: it stops rather than read the wrong place."""
    v = os.environ.get(k)
    if not v:
        raise SystemExit("environment variable %s is not set (the "
                         "deposited copy has no default; see README.md,"
                         " section 'Environment variables')" % k)
    return v

try:
    _o = subprocess.run([sys.executable, "-X", "utf8",
                         os.path.join(SP, "recheck_panelG.py")],
                        capture_output=True, text=True, encoding="utf8", timeout=300).stdout
except Exception as _e:
    _o = ""
if "| element |" in _o:
    blocks = _o.split("| element |")
    g1 = [l for l in blocks[1].splitlines() if l.startswith("| ") and "---" not in l]
    g2 = [l for l in blocks[2].splitlines() if l.startswith("| ") and "---" not in l]
    rows = []
    for L in g1:
        c = [x.strip() for x in L.strip().strip("|").split("|")]
        rows.append([c[0], "high minus low", c[1], c[2], c[3]])
    rows.append(["", "", "", "", ""])
    for L in g2:
        c = [x.strip() for x in L.strip().strip("|").split("|")]
        rows.append([c[0], "state levels (Fig 4 F curves)",
                     "array low %s / high %s" % (c[1], c[2]),
                     "dispersed low %s / high %s" % (c[3], c[4]), ""])
    sheet("PanelG_removed",
          ["Element", "Quantity", "Array (no HG008-T)", "Dispersed sites", "Gap (array - dispersed)"],
          rows, [14, 30, 34, 34, 24],
          note="The per-element high minus low beta, which no panel draws, and the "
               "element means of the four curves that Figure 4 F draws. Computed from "
               "ARRAYPROF_noHG008T_2026-10-16.tsv (md5 a80b86bc5741) and DISP_PROFILES.tsv "
               "(md5 0599c32be6b0) with the original code path. States are defined by the beta cut "
               "0.508, NOT by two modes within a site.")

# ── 6e. 이름이 같은 두 집합 - 섞이면 없는 수가 된다
DISP = os.path.join(_need("RDNA_METH_BLOCKS_DIR"),
                    "dispersed_vs_array_cut0508", "DISP_SITES_TABLE.tsv")
rows = []
if os.path.exists(DISP):
    t = tsv(DISP)
    hh = {c: t[0].index(c) for c in t[0]}
    sel = [r for r in t[1:] if r[hh["cover_pattern"]] == "NTS1+GA+NTS2+prom+gene"]
    rows.append(["Fig 4 F, dispersed side",
                 "NTS1+GA+NTS2+prom+gene (spacer-1 + GA-rich + spacer-2 + 5'-flank + 5S gene)",
                 "%d of %d" % (len(sel), len(t) - 1),
                 "%s of %s" % (format(sum(int(r[hh["n_loci"]]) for r in sel), ","),
                               format(sum(int(r[hh["n_loci"]]) for r in t[1:]), ",")),
                 os.path.basename(DISP), md5(DISP)])
rows.append(["Suppl. Fig S5 A",
             "the same compartment set, in the earlier 32-site analysis",
             "27 of 32", "504,082 of 580,170",
             "Fig4_supplementary_명세_c137457f.md", "produced by that earlier analysis"])
sheet("Dispersed_sets",
      ["Where", "Group", "Sites", "Loci", "Source file", "md5"], rows,
      [26, 66, 14, 24, 40, 14],
      note="WARNING: the two rows carry the SAME group name but come from DIFFERENT site sets. "
           "Putting one row's numbers in the other panel's caption produces a count that does not "
           "exist. Measured 2026-10-20.  --  Note on the column `their_call` in "
           "DISP_SITES_TABLE.tsv: despite the name it is NOT an outside call. It is the "
           "`bimodal` column of an earlier in-house check "
           "(2026-09-17 ... SITE_BIMODAL.tsv, md5 d0f403178721). Its ashmanD values differ from "
           "the current ones for the same sites (4.252/3.574/2.826 against 3.705/2.903/2.281), "
           "so the two must NOT be read as one scale at two thresholds: they are two runs.")

# ── 7. Source tables
STAB = [("FIG4_V3_NUMBERS.tsv", T4, "Fig 4 A, B, E"),
        ("FIG4A_MEAN_CURVE.tsv", T4, "Fig 4 A"),
        ("FIG4B_DENSITY_HPRC.tsv", T4, "Fig 4 B"),
        ("FIG4C_SIX_UNITS.tsv", T4, "Fig 4 D"),
        ("FIG4C2_BLOCK_SPACING.tsv", T4, "Fig 4 E"),
        ("FIG4C3C4_BLOCK_ANATOMY.tsv", T4, "Suppl. Fig S4 E, F"),
        ("FIG4C3C4_BOUNDARY_POSITION.tsv", T4, "Suppl. Fig S4 E, F"),
        ("F2_EXT_STATE_LEVELS_2026-10-11.tsv", T4, "Fig 4 M, N, O"),
        ("ARRAYPROF_noHG008T_2026-10-16.tsv", T4, "Fig 4 F"),
        ("FIG4D_SINGLE_MOLECULES.tsv", T4, "Suppl. Fig S4 B"),
        ("P10_ARRAYS.tsv", P10, "Fig 4 G-L, Suppl. Fig S6"),
        ("P10_FRAME_ALL.tsv", P10, "Fig 4 G-L, Suppl. Fig S6"),
        ("P10_UNIT_STEPS_ABS.tsv", P10, "Fig 4 D, Suppl. Fig S4 K"),
        ("P10_MODEL_TEST.tsv", P10, "Fig 4 F, M-O"),
        ("P10_SUBUNIT_FRAMES.json", P10, "unit frame for every array"),
        ("TRANSCRIBED_FOR_PANEL_A.tsv", ATL, "the 'a' panels"),
        ("ATLAS_TRANSCRIBED_IN_OUR_FRAME.tsv", ATL, "transcript count")]
rows = []
for fn, d, use in STAB:
    p = os.path.join(d, fn)
    n = (len(tsv(p)) - 1) if (os.path.exists(p) and not fn.endswith(".json")) else ""
    rows.append([fn, os.path.relpath(d, H) + "  (02_output_tables = output tables)", use, n, md5(p),
                 "present" if os.path.exists(p) else "MISSING"])
sheet("Source_tables", ["File", "Folder (under the data root)", "Used by", "Rows", "md5", "State"],
      rows, [42, 34, 34, 10, 14, 12],
      note="The tables the figure code reads. Row counts exclude the header.")

del wb["Sheet"]
wb.save(OUT)
print("냈다 : %s" % OUT)
print("   시트 %s" % ", ".join(wb.sheetnames))
print("   md5 %s" % md5(OUT))
