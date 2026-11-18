# -*- coding: utf-8 -*-
# The Source Data workbook for submission.
#   The values are taken from the very table the figure reads; copying them from elsewhere would let them drift from the figure.
#   In English. The README sheet records the source file and its md5.
import io, os, sys, hashlib, collections, statistics
import openpyxl
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter
sys.stdout.reconfigure(encoding="utf-8")

# -- (repo) paths come from the environment --------------------------------
#   The roots are required; the output and table directories fall back to this file's folder.
WORK_ROOT = os.environ["WORK_ROOT"]
DATA_ROOT = os.environ["DATA_ROOT"]
REGISTRY_DIR = os.environ["REGISTRY_DIR"]
TABLE_DIR = os.environ.get(
    "TABLE_DIR",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tables"))
V11 = os.path.join(TABLE_DIR, "V11F4_DATASET_hap_flags_v2.tsv")
REG = os.path.join(REGISTRY_DIR, "cohort_tiers_Fig1_2026-10-08.tsv")
POPC = os.path.join(DATA_ROOT, "02_tables", "POP_COHORT.tsv")
PED = os.path.join(DATA_ROOT, "05_reference",
                   "1kGP_3202_ped_population_2026-09-23.txt")
#   The deposited tables sit side by side in tables/, so the three groups share one directory.
S5 = S6 = S7 = TABLE_DIR
OUT = os.path.join(os.environ.get("SOURCEDATA_DIR",
                                  os.path.join(WORK_ROOT, "04_figures",
                                               "paper_figures_final_2026-10-07")),
                   "SourceData_Fig1_FigS1.xlsx")


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()[:12]


def tsv(p):
    R = [l.rstrip("\n").split("\t") for l in io.open(p, encoding="utf8")]
    return {k: i for i, k in enumerate(R[0])}, R[1:]


wb = openpyxl.Workbook()
wb.remove(wb.active)
SRC = []          # for the README sheet: (sheet, source file, md5, rows)
BOLD = Font(bold=True)


def sheet(name, header, rows, src_path, note=""):
    ws = wb.create_sheet(name[:31])
    ws.append(header)
    for c in ws[1]:
        c.font = BOLD
        c.alignment = Alignment(horizontal="center")
    for r in rows:
        ws.append(list(r))
    ws.freeze_panes = "A2"
    for i, h in enumerate(header, 1):
        w = max(len(str(h)), *(len(str(r[i - 1])) for r in rows)) if rows else len(str(h))
        ws.column_dimensions[get_column_letter(i)].width = min(max(w + 2, 9), 42)
    SRC.append((name, os.path.basename(src_path) if src_path else "(derived)",
                md5(src_path) if src_path else "-", len(rows), note))
    print("   %-26s %5d rows" % (name, len(rows)))


# ── registry ──────────────────────────────────────────────────────────
H, R = tsv(V11)
hp = [r for r in R if r[H["cohort"]] == "HPRC" and r[H["N_units"]].isdigit()]
# CANCER samples are excluded from the main cohort.
# Read from the registry, do not hard-code: 4 haplotypes / 3 samples / 2 people.
CANCER = set(r[H["sample"]].strip().lower() for r in R if r[H["cohort"]] == "cancer")

# Fig 1B: cohort table
Hr, Rr = tsv(REG)
T = dict((r[Hr["name"]], r) for r in Rr)
coh = {}
for key, lab in (("HPRC", "HPRC"), ("GIAB", "GIAB"), ("cancer", "Cancer")):
    rows = [r for r in R if r[H["cohort"]] == key]
    coh[lab] = (len(set(r[H["sample"]] for r in rows)), len(rows))
rows = [
    ("From HPRC samples", "HPRC", "H0", coh["HPRC"][0], coh["HPRC"][1]),
    ("From HPRC samples", "variant reference", "H1",
     int(T["H1"][Hr["individuals"]]), int(T["H1"][Hr["slots"]])),
    ("From HPRC samples", "methylation reference", "H2",
     int(T["H2"][Hr["individuals"]]), int(T["H2"][Hr["slots"]])),
    ("From HPRC samples", "both (H1 and H2)", "H3",
     int(T["H3"][Hr["individuals"]]), int(T["H3"][Hr["slots"]])),
    ("From GIAB samples", "GIAB", "", coh["GIAB"][0], coh["GIAB"][1]),
    ("Both blocks", "Subtotal", "",
     coh["HPRC"][0] + coh["GIAB"][0], coh["HPRC"][1] + coh["GIAB"][1]),
]
sheet("Fig1B_cohort_table", ["Block", "Row", "Tier", "Individuals", "Haplotypes"], rows, V11,
      "H1-H3 are nested subsets of HPRC and are not added to the subtotal. "
      "Cancer references (%d haplotypes / %d samples) exist in the registry but are not shown "
      "in the figure; the cohort table holds two blocks only. The 1kGP (P10) row is not in "
      "this sheet either, for the same reason, that the panel does not draw it."
      % (coh["Cancer"][1], coh["Cancer"][0]))

# Fig 1C: complete units per haplotype (HPRC)
rows = [(r[H["sample"]], r[H["slot"]], int(r[H["N_units"]])) for r in hp]
rows.sort(key=lambda x: (x[0], x[1]))
sheet("Fig1C_units_per_haplotype", ["Sample", "Slot", "Complete_units"], rows, V11,
      "Slot 1 = haplotype A, slot 2 = haplotype B.")
HAP = [x[2] for x in rows]

# Fig 1D: per individual (both haplotypes assembled)
ind = collections.defaultdict(list)
for r in hp:
    ind[r[H["sample"]]].append(int(r[H["N_units"]]))
pair = sorted((s, v) for s, v in ind.items() if len(v) == 2)
rows = [(s, min(v), max(v), sum(v), abs(v[0] - v[1])) for s, v in pair]
sheet("Fig1D_units_per_individual",
      ["Sample", "Units_haplotype_lower", "Units_haplotype_higher",
       "Units_total", "Absolute_difference"], rows, V11,
      "Only individuals with both haplotypes assembled.")
TOT = [x[3] for x in rows]

# Fig 1E: by superpopulation
Hc, Rc = tsv(POPC)
CANON = dict((r[Hc["sample"]], r[Hc["sp"]].strip()) for r in Rc
             if r[Hc["sp"]].strip() not in ("", "NA"))
PED_SP = {}
for i, l in enumerate(io.open(PED, encoding="utf8")):
    f = l.split()
    if i and len(f) >= 7:
        PED_SP[f[1]] = f[6]
lab = {}
for r in hp:
    s = r[H["sample"]]
    if s in CANON:
        lab[s] = CANON[s]
    elif s in PED_SP:
        lab[s] = PED_SP[s]
rows = sorted((lab[r[H["sample"]]], r[H["sample"]], r[H["slot"]], int(r[H["N_units"]]))
              for r in hp if r[H["sample"]] in lab)
sheet("Fig1E_by_superpopulation",
      ["Superpopulation", "Sample", "Slot", "Complete_units"], rows, POPC,
      "Labels from the 1000 Genomes pedigree; 4 individuals have no label and are excluded.")

# ── Supplementary Figure 1 ────────────────────────────────────────────
Hu, Ru = tsv(os.path.join(S6, "FIG2_data_unitcmp.tsv"))
rows = [(r[Hu["sample"]], r[Hu["src"]], float(r[Hu["d_len_unit"]])) for r in Ru]
sheet("FigS1A_length_difference",
      ["Sample", "Source", "Length_difference_units"],
      sorted(rows), os.path.join(S6, "FIG2_data_unitcmp.tsv"),
      "Array reference minus array-spanning read, in units. "
      "The figure clips the axis at +/-2; values outside fall in the overflow bins.")

Hh, Rh = tsv(os.path.join(S6, "FIG2_data_hifi.tsv"))
rows = [(r[Hh["source"]], r[Hh["ref"]], int(r[Hh["slots"]]), float(r[Hh["recall"]]),
         float(r[Hh["precision"]]), int(r[Hh["n_map_var"]])) for r in Rh]
sheet("FigS1E_variant_recall",
      ["Source", "Reference_allele", "Haplotypes", "Recall_percent",
       "Precision_percent", "Mapped_variants"], rows,
      os.path.join(S6, "FIG2_data_hifi.tsv"),
      "The figure shows the MAJ (majority-allele) rows.")

Hd, Rd = tsv(os.path.join(S5, "DC_TABLE.tsv"))
Rd = [r for r in Rd if r[Hd["sample"]].strip().lower() not in CANCER]   # cancer excluded
rows = [(r[Hd["sample"]], r[Hd["slot"]], r[Hd["source"]], int(r[Hd["c_units"]]),
         int(r[Hd["d_longest_read_units"]]), float(r[Hd["d_over_c"]]),
         int(r[Hd["n_reads"]])) for r in Rd]
sheet("FigS1B_longest_read_span",
      ["Sample", "Slot", "Source", "Complete_units", "Units_in_longest_read",
       "Span_ratio", "Reads"], sorted(rows),
      os.path.join(S5, "DC_TABLE.tsv"),
      "5th percentile of Span_ratio = 0.4454, which is where the 0.45 threshold comes from.")

Ht, Rt = tsv(os.path.join(S5, "FIG_data_table.tsv"))
rows = [(r[Ht["sample"]], int(r[Ht["n_hap"]]), int(r[Ht["n_avail"]]), int(r[Ht["n_ok"]]))
        for r in Rt]
sheet("FigS1D_trimming_result",
      ["Sample", "Haplotypes", "Contig_of_correct_length", "Contig_selected"],
      sorted(rows), os.path.join(S5, "FIG_data_table.tsv"),
      "Assembled = sum(Contig_of_correct_length)/sum(Haplotypes); "
      "Selected = sum(Contig_selected)/sum(Contig_of_correct_length).")

Hp2, Rp2 = tsv(os.path.join(S7, "EXAMPLE_PROFILE.tsv"))
rows = [tuple(r) for r in Rp2]
sheet("FigS1I_worked_examples", list(Hp2), rows,
      os.path.join(S7, "EXAMPLE_PROFILE.tsv"),
      "The three panels of (I): positive control, HG00133 haplotype B, HG03270 haplotype A.")

# `NA21133` slot 2 is excluded.
#   At that site the independent read recount returned no reads at all, so it could be neither confirmed nor refuted,
#   and that haplotype attracts reads least well of the fourteen (median local-best 0.571).
#   It is therefore dropped from both numerator and denominator: sites 14 -> 13, haplotypes 13 -> 12.
DROPPED = ("na21133", "2")       # (sample, slot), lower-cased sample
Hs, Rs = tsv(os.path.join(S7, "INTEGRATED_13SITES.tsv"))
Rs = [r for r in Rs
      if (r[Hs["sample"]].strip().lower(), r[Hs["slot"]].strip()) != DROPPED]
# -- Korean cells.
#   They are translated at the point of writing; the input table is left untouched.
#   The class wording follows the summary below, except that the three the summary merges into Mosaic
#     are kept distinct on this sheet; merging them would destroy its purpose.
S1I_EN = {
    "source": {"앵커": "anchor scan"},
    "class": {"read 동의": "True deletion / insertion",
              "모자이크": "Mosaic",
              "반박": "Contradicted, counted as Mosaic in the summary",
              "재현 안 됨": "Not reproduced, counted as Mosaic in the summary"},
    "decision": {"수리 안 함(09-22 모자이크)": "not repaired (09-22, mosaic)",
                 "수리 안 함(09-22)": "not repaired (09-22)",
                 "수리 안 함(재현안됨)": "not repaired (not reproduced)",
                 "수리함 v11-수정1": "repaired in v11-fix1"},
}
_S1I_IX = dict((Hs[c], S1I_EN[c]) for c in S1I_EN if c in Hs)
rows = [tuple(_S1I_IX[j][v] if (j in _S1I_IX and v in _S1I_IX[j]) else v
              for j, v in enumerate(r)) for r in Rs]
sheet("FigS1H_known_sites", list(Hs), rows,
      os.path.join(S7, "INTEGRATED_13SITES.tsv"),
      "%d candidate structural sites in %d haplotypes, flagged by the two tests of this work. "
      "The NA21133 slot 2 site is excluded, from the numerator and the "
      "denominator alike." % (len(rows), len(set((r[Hs["sample"]], r[Hs["slot"]])
                                                 for r in Rs))))

POS_ALL = 56846   # The total number of reference-unit sites; the denominator is written in one place only.
Ha, Ra = tsv(os.path.join(S7, "ALIGNED_POS_RATE.tsv"))
Ra = [r for r in Ra if r[Ha["sample"]].strip().lower() not in CANCER]   # cancer excluded
Ra = [r for r in Ra                                                    # NA21133 slot 2
      if (r[Ha["sample"]].strip().lower(), r[Ha["slot"]].strip()) != DROPPED]
hapset = set((r[Ha["sample"]], r[Ha["slot"]]) for r in Ra)
rows = [("Haplotypes examined", len(hapset)),
        ("Reference unit positions tested", len(Ra))]
cls = collections.Counter()
hapcls = {}
# The class names, in the priority order of S1I_EN above.
#   CLS accepts either spelling, so the sheet is the same whether the table carries the
#   English names of the deposited copy or the Korean names of the upstream copy.
#   A name in neither raises here instead of being counted as zero.
PRI = list(S1I_EN["class"].values())
CLS = dict(S1I_EN["class"])
CLS.update((v, v) for v in PRI)
for r in Rs:
    hapcls.setdefault((r[Hs["sample"]], r[Hs["slot"]]), set()).add(r[Hs["class"]])
# `Contradicted` and `Not reproduced` are merged into `Mosaic`.
#   Three measurements require them to be merged: the decision is the same at all seven sites (no repair),
#   a canonical majority cannot separate the two branches (62.5% against 62.5%, 10 against 10),
#   and merged, the two frames close: sites 6 + 7 = 13, haplotypes 6 + 6 = 12, 579 + 12 = 591.
EN = {"True deletion / insertion": "True deletion / insertion",
      "Mosaic": "Mosaic",
      "Contradicted, counted as Mosaic in the summary": "Mosaic",
      "Not reproduced, counted as Mosaic in the summary": "Mosaic"}
for k, v in hapcls.items():
    seen = set(CLS[c] for c in v)
    for p in PRI:
        if p in seen:
            cls[EN[p]] += 1
            break
for k, v in sorted(cls.items()):
    rows.append((k, v))
rows.append(("Inversion, block move, duplication", 0))   # matching the wording the figure prints
sheet("FigS1H_classification_summary", ["Category", "Haplotypes"], rows,
      os.path.join(S7, "ALIGNED_POS_RATE.tsv"),
      "Counts behind panel (H). The reference set is 612 (HPRC 604 + GIAB 8) out of the "
      "616 rows of the registry, the other 4 being the cancer block. The test reached "
      "%d of 612 haplotypes, and %s of %s reference unit positions were tested = %.1f%%. "
      "Mosaic and refuted are one class, so %d + %d + %d = %d closes."
      % (len(hapset), format(len(Ra), ","), format(POS_ALL, ","),
         len(Ra) / POS_ALL * 100,
         len(hapset) - sum(cls.values()), cls["True deletion / insertion"],
         cls["Mosaic"], len(hapset)))

# ── summary statistics sheet ──────────────────────────────────────────
dc = sorted(float(r[Hd["d_over_c"]]) for r in Rd)
rows = [
    ("Fig 1C", "Complete units per haplotype", "median", round(statistics.median(HAP), 1)),
    ("Fig 1C", "Complete units per haplotype", "minimum", min(HAP)),
    ("Fig 1C", "Complete units per haplotype", "maximum", max(HAP)),
    ("Fig 1C", "Complete units per haplotype", "n (haplotypes)", len(HAP)),
    ("Fig 1D", "Complete units per individual", "median", round(statistics.median(TOT), 1)),
    ("Fig 1D", "Complete units per individual", "minimum", min(TOT)),
    ("Fig 1D", "Complete units per individual", "maximum", max(TOT)),
    ("Fig 1D", "Complete units per individual", "n (individuals)", len(TOT)),
    ("Fig S1B", "Longest-read span ratio", "5th percentile", round(dc[int(len(dc) * 0.05)], 4)),
    ("Fig S1B", "Longest-read span ratio", "n (haplotypes)", len(dc)),
    ("Fig S1B", "Longest-read span ratio", "at or below 0.45",
     sum(1 for x in dc if x <= 0.45)),
]
sheet("Summary_statistics", ["Figure", "Quantity", "Statistic", "Value"], rows, None,
      "Derived from the sheets in this workbook.")

# ── README ────────────────────────────────────────────────────────────
ws = wb.create_sheet("README", 0)
ws.append(["Source data for Figure 1 and Supplementary Figure 1"])
ws["A1"].font = Font(bold=True, size=13)
ws.append([])
ws.append(["Figures", "Fig1_pipeline_EN.pdf / .svg / .png, SFig1_EN and SFig1_KO"])
ws.append([])
ws.append(["Sheet", "Source file", "Source md5 (first 12)", "Rows", "Note"])
for c in ws[7]:
    c.font = BOLD
for s in SRC:
    ws.append(list(s))
for i, w in enumerate((30, 34, 22, 8, 96), 1):
    ws.column_dimensions[get_column_letter(i)].width = w
ws.freeze_panes = "A8"

wb.save(OUT)
print()
print("saved : %s" % OUT)
print("md5 %s · sheets %d" % (md5(OUT), len(wb.sheetnames)))
