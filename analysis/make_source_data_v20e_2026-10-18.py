# -*- coding: utf-8 -*-
# The Source Data workbooks written by this script are an earlier
# generation. They are not the workbooks submitted with the paper.
# The submitted Source Data workbooks for Figure 6 and its
# supplementary figures are the later generation recorded in the
# panel source data ledger
# (59_panel_source_data_ledger_2026-11-07).
# This file is kept, and not deleted, so that the earlier generation
# stays reproducible and so that the two generations can be told
# apart by reading the code rather than by guessing from file names.
# v20e (make_source_data_v20e_2026-10-18.py) = v20d plus the n sheet of main version 21, which holds the real counts of the four real markers. Output goes to rawdata_xlsx_edition21_2026-10-18\
# v20d (make_source_data_v20d_2026-10-17.py) = v20c with only the figure version and md5 in the title changed to v12. Output goes to rawdata_xlsx_edition20_v12_2026-10-17\
# v20c (make_source_data_v20c_2026-10-17.py) = v20b plus the new panel letters of supplementary figure 2 (a and b, c to e, f for marker unit methylation, g, h and i) and the f sheet. Output goes to rawdata_xlsx_edition20_v11_2026-10-17\ and no value changes.
# v20b (make_source_data_v20b_2026-10-17.py) = v20 plus the new panel letters of supplementary figure 2 (a and b, c to e, f, g and h). Output goes to rawdata_xlsx_edition20_v10_2026-10-17\ and no value changes.
# . Source Data xlsx version 20: main version 20, supplementary figure 1 (v6) and supplementary figure 2 (v6).
#   Output = 03_figures\00_Fig6_2026-10-09\rawdata_xlsx_edition20_2026-10-17\ and the older 원자료xlsx_2026-10-16, the version 14 bundle source_data, is left untouched.
#   ① Sheets that did not change are carried over value for value from the version 14 bundle source_data, with only the GM12878 hap1 and hap2 letters renamed to hapA and hapB, which is the standard used in this manuscript.
#   ② Newly built sheets: BDP1 added to main p_ChIP_zone, the end 20 unit values of supplementary 1a, the anchored reads of GM12878 hapB in 1a, BDP1 in supplementary 2d, the interaction table including BDP1, and criteria 1 and 2 of 2h.
#   ③ Check: the numbers printed in the figures (14.7, 47.8, 8.1, the criterion 1 and 2 values and so on) are recomputed from the sheet values and asserted.
import io, os, sys, csv, collections
import openpyxl, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
def _env(name, what):
    """Path from an environment variable. No default, so a missing one stops the run."""
    v = os.environ.get(name)
    if not v:
        sys.exit("set %s to %s" % (name, what))
    return v

R0 = _env("P11_WORK_DIR", "the working folder of this run"); C = os.path.join(R0, "01_code"); T = os.path.join(R0, "02_tables"); F6 = os.path.join(R0, r"03_figures\00_Fig6_2026-10-09")
OLD = os.path.join(F6, "for_paper_tidy_2026-10-17", "source_data"); OUT = os.path.join(F6, "rawdata_xlsx_edition21_2026-10-18"); os.makedirs(OUT, exist_ok=True)
NB = os.path.join(T, "node_original_GM_HOR_anchor_2026-10-17")
REL = {"GM12878 hap1": "GM12878 hapA", "GM12878 hap2": "GM12878 hapB"}
def old_sheet(fn, sh):
    ws = openpyxl.load_workbook(os.path.join(OLD, fn), read_only=True)[sh]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    return [[REL.get(v, v) if isinstance(v, str) else v for v in r] for r in rows]
def df_rows(df): return [list(df.columns)] + df.values.tolist()
def write(fn, title, sheets, notes):
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "README"
    for r in [("item", "description"), ("figure", title), ("made", "2026-10-17, make_source_data_v20_2026-10-17.py")] + notes: ws.append(list(r))
    for nm, rows in sheets:
        w = wb.create_sheet(nm[:31])
        for r in rows: w.append(r)
    p = os.path.join(OUT, fn); wb.save(p); print("  %s : 시트 %s" % (fn, ", ".join(n for n, _ in sheets)))
# ---- shared tables
QG = pd.read_csv(os.path.join(T, "K562_QUALGATE_NEWZONE_2026-10-13.tsv"), sep="\t")
ZRB = pd.read_csv(os.path.join(T, "K562_zoneratio_Q20_NEWZONE_withBDP1BRF1_2026-10-14.tsv"), sep="\t")
BTX = pd.read_csv(os.path.join(T, "K562_interaction_NEWZONE_withBDP1BRF1_2026-10-14.tsv"), sep="\t")
ZEN = {"근위": "proximal", "가운데": "interior", "원위": "distal"}
CALL = {"표적특이": "target-specific", "음성대조": "negative control", "판정 불가(CI 가 1 을 문다)": "inconclusive (CI includes 1)"}
def chip_zone_rows():
    old = old_sheet("SourceData_Fig6-3_main.xlsx", "p_ChIP_zone"); head, body = old[0], old[1:]
    assert head == ["target", "zone", "n_marker", "sum_alt_hq", "sum_exp_hq", "ratio"], head
    b = ZRB[ZRB.target == "BDP1"]; add = [["BDP1", ZEN[r.zone], int(r.n_marker), int(r.sum_alt_hq), float(r.sum_exp_hq), float(r.ratio)] for r in b.itertuples()]
    order = {"POLR3A": 0, "TBP": 1, "BDP1": 2, "SNAPC5": 3}
    rows = sorted(body + add, key=lambda r: (order[r[0]], ["proximal", "interior", "distal"].index(r[1])))
    d = {(r[0], r[1]): r for r in rows}
    for t, want in (("POLR3A", 47.8), ("TBP", 27.9), ("BDP1", 8.1), ("SNAPC5", 2.7)):
        got = round(d[(t, "proximal")][5] / d[(t, "interior")][5], 1); assert got == want, (t, got, want)
    return [head] + rows
def interaction_rows():
    rows = [["target", "nc_interaction_OR", "nc_lo", "nc_hi", "proximal_markers", "interior_markers", "call"]]
    for t in ("POLR3A", "TBP", "BDP1", "SNAPC5", "ATAC"):
        r = BTX[BTX.target == t].iloc[0]; rows.append([t, float(r.nc_interaction_OR), float(r.nc_lo), float(r.nc_hi), int(r["근위표지"]), int(r["가운데표지"]), CALL[r["판정"]]])
    return rows
# ---- main
QG = pd.read_csv(os.path.join(T, "K562_QUALGATE_NEWZONE_2026-10-13.tsv"), sep="\t")
_nr = [["marker", "zone", "element", "WGS_alt", "WGS_depth", "WGS_q", "POLR3A_alt", "POLR3A_depth", "POLR3A_fraction", "POLR3A_fraction_over_q", "ticks_WGS_of_26", "ticks_ChIP_of_26"]]
for v_, z_, el_ in (("2091C", "proximal", "5'-flank"), ("2081A", "proximal", "5'-flank"), ("2128A", "interior", "gene"), ("2226T", "interior", "gene")):
    w_ = QG[(QG.label == "WGSNAT") & (QG.variant == v_)].iloc[0]; c_ = QG[(QG.label == "POLR3A") & (QG.variant == v_)].iloc[0]
    q_ = w_.alt_hq / w_.dep_hq; f_ = c_.alt_hq / c_.dep_hq
    _nr.append([v_, z_, el_, int(w_.alt_hq), int(w_.dep_hq), round(q_, 5), int(c_.alt_hq), int(c_.dep_hq), round(f_, 5), round(f_ / q_, 2), max(1, int(round(q_ * 26 * 20))), max(1, int(round(f_ * 26 * 20)))])
assert [r[10] for r in _nr[1:]] == [1, 1, 3, 11] and [r[11] for r in _nr[1:]] == [7, 7, 1, 1], _nr
m_k = old_sheet("SourceData_Fig6-3_main.xlsx", "m_zone_k")
write("SourceData_Fig6-3_main.xlsx", "Fig 6 panels m–p (main 6-3, version 21; English SVG md5 6115955dd483): m haplotypes and markers · n estimator schematic · o ATAC zone composition · p ChIP zone composition (K562)",
      [("m_markers", old_sheet("SourceData_Fig6-3_main.xlsx", "m_markers")), ("m_zone_k", m_k),
       ("n_example_markers", _nr),
       ("o_ATAC_zone", old_sheet("SourceData_Fig6-3_main.xlsx", "o_ATAC_zone")), ("p_ChIP_zone", chip_zone_rows())],
      [("m", "markers per zone and unit (units = slot:unit, slot 1 = haplotype A, 2 = haplotype B); zone change points k (GM12878 haplotypes relabelled A/B from 1/2)"),
       ("n", "version 21: four real K562 markers in the 5'-flank/gene region where ChIP/ATAC reads are dense (proximal 2091C, 2081A; interior 2128A, 2226T); real marker-allele/total counts in WGS and POLR3A ChIP (Q>=20, K562_QUALGATE_NEWZONE_2026-10-13.tsv); ticks = real fraction x 20 of 26 fragments (same factor for WGS and ChIP); fragment lengths not to scale"),
       ("o", "pre-gate counts (one ruler for both cell lines); proximal/interior = 14.7 (K562), 10.0 (GM12878); shares = sum_alt and sum_exp of proximal vs interior"),
       ("p", "base quality Q>=20; BDP1 rows from K562_zoneratio_Q20_NEWZONE_withBDP1BRF1_2026-10-14.tsv (same ruler; POLR3A identical in both tables); proximal/interior POLR3A 47.8, TBP 27.9, BDP1 8.1, SNAPC5 2.7"),
       ("denominators", "K562 = ENCODE K562 WGS; GM12878 = NA12878 WGS (1000 Genomes ERR3239334), both through the ChIP alignment gate")])
# ---- supplementary figure 1
src = io.open(os.path.join(C, "lowmeth_changepoint.py"), encoding="utf-8").read().split("SETS = []")[0]
src = src.replace('sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")', ""); G = {}; exec(src, G)
UK = G["pooled"]("RU_K562r.tsv.gz"); UG = G["pooled"]("RU_GM12878FS.tsv.gz")
G["V"] = os.path.join(T, "node_original_K562u_methylation_2026-10-17"); UKU = G["pooled"]("RU_K562u.tsv.gz")
A = [["haplotype", "end", "unit_from_end", "unit_in_frame", "frame", "methylated_calls", "total_calls", "beta", "reads"]]
for nm, U, sl, N, frame in (("K562 hapA", UK, "1", 139, "139"), ("K562 hapB", UK, "2", 98, "98 (methylation frame)"), ("GM12878 hapA", UG, "1", 119, "119"), ("GM12878 hapB", UG, "2", 91, "91")):
    for end in ("proximal", "distal"):
        for k in range(1, 21):
            ap = (N - k + 1) if end == "proximal" else k
            UU, fr = (UKU, "UL frame (distal numbering same)") if (nm == "K562 hapA" and end == "distal") else (U, frame)
            d = UU.get((sl, ap)); m, c, r = (d[0], d[1], d[2]) if d else (0, 0, 0)
            A.append([nm, end, k, ap, fr, m, c, round(m / c, 4) if c else None, r])
F = collections.defaultdict(lambda: [0, 0, 0]); seen = set(); FB = [["kind", "unit_from_distal_end", "methylated_calls", "total_calls", "beta", "reads"]]
for kind, fn, keep in (("unique units >= 23 (circle)", "hor_anchor_per_read_unit.tsv", lambda a: "근위닻" in a), ("private variants in repeat (diamond)", "privvar_per_read_unit.tsv", lambda a: True)):
    U2 = collections.defaultdict(lambda: [0, 0, 0])
    for i, l in enumerate(io.open(os.path.join(NB, fn), encoding="utf-8")):
        if i == 0: continue
        r, t, a, q, k, n, m = l.rstrip("\n").split("\t")
        if t != "GM12878_h2" or int(n) < 3 or not keep(a) or (r, int(k)) in seen: continue
        seen.add((r, int(k))); U2[int(k)][0] += int(m); U2[int(k)][1] += int(n); U2[int(k)][2] += 1
    for k in range(1, 21):
        if U2[k][1] and not UG.get(("2", k)): FB.append([kind, k, U2[k][0], U2[k][1], round(U2[k][0] / U2[k][1], 4), U2[k][2]])
ac = old_sheet("SourceData_SuppFig_S09.xlsx", "a_c_m6A_largeMSP"); h = ac[0]; mi = h.index("metric")
write("SourceData_SuppFig1.xlsx", "Supplementary Fig. 1 (version 12 = version 7; English SVG md5 59be44c4b0a4): a gene-region methylation near the array ends and zone boundaries · b methylation · c m6A · d large MSP (cell-line Fiber-seq)",
      [("a_units_1-20", A), ("a_zone_k", m_k), ("a_GM12878_hapB_anchored", FB), ("b_methylation", old_sheet("SourceData_SuppFig_S09.xlsx", "b_methylation")),
       ("c_m6A", [h] + [r for r in ac[1:] if r[mi] == "m6a_gene"]), ("d_largeMSP", [h] + [r for r in ac[1:] if r[mi] == "patch_gene"])],
      [("a", "pooled gene-region (2122–2240) CpG calls (ML>=128 = methylated) from placed reads (n_cpg>=10, >=3 window calls); change point k = two-level binomial MLE, 95% interval by profile likelihood widened by φ (a_zone_k)"),
       ("a GM12878 hapB distal", "open symbols = reads placed by alignment anchored on unique units (unit 23 onward) or by private variants in the 8-unit repeat (units 5–22); shown only where placed reads are absent; the boundary k = 10 is unchanged (not re-estimated)"),
       ("b–d", "values from the cell-line Fiber-seq figure (former Supplementary S09), copied value-for-value; GM12878 hap1/hap2 relabelled hapA/hapB; K562 hapB on the earlier 98-unit frame")])
# ---- supplementary figure 2
gmint = old_sheet("SourceData_SuppFig_S06.xlsx", "c_interaction_GM12878"); gmint[0] = ["target", "nc_interaction_OR", "nc_lo", "nc_hi", "proximal_markers", "interior_markers"]
cr = io.open(os.path.join(C, "supple_fig2_criteria.py"), encoding="utf-8").read().split("fig = plt.figure")[0]; g2 = {}; exec(cr, g2)
H2 = [["target", "criterion_1_expected_proximal_ChIP_reads", "criterion_1_pass(>=3)", "criterion_2_expected_proximal_input_reads", "criterion_2_pass(>=1)", "in_main"]]
for t, main in g2["ROWS"]:
    e = g2["EX"][t]; H2.append([t, round(e, 3), "yes" if e >= 3 else "no", round(g2["EIN"], 3), "yes" if g2["EIN"] >= 1 else "no", "yes" if main else "no"])
want = {"POLR3A": 16.43, "TBP": 219.94, "BDP1": 3.47, "SNAPC5": 6.67, "BRF1": 1.10, "GTF3C2": 1.17}
for r in H2[1:]: assert round(r[1], 2) == want[r[0]], r
assert round(g2["EIN"], 2) == 1.09

# ---- v20c: supplementary 2 panel f, the methylation of units that carry a marker, from K562 reads with no assembly
_MK = {r["variant"]: r for r in csv.DictReader(io.open(os.path.join(T, "P11_K562_MARKERS_NEWZONE_2026-10-13.tsv"), encoding="utf-8"), delimiter="\t")}
_ZN = {"근위(합)": "proximal", "가운데": "interior", "원위(합)": "distal"}
FU = [["read", "q_start", "q_end", "strand", "n_call", "n_meth", "beta", "markers", "marker_zones"]]
_cnt = collections.defaultdict(lambda: [0, 0, 0, 0])
for r in csv.DictReader(io.open(os.path.join(T, "node_original_K562_marker_unit_methylation_2026-10-17", "unit_marker_meth.tsv"), encoding="utf-8"), delimiter="\t"):
    n_, m_ = int(r["n_call"]), int(r["n_meth"])
    if n_ < 3: continue
    vs = [x.split(":")[0] for x in r["markers"].split(";") if x]; zs = sorted({_ZN[_MK[v]["zone"]] for v in vs})
    FU.append([r["read"], int(r["q_start"]), int(r["q_end"]), r["strand"], n_, m_, round(m_ / n_, 4), ";".join(vs), ";".join(zs)])
    b = m_ / n_; k = 0 if b < 0.3 else (2 if b >= 0.5 else 1)
    for key in (["all units"] + (["%s markers" % zs[0]] if len(zs) == 1 else [])): _cnt[key][k] += 1; _cnt[key][3] += 1
FS = [["units", "n", "low_beta_lt_0.3", "intermediate", "high_beta_ge_0.5", "low_%", "high_%", "expected"]]
for key, exp in (("all units", "-"), ("proximal markers", "low"), ("distal markers", "low"), ("interior markers", "high")):
    a_ = _cnt[key]; FS.append([key, a_[3], a_[0], a_[1], a_[2], round(100 * a_[0] / a_[3], 1), round(100 * a_[2] / a_[3], 1), exp])
assert [r[1] for r in FS[1:]] == [2595, 64, 149, 887], [r[1] for r in FS[1:]]
write("SourceData_SuppFig2.xlsx", "Supplementary Fig. 2 (version 12; English SVG md5 229f6c408f36): a, b coverage · c–e marker zone exclusivity · f methylation of units carrying markers (K562 reads, no assembly) · g inclusion criteria · h ATAC zone ratios · i ChIP zone ratios (K562)",
      [("a_b_coverage", old_sheet("SourceData_SuppFig_S06.xlsx", "a_b_coverage")), ("c_e_contigs_K562", old_sheet("SourceData_SuppFig_S01.xlsx", "K562")), ("c_e_contigs_GM12878", old_sheet("SourceData_SuppFig_S01.xlsx", "GM12878")), ("f_marker_unit_meth_summary", FS), ("f_marker_unit_meth_units", FU), ("g_criteria", H2),
       ("h_ATAC_zone", old_sheet("SourceData_SuppFig_S06.xlsx", "c_ATAC_zone")), ("h_interaction_GM12878", gmint), ("i_ChIP_zone", chip_zone_rows()), ("h_i_interaction_K562", interaction_rows())],
      [("f", "K562 Fiber-seq reads (K562r, 874 reads) cut into units by alignment to the canonical unit only (no assembly, no placement); unit beta = gene-region (2122–2240) CpG calls with ML>=128 / calls (>=3 calls); low beta<0.3, high beta>=0.5; marker alleles read at the 70 K562 marker positions"),
       ("h", "pre-gate counts; interaction = proximal/interior ratio divided by that of input (paired marker bootstrap, seed 20261008)"),
       ("i", "base quality Q>=20; BDP1 rows from K562_zoneratio_Q20_NEWZONE_withBDP1BRF1_2026-10-14.tsv; interactions from K562_interaction_NEWZONE_withBDP1BRF1_2026-10-14.tsv"),
       ("c–e", "marker-allele read contigs (former Supplementary S01) copied value-for-value"),
       ("g", "expected reads = sum over proximal markers of depth x WGS marker-allele fraction (Q>=20); all six K562 targets share the same pooled K562 input (ENCFF000YRC, ENCFF000YRN, ENCFF000YQC), so criterion 2 has one value (1.09)")])
print("-> %s" % OUT)
