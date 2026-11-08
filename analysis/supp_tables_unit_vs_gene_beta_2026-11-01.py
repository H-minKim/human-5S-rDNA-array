# -*- coding: utf-8 -*-
# . The two supplementary tables of the cell line Results, Supplementary Table S3 and Supplementary Table S4.
#   Table 1 (second Results paragraph): the fixed window evidence table, units 1 to 10 from each end against the interior, which is from the 21st unit inward, for whole unit beta against gene body beta.
#        This is the same computation, the same seed (20261029) and the same order as the scratch script unit_vs_gene_beta_table_2026-10-29.py, so it must produce the same numbers as the 10-29 evidence table (asserted).
#   Table 2 (third Results paragraph): whole unit beta and gene body beta by change point zone, that is, the k of Figure 6I and Supplementary Figure S9A.
#        The ruler is the read count weighted mean of unit beta, the same ruler as table 1. Supplementary Figure S9B holds the values of Source Data version 21 and may be on a different ruler, so the two are printed side by side.
#   Inputs 02_tables/P11_K562_UNIT_BETA_v2.tsv, P11_GM12878FS_UNIT_BETA_v2.tsv and lowmeth_changepoint_4haplotypes_2026-10-12.tsv.
#   array_pos 1 is the distal end, so the kth unit from the proximal end is N-k+1.
import io, os, sys, csv, hashlib
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8")
def _env(name, what):
    """Path from an environment variable. No default, so a missing one stops the run."""
    v = os.environ.get(name)
    if not v:
        sys.exit("set %s to %s" % (name, what))
    return v

T = _env("P11_TABLE_DIR", "the folder holding the per-unit beta tables")
SD = _env("P11_SOURCE_DATA_SUPPFIG1",
          "the Source Data workbook of the 2026-10-18 generation whose panel B values this script checks against (SourceData_SuppFig1.xlsx of that generation)")
OUTD = _env("P11_TABLE_OUT_DIR", "the folder to write the two supplementary tables into")
O1 = os.path.join(OUTD, "supp_table_last10unit_window_unit_beta_gene_beta_K562_GM12878_2026-11-01.tsv")
O2 = os.path.join(OUTD, "supp_table_changepoint_by_zone_unit_beta_gene_beta_K562_GM12878_2026-11-01.tsv")
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()[:12]
NS = {("K562", 1): 139, ("K562", 2): 98, ("GM12878", 1): 119, ("GM12878", 2): 91}
HN = {1: "A", 2: "B"}
FN = {"K562": "P11_K562_UNIT_BETA_v2.tsv", "GM12878": "P11_GM12878FS_UNIT_BETA_v2.tsv"}
COLS = (("beta_unit", "whole unit"), ("beta_gene_119", "gene body (last 119 bp)"))
def wm(s, col): return np.average(s[col], weights=s.n_read) if len(s) else np.nan
U = {}
for cell, fn in FN.items():
    d = pd.read_csv(os.path.join(T, fn), sep="\t"); print(fn, md5(os.path.join(T, fn)), len(d), "줄")
    for slot in (1, 2):
        t = d[d.slot == slot].copy(); assert t.array_pos.max() <= NS[(cell, slot)] and t.array_pos.is_unique
        U[(cell, slot)] = t
# ── Table 1: the fixed windows, with the same computation, the same seed and the same order as the 10-29 evidence table
rng = np.random.default_rng(20261029)
def boot_diff(e, m, col):
    out = []
    for _ in range(2000):
        a = e.sample(len(e), replace=True, random_state=rng.integers(1e9)); b = m.sample(len(m), replace=True, random_state=rng.integers(1e9))
        out.append(wm(b, col) - wm(a, col))
    return np.percentile(out, [2.5, 97.5])
R1 = []
for cell in ("K562", "GM12878"):
    for slot in (1, 2):
        t = U[(cell, slot)].copy(); N = NS[(cell, slot)]
        t["w"] = np.where(t.array_pos <= 10, "d", np.where(t.array_pos > N - 10, "p", np.where((t.array_pos > 20) & (t.array_pos <= N - 20), "i", "x")))
        for col, cname in COLS:
            tt = t.dropna(subset=[col]); m = tt[tt.w == "i"]; row = [cell, HN[slot], cname]
            vals = {}
            for end in ("p", "d"):
                e = tt[tt.w == end]; lo, hi = boot_diff(e, m, col)
                vals[end] = (wm(e, col), len(e), wm(m, col) - wm(e, col), lo, hi)
            row += ["%.3f" % vals["p"][0], str(vals["p"][1]), "%.3f" % wm(m, col), str(len(m)), "%.3f" % vals["d"][0], str(vals["d"][1])]
            row += ["%.3f [%.3f, %.3f]" % vals[e][2:] for e in ("p", "d")]
            R1.append(row)
# Matched against the 10-29 evidence table (ledger section 1). The 0.03-0.13 and 0.32-0.60 of the main text come from here.
REF1 = {("K562", "A", "whole unit"): ("0.074 [-0.010, 0.155]", "0.031 [-0.062, 0.108]"), ("K562", "A", "gene body (last 119 bp)"): ("0.596 [0.511, 0.668]", "0.552 [0.450, 0.642]"),
        ("K562", "B", "whole unit"): ("0.073 [-0.039, 0.154]", "0.126 [0.002, 0.258]"), ("K562", "B", "gene body (last 119 bp)"): ("0.321 [0.050, 0.532]", "0.541 [0.406, 0.661]"),
        ("GM12878", "A", "whole unit"): ("-0.010 [-0.039, 0.025]", "0.479 [0.365, 0.546]"), ("GM12878", "A", "gene body (last 119 bp)"): ("0.054 [-0.028, 0.155]", "0.477 [0.359, 0.552]"),
        ("GM12878", "B", "whole unit"): ("0.331 [0.253, 0.423]", "0.544 [0.521, 0.560]"), ("GM12878", "B", "gene body (last 119 bp)"): ("0.451 [0.376, 0.526]", "0.530 [0.463, 0.589]")}
for r in R1: assert tuple(r[-2:]) == REF1[tuple(r[:3])], (r, REF1[tuple(r[:3])])
k562 = [float(r[k].split()[0]) for r in R1 if r[0] == "K562" for k in (-2, -1)]
print("표 1 = 10-29 근거표 (8 줄 기계 대조 같음) · K562 유닛 전체 차 %.2f-%.2f · 유전자부 차 %.2f-%.2f" % (
    min(k562[0:2] + k562[4:6]), max(k562[0:2] + k562[4:6]), min(k562[2:4] + k562[6:8]), max(k562[2:4] + k562[6:8])))
# ── Table 2: the change point zones
KK = {}
for r in csv.DictReader(io.open(os.path.join(T, "lowmeth_changepoint_4haplotypes_2026-10-12.tsv"), encoding="utf-8"), delimiter="\t"):
    cell = "K562" if r["집합"].startswith("K562") else "GM12878"
    slot = 1 if ("hapA" in r["집합"] or "hap1" in r["집합"]) else 2
    KK[(cell, slot, r["끝"])] = int(r["k_hat(해당 틀)"])
assert len(KK) == 8
#   The ruler is the one used by Supplementary Figure S9A and B: the pooled beta of the read-unit CpG calls of the RU table summed within the zone, with the same filters as pooled() in the change point code.
#     Gene body means read-units with at least 3 calls in 2122 to 2240 (and n_cpg >= 10). Whole unit means every call of a read-unit with n_cpg >= 10. A call counts as methylated at 128 or above.
#     The gene body values must equal the S9B source data (version 21) to three decimal places (asserted). If they do not, the table is not written.
import gzip, collections
V = _env("V11PIPE_DIR", "the folder holding the per-read unit methylation archives"); G0, G1 = 2122, 2240
def pooled2(ru):
    P = collections.defaultdict(lambda: [0, 0, 0, 0, 0, 0])   # gene body methylated, calls, read-units, then whole unit methylated, calls, read-units
    with gzip.open(os.path.join(V, ru), "rt", encoding="utf-8") as g:
        for r in csv.DictReader(g, delimiter="\t"):
            if int(r["n_cpg"]) < 10: continue
            cs = [(int(p), int(m)) for p, m in (x.split(":") for x in str(r["calls"]).split(",") if ":" in x)]
            k = (r["slot"], int(r["array_pos"])); w = [m for p, m in cs if G0 <= p <= G1]
            if len(w) >= 3: P[k][0] += sum(1 for m in w if m >= 128); P[k][1] += len(w); P[k][2] += 1
            if cs: P[k][3] += sum(1 for p, m in cs if m >= 128); P[k][4] += len(cs); P[k][5] += 1
    return P
PU = {"K562": pooled2("RU_K562r.tsv.gz"), "GM12878": pooled2("RU_GM12878FS.tsv.gz")}
B = pd.read_excel(SD, sheet_name="b_methylation")
HB = {("K562", "A"): "K562 hapA", ("K562", "B"): "K562 hapB (98-unit frame)", ("GM12878", "A"): "GM12878 hapA", ("GM12878", "B"): "GM12878 hapB"}
R2 = []
for cell in ("K562", "GM12878"):
    for slot in (1, 2):
        N = NS[(cell, slot)]; kp, kd = KK[(cell, slot, "근위")], KK[(cell, slot, "원위")]
        zof = lambda ap: "proximal" if ap > N - kp else ("distal" if ap <= kd else "interior")
        for z, kz in (("proximal", kp), ("interior", N - kp - kd), ("distal", kd)):
            a = np.zeros(6, int); nu = 0
            for (sl, ap), v in PU[cell].items():
                if sl == str(slot) and zof(ap) == z: a += np.array(v); nu += (v[4] > 0)
            gb, wb = a[0] / a[1], a[3] / a[4]
            b = B[(B.haplotype == HB[(cell, HN[slot])]) & (B.zone == z)]; assert len(b) == 1
            print("%-8s %s %-9s k %-3d 유닛 %-3d | 유닛전체 %.3f (call %d · read유닛 %d) | 유전자부 %.3f (call %d) | S9B %.3f" % (cell, HN[slot], z, kz, nu, wb, a[4], a[5], gb, a[1], b.value.iloc[0]))
            assert round(gb, 3) == round(float(b.value.iloc[0]), 3), ("S9B 와 자가 다르다", cell, slot, z)
            R2.append([cell, HN[slot], z, str(kz), str(nu), str(a[5]), "%.3f" % wb, str(a[4]), "%.3f" % gb, str(a[1])])
print("표 2 유전자부 = S9B 원자료 12 칸 셋째 자리까지 같음")
H1 = ["# Supplementary Table S3. Computed 2026-11-01 by supp_tables_unit_vs_gene_beta_2026-11-01.py",
      "# Whole-unit and gene-body CpG methylation (beta) at the array ends against the interior in K562 and GM12878, in fixed windows that do not depend on the zones.",
      "# Windows: proximal and distal = units 1-10 from each end; interior = units from the 21st from each end inward. beta = mean of unit beta weighted by read count.",
      "# Difference = interior minus end, with 95% bootstrap interval (2,000 resamples of units). K562 haplotype A in its 139-unit frame, haplotype B in its 98-unit frame.",
      "\t".join(["cell_line", "haplotype", "beta", "proximal_1-10", "n_units_proximal", "interior", "n_units_interior", "distal_1-10", "n_units_distal", "interior_minus_proximal [95% CI]", "interior_minus_distal [95% CI]"])]
H2 = ["# Supplementary Table S4. Computed 2026-11-01 by supp_tables_unit_vs_gene_beta_2026-11-01.py",
      "# Whole-unit and gene-body CpG methylation (beta) by zone in K562 and GM12878. Zones as in Figure 6I: proximal and distal = units from each end to the change point k (Supplementary Figure S9A), interior = the rest.",
      "# beta = methylated CpG calls over all CpG calls pooled over the zone, from read segments of one unit with at least 10 CpG. Gene body = positions 2,122-2,240, from read segments with at least 3 calls there, as in Supplementary Figure S9A, B (gene-body values equal S9B).",
      "# Whole unit = all CpG calls of the read segment. K562 haplotype A in its 139-unit frame, haplotype B in its 98-unit frame. GM12878 haplotype A has a proximal zone of one unit.",
      "\t".join(["cell_line", "haplotype", "zone", "units_in_zone", "units_with_calls", "read_segments", "whole_unit_beta", "whole_unit_CpG_calls", "gene_body_beta", "gene_body_CpG_calls"])]
if "--write" not in sys.argv: print("DRY"); sys.exit()
os.makedirs(OUTD, exist_ok=True)
for p, H, R in ((O1, H1, R1), (O2, H2, R2)):
    assert not os.path.exists(p)
    io.open(p, "w", encoding="utf-8", newline="\n").write("\n".join(H + ["\t".join(r) for r in R]) + "\n"); print(md5(p), p)
