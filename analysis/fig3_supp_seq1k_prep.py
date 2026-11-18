# -*- coding: utf-8 -*-
# Re-bins the CpG line of supplementary panel A onto a single window width of 1,000 bp.
#
#   why
#     With a window width that varies with position (50 bp over 0-2 kb, 200 bp over 2-8 kb
#     and 1,000 bp beyond 8 kb), the variance looks larger near the junction.
#     That is a difference of resolution, not of sequence: measuring a perfectly uniform synthetic array,
#     20 copies of the canonical unit concatenated, through the same windows reproduces the same pattern,
#     200 at 50 bp, 95 at 200 bp and 15 at 1 kb. The density itself, the median, is the same in every interval.
#   Drawing one line through one window width removes that artefact, and what the panel has to say
#     (the array carries twice the flank, and the proximal and distal ends are nearly symmetric) survives at 1 kb.
#
#   Nothing is recomputed: bp, n_cpg and n_alu from JUNCTION_SEQ_ANNOT.tsv are simply re-summed into 1 kb bins.
#     The coordinate u, the distance from the nearest complete unit boundary, uses the same expression as fig1_v5_prep.py.
#   One limitation: each original bin (50 or 200 bp) goes whole into the 1 kb bin its start falls in.
#     Bins straddling a boundary are not split; this is the same approximation the current table makes, and the bp total is preserved.
import io, os, sys
import numpy as np, pandas as pd

# The data roots come from environment variables.
#   They have no default: the script stops with a message if one is not set.
DATA_ROOT = os.environ["DATA_ROOT"]
V_ROOT = os.environ["V_ROOT"]
S = DATA_ROOT + "/"
SUMF = V_ROOT + "/03_maps/g8canonv10/CONSENSUS_SUM_G8_full.tsv"
# Where the tables are read and written. The tables are not deposited next to the code, so these paths come from the environment.
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
W = 1000.0
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

D = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "data_source", "gate"])
D = D[(D.gate == 1) & (D.data_source.astype(str) == "HPRC_ONT_submissions")]
HP = set(zip(D["sample"].astype(str), D.slot.astype(str)))
print("[cohort] haplotypes %d" % len(HP))

SUM = pd.read_csv(SUMF, sep="\t", dtype={"slot": str})
SUM["off_원위"] = SUM.off_prox
SUM["off_근위"] = SUM.off_dist
off = {}
for r in SUM.itertuples(index=False):
    off[(r.sample, str(r.slot), "원위")] = float(r.off_원위)
    off[(r.sample, str(r.slot), "근위")] = float(r.off_근위)

A = pd.read_csv(S + "02_tables/JUNCTION_SEQ_ANNOT_V145.tsv", sep="\t", dtype={"slot": str})
#   The input is V145.
#     v143 attaches the flank fasta by the slot in the file name, which gives 50.1% of haplotypes
#       the sequence of a different haplotype. Against the original FASTA over the 457 overlapping haplotypes,
#       the number matching at every bin is 263 of 457 for v143 and 457 of 457 for V145.
#     V145 holds 487 haplotypes where v143 holds 470.
A = A[[t in HP for t in zip(A["sample"].astype(str), A.slot.astype(str))]].copy()
o = np.array([off.get((s, sl, sd), np.nan) for s, sl, sd in
              zip(A["sample"].astype(str), A.slot.astype(str), A.side)])
b = A.b.values.astype(float)
A["u"] = np.where(b <= 0, o - np.abs(b), b + o)      # the same expression as fig1_v5_prep
A = A.dropna(subset=["u"])
print("[recoordinate] sequence rows %s, haplotypes %d" % (format(len(A), ","), A.groupby(["sample", "slot"]).ngroups))

# ---- re-bin onto a single window width of 1 kb ----
a = np.abs(A.u.values)
lo = (a // W) * W
A["u1k"] = np.where(A.u.values > 0, lo, -lo)
A = A[a <= 52000]
g = (A.groupby(["side", "u1k"])
       .agg(bp=("bp", "sum"), n_cpg=("n_cpg", "sum"), n_alu=("n_alu", "sum"),
            n_hap=("sample", "size")).reset_index())
g["cpg_per_kb"] = 1000.0 * g.n_cpg / g.bp
g["alu_frac"] = g.n_alu / g.bp
g = g.rename(columns={"u1k": "u"}).sort_values(["side", "u"])
OUT = os.path.join(HERE, "FIG3_SUPP_A_SEQ1K.tsv")
g.to_csv(OUT, sep="\t", index=False, float_format="%.5f")

# ---- self-check ----
old = pd.read_csv(os.path.join(HERE, "FIG1_V5_A_bp.tsv"), sep="\t")
print()
# FIG1_V5_A_bp.tsv, the table this is compared against, is built from v143. The window is the same 1 kb, so
#   the difference is exactly the difference between the input tables.
#   It should not be 0; a difference of 0 would mean V145 had not taken effect.
print("check (1): how much the deep side (1 kb windows) differs from the v143 table")
n_same = n_diff = 0
for side in ("근위", "원위"):
    o2 = old[(old.side == side) & (old.u < -8000) & old.cpg_per_kb.notna()]
    n2 = g[(g.side == side)].set_index("u")
    for r in o2.itertuples(index=False):
        if r.u in n2.index:
            d = abs(float(n2.loc[r.u, "cpg_per_kb"]) - r.cpg_per_kb)
            if d < 0.05:
                n_same += 1
            else:
                n_diff += 1
                if n_diff <= 3:
                    print("   differs %s u %.0f: original %.2f vs new %.2f" % (side, r.u, r.cpg_per_kb, n2.loc[r.u, "cpg_per_kb"]))
print("   bins equal to v143 %d, bins that differ %d; the windows are the same, so this difference is the amount of the change"
      % (n_same, n_diff))
if n_diff == 0:
    sys.exit("stopping: no bin differs. V145 was not read.")
print()
print("check (2): fluctuation per interval (median difference between neighbouring bins)")
for side in ("근위", "원위"):
    for lo_, hi_, nm in ((0, 8000, "0~8 kb (전에는 50·200 bp 창)"), (8000, 43000, "8~43 kb (전에도 1 kb 창)")):
        for T, tag in ((old, "옛"), (g, "새")):
            d = T[(T.side == side) & (T.u < 0) & T.cpg_per_kb.notna()].copy()
            d["dist"] = -d.u
            q = d[(d.dist >= lo_) & (d.dist < hi_)].sort_values("dist")
            if len(q) > 2:
                j = np.median(np.abs(np.diff(q.cpg_per_kb.values)))
                print("   %-4s %-28s %s sawtooth %6.2f, median %6.2f, bins %d"
                      % (side, nm, tag, j, q.cpg_per_kb.median(), len(q)))
print()
print("table written %s, bins %d (proximal %d, distal %d)"
      % (os.path.basename(OUT), len(g), (g.side == "근위").sum(), (g.side == "원위").sum()))
