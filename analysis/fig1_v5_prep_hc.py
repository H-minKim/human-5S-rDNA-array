# -*- coding: utf-8 -*-
# Fig 1 (v5) panel A data: the start of the array is moved to the first complete unit.
#   The incomplete (partial) units at either end do not belong to the array.
#   Partial unit lengths from the G8 map: PARTIAL_PROX lies before array_pos 1, at the distal end (median 1,635 bp),
#                                  and PARTIAL_DIST lies after the last unit, at the proximal end (median 1,928 bp).
#     (In L1, array_pos equals ad, so array_pos 1 in the map frame is the distal end.)
#   The new coordinate u is the distance in bp from the boundary of the first complete unit: u > 0 lies outside, in the partial unit and then the flank, and u < 0 lies inside the array, in complete units.
#   Outputs: FIG1_V5_A_bp.tsv and FIG1_V5_CHECK.txt.
#
#   The hard-call version.
#     The window beta as a probability mean, sum_p / n_call, and the unit axis of the same figure,
#     a hard-call fraction, are on different rulers (+0.103 systematically at the distal end).
#     So sum_h, the count of calls at ML >= 128, is read from JUNCTION_FINE_METH_HC.tsv and
#       beta_hc reported alongside. The existing beta column, the probability mean, is kept.
#     Output: FIG1_V5_A_bp_HC.tsv; the existing file is not overwritten.
import io, os, hashlib
import numpy as np, pandas as pd
# The data roots come from environment variables.
#   They have no default: the script stops with a message if one is not set.
DATA_ROOT = os.environ["DATA_ROOT"]
V_ROOT = os.environ["V_ROOT"]
S = DATA_ROOT + "/"
SUMF = V_ROOT + "/03_maps/g8canonv10/CONSENSUS_SUM_G8_full.tsv"
# Where the tables are read and written. The tables are not deposited next to the code, so these paths come from the environment.
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
LOG = io.open(os.path.join(HERE, "FIG1_V5_PREP_HC.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
UNIT = 2240.0

D = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "data_source", "ap", "ad", "N", "b_call", "gate"])
D = D[(D.gate == 1) & (D.data_source.astype(str) == "HPRC_ONT_submissions")].copy()
HP = set(zip(D["sample"].astype(str), D.slot.astype(str)))
u1_prox = D[D.ap == 1].b_call.median(); u1_dist = D[D.ad == 1].b_call.median()
u2_prox = D[D.ap == 2].b_call.median(); u2_dist = D[D.ad == 2].b_call.median()
P("[unit-axis reference value (L1)] proximal first complete unit β %.3f (second %.3f), distal first complete unit β %.3f (second %.3f)"
  % (u1_prox, u2_prox, u1_dist, u2_dist))

SUM = pd.read_csv(SUMF, sep="\t", dtype={"slot": str})
SUM["off_원위"] = SUM.off_prox          # the partial unit at the start of the map is the distal one
SUM["off_근위"] = SUM.off_dist
P("[partial unit length] distal median %d bp (10-90%% %d-%d), proximal median %d bp (10-90%% %d-%d)"
  % (SUM.off_원위.median(), SUM.off_원위.quantile(.1), SUM.off_원위.quantile(.9),
     SUM.off_근위.median(), SUM.off_근위.quantile(.1), SUM.off_근위.quantile(.9)))

J = pd.read_csv(S + "02_tables/JUNCTION_FINE_METH_HC.tsv", sep="\t", dtype={"slot": str})
J = J[J["sample"].astype(str) != "HG008-T"].copy()   # The flank also uses HPRC data only (common rule 3); HG008-T is excluded as a cancer sample with mixed PacBio data.
J = J[[t in HP for t in zip(J["sample"].astype(str), J.slot.astype(str))]].copy()
A = pd.read_csv(S + "02_tables/JUNCTION_SEQ_ANNOT.tsv", sep="\t", dtype={"slot": str})
A = A[[t in HP for t in zip(A["sample"].astype(str), A.slot.astype(str))]].copy()
off = {}
for r in SUM.itertuples(index=False):
    off[(r.sample, str(r.slot), "원위")] = float(r.off_원위)
    off[(r.sample, str(r.slot), "근위")] = float(r.off_근위)
def newu(df):
    o = np.array([off.get((s, sl, sd), np.nan) for s, sl, sd in
                  zip(df["sample"].astype(str), df.slot.astype(str), df.side)])
    b = df.b.values.astype(float)
    # b <= 0 lies inside the array, with |b| the distance from the junction, so u = o - |b|; u > 0 is still within the partial unit.
    # b>0  : flank                             → u = b + o
    return np.where(b <= 0, o - np.abs(b), b + o)
for T in (J, A):
    T["u"] = newu(T)
J = J.dropna(subset=["u"]); A = A.dropna(subset=["u"])
P("[recoordinate] methylation rows %d, sequence rows %d, haplotypes %d" % (len(J), len(A), J.groupby(["sample", "slot"]).ngroups))

def rb(u):
    a = abs(u)
    if a < 2000: w = 50 if u < 0 else 100
    elif a < 10000: w = 200 if u < 0 else 500
    elif a <= 52000: w = 1000 if u < 0 else 5000
    else: return np.nan
    lo = (a // w) * w
    return lo if u > 0 else -lo
J["u2"] = [rb(x) for x in J.u]; A["u2"] = [rb(x) for x in A.u]
J = J.dropna(subset=["u2"]); A = A.dropna(subset=["u2"])
J["beta"] = J.sum_p / J.n_call
J["beta_hc"] = J.sum_h / J.n_call          # the hard-call fraction
g = J.groupby(["side", "u2"]).apply(lambda x: pd.Series(
    {"beta": x.sum_p.sum() / x.n_call.sum(),
     "beta_hc": x.sum_h.sum() / x.n_call.sum(),          # the same window, on the hard-call ruler
     "q1": x.beta.quantile(.25), "q3": x.beta.quantile(.75),
     "q1_hc": x.beta_hc.quantile(.25), "q3_hc": x.beta_hc.quantile(.75),
     "n_call": x.n_call.sum(), "n_hap": x.groupby(["sample", "slot"]).ngroups})).reset_index()
s = A.groupby(["side", "u2"]).agg(bp=("bp", "sum"), cg=("n_cpg", "sum"), al=("n_alu", "sum")).reset_index()
s["cpg_per_kb"] = 1000 * s.cg / s.bp; s["alu_frac"] = s.al / s.bp
g = g.merge(s[["side", "u2", "cpg_per_kb", "alu_frac"]], on=["side", "u2"], how="left").rename(columns={"u2": "u"})
g.to_csv(os.path.join(HERE, "FIG1_V5_A_bp_HC.tsv"), sep="\t", index=False, float_format="%.5f")
P("[saved] FIG1_V5_A_bp_HC.tsv, windows %d" % len(g))

# ---- check: in the new coordinate, does beta over the first complete unit match the unit-axis value? ----
P()
P("[check] β of the 'first complete unit' interval measured in the new coordinate u vs L1 unit-axis value")
def band(side, lo, hi):
    v = J[(J.side == side) & (-J.u2 >= lo) & (-J.u2 < hi)]      # u < 0 is inside the array
    return v.sum_p.sum() / v.n_call.sum(), v.n_call.sum()
for side, ref1, ref2 in (("근위", u1_prox, u2_prox), ("원위", u1_dist, u2_dist)):
    b1, n1 = band(side, 0, UNIT); b2, n2 = band(side, UNIT, 2 * UNIT)
    P("   %s: first unit interval β %.3f (L1 %.3f), second unit interval β %.3f (L1 %.3f)" % (side, b1, ref1, b2, ref2))
    # for reference, the same span in the junction coordinate b, which includes the partial unit
    v = J[(J.side == side) & (J.b <= 0) & (-J.b < UNIT)]
    P("      (note) β of the same interval in coordinates that count the partial unit as array %.3f: there the proximal first unit looks high-methylation"
      % (v.sum_p.sum() / v.n_call.sum()))
LOG.close()
