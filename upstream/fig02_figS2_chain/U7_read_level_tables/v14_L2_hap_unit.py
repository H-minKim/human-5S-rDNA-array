# -*- coding: utf-8 -*-
# L2 = haplotype x position combined table.
#   The L1 values of the reads covering the same (sample, slot, array_pos) are summarised without losing bimodality.
#   With only mean and variance, "40% low methylation + 60% high methylation" and "all intermediate" cannot be told apart.
#   -> f_hi (mixing fraction) and per-state values (b_lo_mean / b_hi_mean) are reported separately.
#
#   (C) haplotype constants and (D) the opposite slot of the same individual (p_ prefix) are also attached.
import io, os
import numpy as np
import pandas as pd

OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"
T = OUT + "/02_tables"
LOGD = OUT + "/04_logs"
THR = 0.508
LO = 0.35
HI = 0.65

LOG = io.open(os.path.join(LOGD, "V14_L2.txt"), "w", encoding="utf8", newline="\n")


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


P("=" * 100)
P("L2 = haplotype x position integrated table (bimodality-preserving summary)")
P("=" * 100)
P()

L1 = pd.read_parquet(os.path.join(T, "L1_READ_UNIT.parquet"))
P("  L1 %s rows" % format(len(L1), ","))
L1 = L1[L1.gate == 1].copy()
P("  only read-units passing the gate are used : %s rows (%.4f)"
  % (format(len(L1), ","), len(L1) / 787885))
P()

L1["lo"] = (L1.b_call < LO).astype("float32")      # hard call
L1["mid"] = ((L1.b_call >= LO) & (L1.b_call <= HI)).astype("float32")   # hard call
L1["hi"] = (L1.b_call > HI).astype("float32")      # hard call
L1["noise_split"] = ((L1.bA - L1.bB) ** 2 / 4.0).astype("float32")
L1["noise_binom"] = (L1.b_call * (1 - L1.b_call) / np.maximum(L1.n_canon, 1)).astype("float32")  # hard call
L1["switch_next"] = np.where(L1.st_next.notna(), (L1.st != L1.st_next).astype(float), np.nan)

# split halves: read_id hash parity (reproducible)
L1["half"] = (pd.util.hash_pandas_object(L1.read_id, index=False) % 2).astype("int8")

K = ["sample", "slot", "array_pos"]
G = L1.groupby(K, observed=True, sort=False)

P("  positions %s" % format(G.ngroups, ","))
P()

agg = G.agg(
    n_read=("read_id", "size"),
    f_hi_thr=("st", "mean"),
    f_lo=("lo", "mean"),
    f_mid=("mid", "mean"),
    f_hi=("hi", "mean"),
    b_mean=("b_call", "mean"),          # hard call
    b_sd=("b_call", "std"),             # hard call
    b_call_mean=("b_call", "mean"),
    noise_split=("noise_split", "mean"),
    noise_binom=("noise_binom", "mean"),
    b_alu_mean=("b_alu", "mean"),
    b_nts1_mean=("b_nts1", "mean"),
    b_ga_mean=("b_ga", "mean"),
    b_nts2_mean=("b_nts2", "mean"),
    b_ups_mean=("b_ups", "mean"),
    b_gene_mean=("b_gene", "mean"),
    b_sp_mean=("b_sp", "mean"),
    agree_prev=("agree_prev", "mean"),
    agree_next=("agree_next", "mean"),
    p_switch=("switch_next", "mean"),
    blk_len_mean=("blk_len", "mean"),
    n_canon_mean=("n_canon", "mean"),
    frac_jx=("kind", lambda s: float((s.astype(str) == "junction").mean())),
).reset_index()
agg["b_var"] = agg.b_sd ** 2
agg["b_var_true"] = np.maximum(agg.b_var - agg.noise_split, 0.0)

# per-state values (the core of preserving bimodality)
LOW = L1[L1.st == 0].groupby(K, observed=True, sort=False).agg(
    n_lo=("b_call", "size"), b_lo_mean=("b_call", "mean"), b_lo_sd=("b_call", "std"),   # hard call
    b_lo_gene=("b_gene", "mean"), b_lo_sp=("b_sp", "mean"),
    b_lo_ups=("b_ups", "mean"), b_lo_alu=("b_alu", "mean")).reset_index()
HIG = L1[L1.st == 1].groupby(K, observed=True, sort=False).agg(
    n_hi=("b_call", "size"), b_hi_mean=("b_call", "mean"), b_hi_sd=("b_call", "std"),   # hard call
    b_hi_gene=("b_gene", "mean"), b_hi_sp=("b_sp", "mean"),
    b_hi_ups=("b_ups", "mean"), b_hi_alu=("b_alu", "mean")).reset_index()
agg = agg.merge(LOW, on=K, how="left").merge(HIG, on=K, how="left")
agg["r_lo_gene_sp"] = agg.b_lo_gene / agg.b_lo_sp.replace(0, np.nan)
agg["r_hi_gene_sp"] = agg.b_hi_gene / agg.b_hi_sp.replace(0, np.nan)

# split-half reproducibility
HA = L1[L1.half == 0].groupby(K, observed=True, sort=False).agg(
    n_A=("st", "size"), f_A=("st", "mean")).reset_index()
HB = L1[L1.half == 1].groupby(K, observed=True, sort=False).agg(
    n_B=("st", "size"), f_B=("st", "mean")).reset_index()
agg = agg.merge(HA, on=K, how="left").merge(HB, on=K, how="left")

# fixed values per position and unit (independent of reads)
FIX = L1.drop_duplicates(K)[
    K + ["ap", "ad", "side", "ord_near", "pos_label", "rel", "N",
         "unit_len", "n_var", "n_var_gene", "varlist",
         "n_var_common", "n_var_rare", "n_var_unmatched",
         "n_snv", "n_ins", "n_del", "n_var_cpg", "n_var_coding",
         "nv_AluY", "nv_NTS1", "nv_GArich", "nv_NTS2", "nv_coding", "nv_NA",
         "vsig", "vsig_n", "vsig_minfreq", "chain_ok", "gap_to_next",
         "n_read_slot", "gate_rate", "cov_frac", "depth", "n_break", "revcomp",
         "data_source", "pop", "super_pop"]]
L2 = FIX.merge(agg, on=K, how="left")

# (C) haplotype constants: boundary and number of low-methylation units
HB2 = pd.read_csv(os.path.join(T, "HAP_LENGTH_BOUNDARY.tsv"), sep="\t",
                  dtype={"sample": str, "slot": str})
HB2 = HB2.rename(columns={"f": "f_hap"})[["sample", "slot", "f_hap", "Bp", "Bd", "n_low"]]
L2 = L2.merge(HB2, on=["sample", "slot"], how="left")
P("  fraction with haplotype constants (f_hap, Bp, Bd, n_low) attached %.4f" % L2.f_hap.notna().mean())
P()

# ---------------------------------------------------------------- (D) opposite slot
PCOLS = ["n_read", "f_hi", "f_lo", "f_mid", "b_mean", "b_var", "b_var_true",
         "b_lo_mean", "b_hi_mean", "b_gene_mean", "b_sp_mean",
         "r_lo_gene_sp", "r_hi_gene_sp", "agree_next", "p_switch",
         "n_var_rare", "vsig_n", "N", "f_hap", "Bp", "Bd", "n_low"]
OTH = {"1": "2", "2": "1"}
Pt = L2[["sample", "slot", "ap", "ad", "rel"] + PCOLS].copy()
Pt["slot"] = Pt.slot.astype(str).map(OTH)

# (D-1) matched by the same proximal order (ap)
p_ap = Pt.drop(columns=["ad", "rel"]).rename(
    columns={c: "p_" + c for c in PCOLS})
L2 = L2.merge(p_ap, on=["sample", "slot", "ap"], how="left")
# (D-2) also by the same distal order (ad), f_hi only
p_ad = Pt[["sample", "slot", "ad", "f_hi", "b_mean"]].rename(
    columns={"f_hi": "pd_f_hi", "b_mean": "pd_b_mean"})
L2 = L2.merge(p_ad, on=["sample", "slot", "ad"], how="left")
# (D-3) individual level: N_other, dN
NT = L2.drop_duplicates(["sample", "slot"])[["sample", "slot", "N"]].copy()
NT["slot"] = NT.slot.astype(str).map(OTH)
NT = NT.rename(columns={"N": "N_other"})
L2 = L2.merge(NT, on=["sample", "slot"], how="left")
L2["dN"] = (L2.N - L2.N_other).abs()
L2["both_eligible"] = L2.N_other.notna().astype("int8")

P("  positions with the opposite slot present %.4f (individuals %d)"
  % (L2.N_other.notna().mean(),
     L2.loc[L2.N_other.notna(), "sample"].nunique()))
P("  of these, the partner position exists at the same ap %.4f" % L2.p_f_hi.notna().mean())
P()

L2.to_parquet(os.path.join(T, "L2_HAP_UNIT.parquet"), index=False)
P("  saved : 02_tables/L2_HAP_UNIT.parquet  %s rows x %d columns"
  % (format(len(L2), ","), L2.shape[1]))
P()

# ---------------------------------------------------------------- L3 haplotype
L3 = L2.groupby(["sample", "slot"], observed=True).agg(
    N=("N", "first"), n_pos=("array_pos", "size"),
    n_read_slot=("n_read_slot", "first"),
    depth=("depth", "first"), cov_frac=("cov_frac", "first"),
    gate_rate=("gate_rate", "first"), n_break=("n_break", "first"),
    f_hap=("f_hap", "first"), Bp=("Bp", "first"), Bd=("Bd", "first"),
    n_low=("n_low", "first"),
    f_hi_mean=("f_hi", "mean"),
    b_mean=("b_mean", "mean"),
    n_var_rare_mean=("n_var_rare", "mean"),
    vsig_n_mean=("vsig_n", "mean"),
    N_other=("N_other", "first"), dN=("dN", "first"),
    pop=("pop", "first"), super_pop=("super_pop", "first"),
).reset_index()
for lab, m in (("prox8", L2.ap <= 8), ("prox9_14", (L2.ap >= 9) & (L2.ap <= 14)),
               ("mid", (L2.ap > 20) & (L2.ad > 20)),
               ("dist13_17", (L2.ad >= 13) & (L2.ad <= 17)),
               ("dist8", L2.ad <= 8)):
    s = L2[m].groupby(["sample", "slot"], observed=True)["f_hi"].mean().rename("f_hi_" + lab)
    L3 = L3.merge(s.reset_index(), on=["sample", "slot"], how="left")
L3.to_csv(os.path.join(T, "L3_HAP.tsv"), sep="\t", index=False, float_format="%.6g")
P("  saved : 02_tables/L3_HAP.tsv  %d slots x %d columns" % (len(L3), L3.shape[1]))
P()

# ---------------------------------------------------------------- check
P("=" * 100)
P("check: does this table reproduce the conclusions")
P("=" * 100)
P()
P("  [1] is bimodality preserved (per L2 position)")
P("      distribution of f_hi at positions with n_read>=10")
s = L2[L2.n_read >= 10]
q = s.f_hi.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).values
P("      5%% %.3f, 25%% %.3f, median %.3f, 75%% %.3f, 95%% %.3f" % tuple(q))
P("      f_hi below 0.05 (almost all low-methylation) %.3f, above 0.95 (almost all high-methylation) %.3f"
  % ((s.f_hi < 0.05).mean(), (s.f_hi > 0.95).mean()))
P("      = values pile up at both ends at the position level too (not in the middle)")
P()
P("  [2] is the state value independent of position (two-state model)")
for lab, m in (("근위 1-8", s.ap <= 8), ("가운데", (s.ap > 20) & (s.ad > 20)),
               ("원위 1-8", s.ad <= 8)):
    t = s[m]
    P("      %-8s  f_hi %.3f | low-methylation state value %.4f (n %s) | high-methylation state value %.4f (n %s)"
      % (lab, t.f_hi.mean(), t.b_lo_mean.mean(), format(int(t.n_lo.sum()), ","),
         t.b_hi_mean.mean(), format(int(t.n_hi.sum()), ",")))
P("      -> state values stay nearly the same across positions, only f_hi changes")
P()
P("  [3] true variance vs noise")
P("      b_var mean %.5f, split noise %.5f, binomial noise %.5f, true variance %.5f (%.1f%%)"
  % (s.b_var.mean(), s.noise_split.mean(), s.noise_binom.mean(),
     s.b_var_true.mean(), 100 * s.b_var_true.mean() / s.b_var.mean()))
P()
P("  [4] split-half reproducibility (f_A vs f_B)")
t = s[(s.n_A >= 5) & (s.n_B >= 5)]
r = np.corrcoef(t.f_A, t.f_B)[0, 1]
P("      r(f_A, f_B) = %.4f (n=%s), Spearman-Brown ceiling %.4f"
  % (r, format(len(t), ","), 2 * r / (1 + r)))
P()
P("  [5] both slots : is f_hi similar at the same position")
t = L2[(L2.n_read >= 10) & L2.p_f_hi.notna() & (L2.p_n_read >= 10)]
P("      r(f_hi, p_f_hi) = %.4f (positions %s, individuals %d)"
  % (np.corrcoef(t.f_hi, t.p_f_hi)[0, 1], format(len(t), ","), t["sample"].nunique()))
P("      note: the position effect (both at the same ap) is shared, so a positive value is expected; look at it after subtracting the per-position mean")
tt = t.copy()
mu = tt.groupby("ap")["f_hi"].transform("mean")
mp = tt.groupby("ap")["p_f_hi"].transform("mean")
P("      r after subtracting the position mean = %.4f"
  % np.corrcoef(tt.f_hi - mu, tt.p_f_hi - mp)[0, 1])
P()
LOG.close()
print("done")
