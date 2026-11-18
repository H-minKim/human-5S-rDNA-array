# -*- coding: utf-8 -*-
# Fig 1 (v4) panel A data: HPRC ONT only, with the flank and the array on one continuous bp axis.
#  Inputs: 02_tables/JUNCTION_FINE_METH.tsv (haplotype x side x bp window) and 02_tables/L1_READ_UNIT.parquet for the interior, on the unit axis.
#         02_tables/JUNCTION_SEQ_ANNOT.tsv, the CpG density and Alu-like fraction over the same windows.
#  Outputs: FIG1_V4_A_bp.tsv (side x window x methylation and sequence), FIG1_V4_A_mid.tsv (the 12 length-normalised interior bins) and FIG1_V4_UNITS.tsv.
import io, os, hashlib
import numpy as np, pandas as pd
S = os.environ["DATA_ROOT"].rstrip("/") + "/"   # required; no default
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))   # 
LOG = io.open(os.path.join(HERE, "FIG1_V4_PREP.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
INP = [S + "02_tables/JUNCTION_FINE_METH.tsv", S + "02_tables/L1_READ_UNIT.parquet", S + "02_tables/JUNCTION_SEQ_ANNOT.tsv"]
pd.DataFrame([(p, md5(p), os.path.getsize(p)) for p in INP], columns=["path", "md5", "bytes"]).to_csv(
    os.path.join(HERE, "FIG1_V4_INPUTS.tsv"), sep="\t", index=False)
P("[input]")
for p in INP: P("   %-30s %s %12d" % (os.path.basename(p), md5(p), os.path.getsize(p)))

D = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "read_id", "data_source", "ap", "ad", "N", "b_call", "gate"])
D = D[(D.gate == 1) & (D.data_source.astype(str) == "HPRC_ONT_submissions")].copy()
D["hap"] = D["sample"] + "|" + D.slot.astype(str)
HP = set(zip(D["sample"].astype(str), D.slot.astype(str)))
P()
P("[HPRC] haplotypes %d, units %d" % (len(HP), len(D)))

J = pd.read_csv(S + "02_tables/JUNCTION_FINE_METH.tsv", sep="\t", dtype={"slot": str})
J = J[J["sample"].astype(str) != "HG008-T"].copy()   # The flank also uses HPRC data only; HG008-T is excluded as a cancer sample with mixed PacBio data.
J["key"] = list(zip(J["sample"].astype(str), J.slot.astype(str)))
J = J[J.key.isin(HP)]
P("[bp windows] HPRC haplotypes %d, rows %d" % (J.groupby(["sample", "slot"]).ngroups, len(J)))
J["beta"] = J.sum_p / J.n_call
agg = J.groupby(["side", "b"]).apply(
    lambda x: pd.Series({"beta": x.sum_p.sum() / x.n_call.sum(),
                         "q1": x.beta.quantile(.25), "q3": x.beta.quantile(.75),
                         "n_call": x.n_call.sum(), "n_hap": x.groupby(["sample", "slot"]).ngroups})).reset_index()
try:
    A = pd.read_csv(S + "02_tables/JUNCTION_SEQ_ANNOT.tsv", sep="\t", dtype={"slot": str})
    A["key"] = list(zip(A["sample"].astype(str), A.slot.astype(str)))
    A = A[A.key.isin(HP)]
    sq = A.groupby(["side", "b"]).agg(bp=("bp", "sum"), cg=("n_cpg", "sum"), al=("n_alu", "sum")).reset_index()
    sq["cpg_per_kb"] = 1000 * sq.cg / sq.bp
    sq["alu_frac"] = sq.al / sq.bp
    agg = agg.merge(sq[["side", "b", "cpg_per_kb", "alu_frac", "bp"]], on=["side", "b"], how="left")
    P("[sequence] annotation attached, haplotypes %d" % A.groupby(["sample", "slot"]).ngroups)
except FileNotFoundError:
    P("[sequence] none yet")
agg.to_csv(os.path.join(HERE, "FIG1_V4_A_bp.tsv"), sep="\t", index=False, float_format="%.5f")
P("   windows %d (windows with 50 or more haplotypes %d)" % (len(agg), (agg.n_hap >= 50).sum()))

# the interior, on the unit axis: 12 length-normalised bins, outside the 20 units at each end
NED, MIDB = 20, 12
d = D[D.N >= 2 * NED + MIDB].copy()
m = (d.ap > NED) & (d.ad > NED)
mid = d[m].copy()
u = (mid.ap - NED - 1) / np.maximum(mid.N - 2 * NED - 1, 1)
mid["k"] = np.clip((u * MIDB).astype(int), 0, MIDB - 1)
g = mid.groupby("k").b_call.agg(["median", "mean", "size", lambda v: v.quantile(.25), lambda v: v.quantile(.75)])
g.columns = ["med", "mean", "n", "q1", "q3"]
g.reset_index().to_csv(os.path.join(HERE, "FIG1_V4_A_mid.tsv"), sep="\t", index=False, float_format="%.4f")
P("[interior] haplotypes with N>=%d: inner units %d, 12-bin median β %s"
  % (2 * NED + MIDB, len(mid), " ".join("%.2f" % v for v in g["med"])))

# unit-axis values for the 20 units at each end, used to place unit ticks on the bp axis
un = []
for side, col in (("근위", "ap"), ("원위", "ad")):
    v = D[D[col] <= 20].groupby(col).b_call.agg(["median", "size"])
    for i, r in v.iterrows(): un.append((side, int(i), r["median"], int(r["size"])))
pd.DataFrame(un, columns=["side", "idx", "med", "n"]).to_csv(
    os.path.join(HERE, "FIG1_V4_UNITS.tsv"), sep="\t", index=False, float_format="%.4f")
P("[unit] median of the end 20 units saved")
LOG.close()
