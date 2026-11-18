# -*- coding: utf-8 -*-
# Fig 1-2 data: the orange line is recomputed on the same ruler as the black one, sum(methylated) / sum(calls).
#   The two lines are put on the same statistic. The values go into a new table and the existing table is left alone.
#   Gate: the new columns are added only after the existing table (median and unweighted mean) has been reproduced exactly.
#   The same filter as fig1_v4_prep: gate == 1 and HPRC_ONT_submissions.
#   The weight is n_canon, the number of CpG calls in that read-unit.
#   This is not the same quantity as the black line's sum_p / n_call.
#     The weight is the same but the quantity averaged differs: on the unit axis here the hard call b_call is averaged,
#     while the black line (fig1_v5_prep.py, FIG1_V5_A_bp.tsv) averages the probability sum sum_p.
#     They can be told apart: b_call * n_canon is an integer, being a hard call, while sum_p is not.
#     Size of the difference, over the same read x unit: at distal unit 1 the hard call gives 0.028 against 0.149 for the probability mean (+0.121),
#              and at proximal unit 10, 0.724 against 0.720 (-0.004). The gap opens only where beta is near 0.
#     Overlaying the two lines therefore leaves a systematic offset at the distal end; the measured median difference is +0.103.
#       The legend of Figure 3A records this.
import io, os, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd

HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))   # 
S = os.environ["DATA_ROOT"].rstrip("/") + "/"   # required; no default
md5 = lambda p: hashlib.md5(io.open(p, "rb").read()).hexdigest()
NED, MIDB = 20, 12

D = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "ap", "ad", "N", "b_call", "gate",
                             "data_source", "n_canon"])
D = D[(D.gate == 1) & (D.data_source.astype(str) == "HPRC_ONT_submissions")].copy()
print("read-units %d, haplotypes %d" % (len(D), D.groupby(["sample", "slot"]).ngroups))

# ---------- the 12 interior bins ----------
d = D[D.N >= 2 * NED + MIDB].copy()
mid = d[(d.ap > NED) & (d.ad > NED)].copy()
u = (mid.ap - NED - 1) / np.maximum(mid.N - 2 * NED - 1, 1)
mid["k"] = np.clip((u * MIDB).astype(int), 0, MIDB - 1)
gm = mid.groupby("k").apply(lambda v: pd.Series({
    "med": v.b_call.median(), "mean_simple": v.b_call.mean(),
    "mean_pooled": (v.b_call * v.n_canon).sum() / v.n_canon.sum(),
    "n": len(v), "n_call": int(v.n_canon.sum())}), include_groups=False).reset_index()

OLDM = pd.read_csv(HERE + "/FIG1_V4_A_mid.tsv", sep="\t")
m = OLDM.merge(gm[["k", "med", "mean_simple"]], on="k", suffixes=("_old", "_new"))
d1 = (m.med_old - m.med_new.round(4)).abs().max()
d2 = (m["mean"] - m.mean_simple.round(4)).abs().max()
print()
print("gate (1) interior 12 bins: max median difference %.6f, max simple-mean difference %.6f" % (d1, d2))
if max(d1, d2) > 1e-9:
    sys.exit("reproduction failed: stopping")
print("   -> reproduced exactly")
gm.to_csv(HERE + "/FIG1_2_MID.tsv", sep="\t", index=False, float_format="%.5f")

# ---------- the 20 units at each end ----------
rows = []
for side, col in (("근위", "ap"), ("원위", "ad")):
    g = D[D[col] <= NED]
    for i, v in g.groupby(col):
        w = v.n_canon.astype(float)
        rows.append(dict(side=side, idx=int(i), med=float(v.b_call.median()),
                         mean_simple=float(v.b_call.mean()),
                         mean_pooled=float((v.b_call * w).sum() / w.sum()),
                         n=int(len(v)), n_call=int(w.sum())))
gu = pd.DataFrame(rows)
OLDU = pd.read_csv(HERE + "/FIG1_V4_UNITS.tsv", sep="\t")
mu = OLDU.merge(gu[["side", "idx", "med"]], on=["side", "idx"], suffixes=("_old", "_new"))
d3 = (mu.med_old - mu.med_new.round(4)).abs().max()
print("gate (1) end 20 units: max median difference %.6f" % d3)
if d3 > 1e-9:
    sys.exit("reproduction failed: stopping")
print("   -> reproduced exactly")
gu.to_csv(HERE + "/FIG1_2_UNITS.tsv", sep="\t", index=False, float_format="%.5f")
print()
print("saved: FIG1_2_UNITS.tsv (md5 %s), FIG1_2_MID.tsv (md5 %s)"
      % (md5(HERE + "/FIG1_2_UNITS.tsv")[:12], md5(HERE + "/FIG1_2_MID.tsv")[:12]))
print()
print("interior 12 bins: median %s" % " ".join("%.2f" % v for v in gm.med))
print("             weighted mean %s" % " ".join("%.2f" % v for v in gm.mean_pooled))
