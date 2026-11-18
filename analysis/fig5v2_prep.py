# -*- coding: utf-8 -*-
# Fig 5 (v2) data: the zones are recounted as the canonical eight, with proximal on the left.
#   1) First reproduce the original table (TRIO_SRCDEF_PAIRZONE) from the merged READMAP, to confirm the orientation and the counting.
#   2) Then build the eight zones with that same counting.
#   The canonical zones: Z1 proximal 1-3, Z2 proximal 4-10, Z3 proximal 11-17, ZI interior, Z4 distal 11-22, Z5 distal 9-10, Z6 distal 3-8, Z7 distal 1-2.
import io, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats
# The data roots come from environment variables.
#   They have no default: the script stops with a message if one is not set.
DATA_ROOT = os.environ["DATA_ROOT"]
T = DATA_ROOT + "/02_tables/"
# Where the tables are read and written. The tables are not deposited next to the code, so these paths come from the environment.
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
LOG = io.open(os.path.join(HERE, "FIG5V2_PREP.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()

PZ0 = pd.read_csv(T + "TRIO_SRCDEF_PAIRZONE.tsv", sep="\t")
pairs = PZ0[["child", "child_slot", "parent", "role", "trans_slot", "nontrans_slot", "call"]].drop_duplicates()
P("[pairs] %d (children %d, families %d)" % (len(pairs), pairs.child.nunique(), pairs.parent.nunique()))

ids = set(pairs.child) | set(pairs.parent)
R = pd.read_csv(DATA_ROOT + "/02_tables/READMAP_METH_ALL_v2_2026-09-27.tsv", sep="\t",
                usecols=["sample", "data_source", "slot", "array_pos", "n_reads", "beta_call"])
R = R[R["sample"].isin(ids)].copy()
P("[per-position β] from the merged set: trio %d individuals, rows %s, data source %s"
  % (R["sample"].nunique(), format(len(R), ","), sorted(R.data_source.unique())))
R["hap"] = R["sample"].astype(str) + "|" + R.slot.astype(str)
NM = R.groupby("hap").array_pos.max().to_dict()
R["N"] = R.hap.map(NM)

# ---------- 1) orientation check: the distal end must be the more depleted one ----------
a1 = R[R.array_pos == 1].beta_call.median()
aN = R[R.array_pos == R.N].beta_call.median()
P()
P("[direction] β median at array_pos 1 %.3f, β median at array_pos N %.3f -> %s"
  % (a1, aN, "array_pos 1 = distal (same as the canonical rule)" if a1 < aN else "반대 — 아래 셈법을 뒤집어야 한다"))
FLIP = a1 >= aN
R["ad"] = R.N - R.array_pos + 1 if FLIP else R.array_pos      # index counted from the distal end
R["ap"] = R.array_pos if FLIP else R.N - R.array_pos + 1      # index counted from the proximal end

# ---------- reproducing the original table (4 zones) ----------
def zone4(ap, ad):
    if ap <= 10: return "근위1-10"
    if ap <= 20: return "근위11-20"
    if ad <= 10: return "원위1-10"
    if ad <= 20: return "원위11-20"
    return None
def zone8(ap, ad):
    if ap <= ad:
        if ap <= 3: return "Z1 근위1-3"
        if ap <= 10: return "Z2 근위4-10"
        if ap <= 17: return "Z3 근위11-17"
        return "ZI 내부"
    if ad <= 2: return "Z7 원위1-2"
    if ad <= 8: return "Z6 원위3-8"
    if ad <= 10: return "Z5 원위9-10"
    if ad <= 22: return "Z4 원위11-22"
    return "ZI 내부"

R["z4"] = [zone4(a, d) for a, d in zip(R.ap, R.ad)]
def zone5(ap, ad):
    """Five common zones = the eight canonical zones merged
       Z1 / Z2+Z3 / ZI / Z4 / Z5+Z6+Z7. No new boundary is introduced.
       The condition must be the one used by main (E)(F)(G) and FigX26."""
    if ap <= ad:
        if ap <= 3:  return "근위 1-3"
        if ap <= 17: return "근위 4-17"
        return "내부"
    if ad <= 10: return "원위 1-10"
    if ad <= 22: return "원위 11-22"
    return "내부"

R["z8"] = [zone8(a, d) for a, d in zip(R.ap, R.ad)]
R["z5"] = [zone5(a, d) for a, d in zip(R.ap, R.ad)]
MINREAD = 2
Rq = R[R.n_reads >= MINREAD].copy()

def zone_beta(hap, zcol):
    v = Rq[Rq.hap == hap]
    return v.groupby(zcol).beta_call.median()

def pair_table(zcol, zorder):
    rows = []
    for _, p in pairs.iterrows():
        kid = "%s|%d" % (p.child, p.child_slot)
        tr = "%s|%s" % (p.parent, p.trans_slot)
        nt = None if str(p.nontrans_slot).strip() == "-" else "%s|%s" % (p.parent, p.nontrans_slot)
        bk, bt = zone_beta(kid, zcol), zone_beta(tr, zcol)
        bn = zone_beta(nt, zcol) if nt else None
        for z in zorder:
            if z not in bk.index or z not in bt.index: continue
            r = dict(child=p.child, child_slot=p.child_slot, parent=p.parent, role=p.role,
                     call=p.call, zone=z, b_child=bk[z], b_trans=bt[z],
                     d_trans=abs(bk[z] - bt[z]))
            if bn is not None and z in bn.index:
                r["b_non"] = bn[z]; r["d_non"] = abs(bk[z] - bn[z])
            rows.append(r)
    return pd.DataFrame(rows)

Z4ORD = ["근위1-10", "근위11-20", "원위11-20", "원위1-10"]
P4 = pair_table("z4", Z4ORD)
P()
P("[(1) reproduction] source table (4 zones) vs recomputation from the merged set: median |Δβ|")
P("   %-12s %10s %10s %10s %10s" % ("구역", "원표 전달", "reproduced transmitted", "원표 비전달", "reproduced non-transmitted"))
for z in Z4ORD:
    o = PZ0[PZ0.zone == z]; n = P4[P4.zone == z]
    P("   %-12s %10.4f %10.4f %10.4f %10.4f"
      % (z, o.d_trans.median(), n.d_trans.median(),
         o.d_non.median(), n.d_non.median() if "d_non" in n else np.nan))
P("   -> close values mean the same computation (median comparison, pairs %d vs %d)"
  % (PZ0[["child", "child_slot", "parent"]].drop_duplicates().shape[0], len(pairs)))

# ---------- 2) the eight zones ----------
Z8ORD = ["Z1 근위1-3", "Z2 근위4-10", "Z3 근위11-17", "ZI 내부",
         "Z4 원위11-22", "Z5 원위9-10", "Z6 원위3-8", "Z7 원위1-2"]
P8 = pair_table("z8", Z8ORD)
P8.to_csv(os.path.join(HERE, "FIG5V2_PAIRZONE8.tsv"), sep="\t", index=False, float_format="%.5f")
rows = []
P()
P("[(2) 8 zones] child to transmitted parental haplotype  vs  child to non-transmitted haplotype of the same parent")
P("   %-14s %5s %10s %10s %8s %9s" % ("구역", "짝", "전달 |Δβ|", "non-transmitted |Δβ|", "closer", "Wilcoxon p"))
for z in Z8ORD:
    g = P8[(P8.zone == z)].dropna(subset=["d_trans", "d_non"])
    if len(g) < 3:
        P("   %-14s %5d  (too few pairs, not tested)" % (z, len(g)))
        rows.append(dict(zone=z, n_pair=len(g), med_trans=np.nan, med_non=np.nan,
                         closer=np.nan, p=np.nan)); continue
    w = stats.wilcoxon(g.d_trans, g.d_non)
    rows.append(dict(zone=z, n_pair=len(g), med_trans=g.d_trans.median(), med_non=g.d_non.median(),
                     closer=int((g.d_trans < g.d_non).sum()), p=w.pvalue))
    P("   %-14s %5d %10.4f %10.4f %5d/%-3d %9.4f"
      % (z, len(g), g.d_trans.median(), g.d_non.median(),
         int((g.d_trans < g.d_non).sum()), len(g), w.pvalue))
SM = pd.DataFrame(rows)
SM.to_csv(os.path.join(HERE, "FIG5V2_SUMMARY8.tsv"), sep="\t", index=False, float_format="%.6g")

# ---------- 2b) the common five zones, as used by Supplementary J ----------
#   These group the canonical eight: Z1 / Z2+Z3 / ZI / Z4 / Z5+Z6+Z7. No new boundary is introduced.
Z5ORD = ["근위 1-3", "근위 4-17", "내부", "원위 11-22", "원위 1-10"]
P5 = pair_table("z5", Z5ORD)
P5.to_csv(os.path.join(HERE, "FIG5V2_PAIRZONE5.tsv"), sep="\t", index=False, float_format="%.5f")
rows5 = []
P()
P("[(2)-b 5 zones] same computation, zones reduced to the common 5")
P("   %-12s %5s %10s %10s %8s %9s" % ("구역", "짝", "전달 |Δβ|", "non-transmitted |Δβ|", "closer", "Wilcoxon p"))
for z in Z5ORD:
    g = P5[(P5.zone == z)].dropna(subset=["d_trans", "d_non"])
    if len(g) < 3:
        rows5.append(dict(zone=z, n_pair=len(g), med_trans=np.nan, med_non=np.nan,
                          n_closer=np.nan, p=np.nan))
        P("   %-12s %5d  (too few pairs, not tested)" % (z, len(g))); continue
    w = stats.wilcoxon(g.d_trans, g.d_non)
    rows5.append(dict(zone=z, n_pair=len(g), med_trans=g.d_trans.median(),
                      med_non=g.d_non.median(),
                      n_closer=int((g.d_trans < g.d_non).sum()), p=w.pvalue))
    P("   %-12s %5d %10.4f %10.4f %5d/%-3d %9.4f"
      % (z, len(g), g.d_trans.median(), g.d_non.median(),
         int((g.d_trans < g.d_non).sum()), len(g), w.pvalue))
SM5 = pd.DataFrame(rows5)
SM5.to_csv(os.path.join(HERE, "FIG5V2_ZONE5.tsv"), sep="\t", index=False, float_format="%.6g")
P("   [saved] FIG5V2_PAIRZONE5.tsv, FIG5V2_ZONE5.tsv")

# ---------- 3) can transmission in the flank be measured ----------
P()
P("[(3) flank transmission] whether 0-2 kb / 2-50 kb is transmitted")
FH = pd.read_csv(T + "FLANK_METH_HAP_A0_CLEAN.tsv", sep="\t")
have = sorted(set(FH["sample"]) & ids)
miss = sorted(ids - set(FH["sample"]))
P("   trio individuals present in the flank β table (FLANK_METH_HAP_A0_CLEAN): %s" % have)
P("   individuals not present: %s" % miss)
P("   -> there is no flank β for the parents (parents are 1kGP-ONT and GIAB; flank windows are computed only for HPRC assemblies)")
P("      so 'transmitted haplotype vs non-transmitted haplotype' cannot be compared in the flank with the present data.")
P("      it requires running the flank window pipeline on the parental ONT data (node job).")
pd.DataFrame(dict(있음=[", ".join(have)], 없음=[", ".join(miss)])).to_csv(
    os.path.join(HERE, "FIG5V2_FLANK_AVAIL.tsv"), sep="\t", index=False)
P()
P("[saved] FIG5V2_PAIRZONE8.tsv, FIG5V2_SUMMARY8.tsv, FIG5V2_FLANK_AVAIL.tsv")
LOG.close()
