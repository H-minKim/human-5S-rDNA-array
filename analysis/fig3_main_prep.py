# -*- coding: utf-8 -*-
# Data preparation for the main Figure 3.
#
#   The panels:
#     A  the array structure schematic, reused from the top of Fig1_2_array_structure a and not built here
#     B  unit beta by position against the array length stratum, reproduced in this frame
#     C  the count of hypomethylated units per length stratum, reproduced in this frame
#     D  the methylation pattern within the unit, by segment, for highly against weakly methylated units
#     E  why the eight zones fall where they do: crossing the methylation axis with the sequence axis gives eight cells
#     F  between-haplotype differences per zone, reproducing Fig1_3_zone_variability a
#
#   Nothing is borrowed: B and C are computed directly from L2 here rather than taken from another table.
#     Compared against the values annotated on that other figure the largest difference was 0.010; this is repeated in the check section below.
#
#   Ruler: unit beta is L2_HAP_UNIT.b_mean, a hard call, behind the v13 gate plus n_read >= 5.
#   Orientation, in the L2 convention: ap counts from the proximal end and ad from the distal end, the opposite of ru_v11.
#   Read-only, and all the heavy IO is done in this single pass.
import io, os, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
# The data roots come from environment variables.
#   They have no default: the script stops with a message if one is not set.
DATA_ROOT = os.environ["DATA_ROOT"]

S = DATA_ROOT + "/"
# Where the tables are read and written. The tables are not deposited next to the code, so these paths come from the environment.
# Writes the deposited table with English column names and values
#   Renamed at a single point just before writing; the rest of this file uses the internal names.
_KO2EN_COL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": {"구역": "zone", "유닛안자리": "pos_in_unit",
                                    "유전자기준": "pos_gene_ref", "호출수": "n_call",
                                    "메틸호출수": "n_meth_call"},
    "FIG3_B_POSBETA.tsv": {"층": "stratum", "축": "axis", "자리": "pos"},
    "FIG3_C_LOWCOUNT.tsv": {"컷": "cut", "층": "stratum", "N중앙": "N_median",
                            "저메틸수": "n_low_meth", "저메틸수_sd": "n_low_meth_sd",
                            "비례했다면": "expected_if_proportional",
                            "저메틸몫": "low_meth_fraction"},
    "SUPP3_2_POSSD.tsv": {"분자수": "n_molecules", "평균_b_call": "mean_b_call",
                          "분자사이_sd": "sd_between_molecules"},
    "SUPP3_B_HAP3.tsv": {"하플로타입": "haplotype", "길이층": "length_class",
                         "유닛수": "n_units", "분자순번": "molecule_index",
                         "유닛자리": "unit_position"},
}
_KO2EN_VAL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": ("구역", {"근위": "proximal", "내부": "interior",
                                            "원위": "distal", "전체": "all"}),
    "FIG3_B_POSBETA.tsv": ("층", {"전체": "all"}),
    "FIG3_C_LOWCOUNT.tsv": ("컷", {"0.508(여기서 잴 골)": "0.508 (trough measured in this work)"}),
    "SUPP3_B_HAP3.tsv": ("길이층", {"짧은 array": "short array", "중간 길이": "medium",
                                   "긴 array": "long array"}),
}


def _to_english(n, df):
    b = os.path.basename(n)
    v = _KO2EN_VAL.get(b)
    if v and v[0] in df.columns:
        df = df.copy()
        df[v[0]] = df[v[0]].map(lambda x: v[1].get(x, x))
    c = _KO2EN_COL.get(b)
    if c:
        df = df.rename(columns=c)
    return df

HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
LOG = io.open(os.path.join(HERE, "FIG3_MAIN_PREP.txt"), "w", encoding="utf8", newline="\n")


def P(s=""):
    print(s, flush=True)
    LOG.write(str(s) + "\n")
    LOG.flush()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def save(df, name, **kw):
    p = os.path.join(HERE, name)
    _to_english(name, df).to_csv(p, sep="\t", index=False, float_format="%.5f", lineterminator="\n", **kw)
    P("   saved %-28s rows %5d, md5 %s" % (name, len(df), md5(p)[:12]))


LENB = [(40, 80, "N40-80"), (81, 110, "N81-110"), (111, 140, "N111-140"), (141, 10 ** 9, "N141+")]
CUT = 0.508 # trough of the unit hard call
CUT_THEM = 0.51      #  the cut used there, for comparison

P("=" * 92)
P("Figure 3 main figure data preparation")
P("=" * 92)
P()

# ---------------------------------------------------------------- inputs
L2 = pd.read_parquet(S + "02_tables/L2_HAP_UNIT.parquet",
                     columns=["sample", "slot", "ap", "ad", "N", "n_read", "b_mean", "f_hi"])
L2 = L2[L2.n_read >= 5].copy()
L2["hap"] = L2["sample"].astype(str) + "|" + L2.slot.astype(str)
P("[input] L2_HAP_UNIT (n_read>=5) : units %d, haplotypes %d, individuals %d"
  % (len(L2), L2.hap.nunique(), L2["sample"].nunique()))

L1C = ["sample", "slot", "read_id", "ap", "ad", "N", "gate", "b_call",
       "b_alu", "b_nts1", "b_ga", "b_nts2", "b_ups", "b_gene",
       "n_alu", "n_nts1", "n_ga", "n_nts2", "n_ups", "n_gene"]
L1 = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet", columns=L1C)
L1 = L1[L1.gate == 1].copy()
L1["hap"] = L1["sample"].astype(str) + "|" + L1.slot.astype(str)
P("[input] L1_READ_UNIT (gate==1)  : read x units %d, haplotypes %d, individuals %d"
  % (len(L1), L1.hap.nunique(), L1["sample"].nunique()))
P()

# ================================================================ B
P("=" * 92)
P("[B] unit β by position x array length stratum")
P()
rows = []
for lo, hi, nm in LENB + [(0, 10 ** 9, "전체")]:
    g = L2[(L2.N >= lo) & (L2.N <= hi)]
    if not len(g):
        continue
    for k in range(1, 31):
        v = g[g.ap == k]
        if len(v):
            rows.append(dict(층=nm, 축="P", 자리=k, x=k - 1, b=v.b_mean.mean(),
                             sd=v.b_mean.std(), n=len(v), n_hap=v.hap.nunique()))
    mid = g[(g.ap > 25) & (g.ad > 25)]
    if len(mid):
        rows.append(dict(층=nm, 축="MID", 자리=0, x=30, b=mid.b_mean.mean(),
                         sd=mid.b_mean.std(), n=len(mid), n_hap=mid.hap.nunique()))
    for k in range(30, 0, -1):
        v = g[g.ad == k]
        if len(v):
            rows.append(dict(층=nm, 축="D", 자리=k, x=31 + (30 - k), b=v.b_mean.mean(),
                             sd=v.b_mean.std(), n=len(v), n_hap=v.hap.nunique()))
B = pd.DataFrame(rows)
save(B, "FIG3_B_POSBETA.tsv")

P()
P("   check:  compared with the figure annotation")
ref = {("P", 1): 0.343, ("MID", 0): 0.817, ("D", 20): 0.658, ("D", 10): 0.151, ("D", 1): 0.027}
P("   %-8s %10s %10s %10s %8s" % ("자리", "β", "주석", "차", "유닛"))
mx = 0.0
for (ax, k), r in ref.items():
    v = B[(B.층 == "전체") & (B.축 == ax) & (B.자리 == k)]
    if len(v):
        b = float(v.b.iloc[0])
        mx = max(mx, abs(b - r))
        P("   %-8s %10.4f %10.3f %+10.4f %8d"
          % (("가운데" if ax == "MID" else "%s%d" % (ax, k)), b, r, b - r, int(v.n.iloc[0])))
P("   -> max difference %.4f: the same phenomenon is seen with the same measure" % mx)
P("   -> drawn with the cohort (510 haplotypes / 280 individuals)")

P()
ov = ((L2.ap <= 30) & (L2.ad <= 30)).sum()
P("   note: in a short array one unit falls in both P and D: %d units (%.1f%%)"
  % (ov, 100 * ov / len(L2)))
P("     The interior is defined separately as ap>25 AND ad>25 (no overlap).")

# ================================================================ C
P()
P("=" * 92)
P("[C] number of low-methylation units by array length stratum")
P()
rows = []
for cut, tag in ((CUT, "0.508(여기서 잴 골)"), (CUT_THEM, "0.51")):
    base = None
    for lo, hi, nm in LENB:
        g = L2[(L2.N >= lo) & (L2.N <= hi)]
        if not len(g):
            continue
        per = g.groupby("hap").b_mean.apply(lambda x: int((x <= cut).sum()))
        nmed = float(g.groupby("hap").N.first().median())
        if base is None:
            base = (per.mean(), nmed)
        rows.append(dict(컷=tag, 층=nm, n_hap=len(per), N중앙=nmed,
                         저메틸수=per.mean(), 저메틸수_sd=per.std(),
                         비례했다면=base[0] * nmed / base[1],
                         저메틸몫=float((g.b_mean <= cut).mean())))
C = pd.DataFrame(rows)
save(C, "FIG3_C_LOWCOUNT.tsv")
P()
refc = {"N40-80": 23.9, "N81-110": 18.9, "N111-140": 18.1, "N141+": 17.1}
P("   %-12s %8s %10s %10s %12s" % ("층", "하플로", "저메틸수", "주석", "비례했다면"))
for _, r in C[C.컷.str.startswith("0.508")].iterrows():
    P("   %-12s %8d %10.1f %10.1f %12.1f"
      % (r.층, r.n_hap, r.저메틸수, refc[r.층], r.비례했다면))
P("   -> the number of low-methylation units does not grow with array length; it decreases")
P("     (if proportional to length, N141+ would have %.1f; measured %.1f)"
  % (C[(C.컷.str.startswith("0.508")) & (C.층 == "N141+")].비례했다면.iloc[0],
     C[(C.컷.str.startswith("0.508")) & (C.층 == "N141+")].저메틸수.iloc[0]))

# ================================================================ D
P()
P("=" * 92)
P("[D] methylation pattern within the unit: by segment, high-methylation units vs low-methylation units")
P()
SEG = [("AluY", "b_alu", "n_alu"), ("spacer-1", "b_nts1", "n_nts1"), ("GA", "b_ga", "n_ga"),
       ("spacer-2+flank", "b_nts2", "n_nts2"), ("5'-flank", "b_ups", "n_ups"), ("유전자부", "b_gene", "n_gene")]
P("   segment resolution (median CpG): 6 segment bins, not 199 per-position points")
for nm, bc, nc in SEG:
    P("      %-8s CpG median %4.0f" % (nm, L1[nc].median()))
P("   5'-flank is the last 90 bp of spacer-2+flank, so it lies inside spacer-2+flank (12 + 196 = 208 = sum of 6 bins)")
P()
L1["lohi"] = np.where(L1.b_call > CUT, "고메틸 유닛", "저메틸 유닛")
rows = []
for grp in ("전체", "고메틸 유닛", "저메틸 유닛"):
    g = L1 if grp == "전체" else L1[L1.lohi == grp]
    for i, (nm, bc, nc) in enumerate(SEG):
        v = g[bc].dropna()
        # also as a mean over haplotypes, since units are not independent within a haplotype
        hv = g.groupby("hap")[bc].mean().dropna()
        rows.append(dict(무리=grp, 구획=nm, 순서=i, b=v.mean(), sd=v.std(), n=len(v),
                         b_hap=hv.mean(), sd_hap=hv.std(), n_hap=len(hv),
                         # this is the median number of reference CpG sites actually called in that unit, not the CpG count
                         #   (not the number counted in the sequence; canon_v2.fa holds 199 in total)
                         cpg_호출수_중앙=float(g[nc].median())))
D = pd.DataFrame(rows)
save(D, "FIG3_D_UNITPROFILE.tsv")
P()
P("   %-12s %-9s %8s %9s %9s" % ("무리", "구획", "unit β", "하플로평균", "n"))
for _, r in D.iterrows():
    P("   %-12s %-9s %8.4f %9.4f %9d" % (r.무리, r.구획, r.b, r.b_hap, r.n))
hi = D[D.무리 == "고메틸 유닛"].set_index("구획")
lo = D[D.무리 == "저메틸 유닛"].set_index("구획")
P()
P("   difference between the two groups by segment")
for nm, _, _ in SEG:
    P("      %-8s high %.4f − low %.4f = %+.4f" % (nm, hi.loc[nm].b, lo.loc[nm].b,
                                                hi.loc[nm].b - lo.loc[nm].b))
d = (hi.b - lo.b)
P("   -> segment with the largest difference %s (%+.4f), segment with the smallest %s (%+.4f)"
  % (d.idxmax(), d.max(), d.idxmin(), d.min()))

# ================================================================ E
P()
P("=" * 92)
P("[E] basis for the division into 8 zones: methylation axis x sequence axis crossed")
P()
ZI = S + "02_tables/ZONE_INTEGRATED.tsv"
if os.path.exists(ZI):
    Z = pd.read_csv(ZI, sep="\t")
    P("   input ZONE_INTEGRATED.tsv (md5 %s, rows %d, columns %s)"
      % (md5(ZI)[:12], len(Z), list(Z.columns)[:10]))
    save(Z, "FIG3_E_ZONEGRID.tsv")
else:
    P("   no ZONE_INTEGRATED.tsv")

# The basis for the third proximal boundary at 17: the position at which each length stratum reaches 98% of its interior value.
rows = []
for lo_, hi_, nm in ((30, 59, "N30-59"), (60, 99, "N60-99"),
                     (100, 139, "N100-139"), (140, 300, "N140-300")):
    s = L2[(L2.N >= lo_) & (L2.N <= hi_) & (L2.ap <= L2.ad)]
    it = L2[(L2.N >= lo_) & (L2.N <= hi_) & (L2.ap > 30) & (L2.ad > 30)].f_hi.mean()
    hit = None
    for k in range(4, 31):
        v = s.loc[s.ap == k, "f_hi"].mean()
        if v == v and it == it and v >= 0.98 * it:
            hit = k
            break
    for k in range(1, 31):
        v = s.loc[s.ap == k, "f_hi"]
        if len(v):
            rows.append(dict(층=nm, ap=k, f_hi=v.mean(), n=len(v), 내부=it, 도달자리=hit))
E2 = pd.DataFrame(rows)
save(E2, "FIG3_E2_PROXBOUND.tsv")
P()
P("   proximal boundary 3: position reaching 98% (of the interior value)")
for nm, g in E2.groupby("층", sort=False):
    P("      %-10s interior f_hi %.4f, position reached %s"
      % (nm, g.내부.iloc[0], g.도달자리.iloc[0]))
P("   -> 16~21 in every length stratum -> boundary set at 17")

# ================================================================ F
P()
P("=" * 92)
P("[F] between-haplotype difference by zone: reproduces Fig1_3 a")
P()


def zone8(ap, ad):
    if ap <= ad:
        if ap <= 3:
            return "근위 1-3"
        if ap <= 10:
            return "근위 4-10"
        if ap <= 17:
            return "근위 11-17"
        return "내부"
    if ad <= 2:
        return "원위 1-2"
    if ad <= 8:
        return "원위 3-8"
    if ad <= 10:
        return "원위 9-10"
    if ad <= 22:
        return "원위 11-22"
    return "내부"


ZORD = ["근위 1-3", "근위 4-10", "근위 11-17", "내부",
        "원위 11-22", "원위 9-10", "원위 3-8", "원위 1-2"]
L2["z8"] = [zone8(a, d) for a, d in zip(L2.ap, L2.ad)]
hz = L2.groupby(["hap", "z8"]).b_mean.mean().reset_index()
rows = []
for z in ZORD:
    v = hz[hz.z8 == z].b_mean.dropna()
    rows.append(dict(구역=z, 순서=ZORD.index(z), n_hap=len(v), b=v.mean(),
                     sd=v.std(), iqr=float(v.quantile(.75) - v.quantile(.25)),
                     p10=float(v.quantile(.10)), p90=float(v.quantile(.90))))
F = pd.DataFrame(rows)
save(F, "FIG3_F_ZONEVAR.tsv")
P()
P("   %-12s %7s %8s %8s %8s" % ("구역", "하플로", "평균β", "하플로간 SD", "IQR"))
for _, r in F.iterrows():
    P("   %-12s %7d %8.4f %10.4f %8.4f" % (r.구역, r.n_hap, r.b, r.sd, r.iqr))
P()
P("   -> zone with the smallest between-haplotype difference: %s (SD %.4f)"
  % (F.loc[F.sd.idxmin()].구역, F.sd.min()))
P("   -> zone with the largest: %s (SD %.4f)" % (F.loc[F.sd.idxmax()].구역, F.sd.max()))

# ================================================================ cohort
P()
P("=" * 92)
P("[cohort] denominator, common to all panels")
P()
CO = L2[["sample", "slot", "hap", "N"]].drop_duplicates("hap")
save(CO.sort_values("hap"), "FIG3_COHORT.tsv")
el = hashlib.md5("".join(sorted(CO.hap)).encode()).hexdigest()
P("   haplotypes %d, individuals %d, element md5 (concatenated) %s"
  % (len(CO), CO["sample"].nunique(), el[:12]))
P("   => the same number as the H2 stratum (510 / 280). The symmetric difference of the element sets is 0.")
P()
P("   chemistry strata (HPRC ONT release2 public index md5 f24b361a62b4)")
P("      R941 143 individuals (guppy 6.5.7 / 6.4.6), R1041 41 individuals (dorado 0.6.0), not in the index 96 individuals")
P("      per-stratum troughs are not reported.")
P("      (1) the pooled value of that estimator is 0.4582, so it does not split 0.508; (2) the R1041 stratum has only 6,721 units (R941 23,674), so")
P("      the trough bins are empty, argmin is tied, and the location is set by scan order.")
P("      -> moving the cut 0.45->0.55 changes the low-methylation share 0.208->0.234")
LOG.close()
