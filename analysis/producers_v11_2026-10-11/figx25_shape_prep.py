# -*- coding: utf-8 -*-
# FigX25 data preparation, per data source: 1) is the unit beta bimodal, 2) do the distal-to-proximal position profiles resemble each other.
#   
#
#   Only these two are examined: the absolute level cannot be separated from platform and basecaller confounding,
#   and 1kGP holds only one fully assembled array, so a comparison of means does not stand.
#   Therefore only the shape is used: bimodality and the per-position profile, both of which are less sensitive to a difference in level.
#
#   The orientation comes from measurement, not from the column name:
#     unit beta at ap = 1 is 0.030 (0.080 over the first 10 units), the other end is 0.354 and the interior is 0.815.
#     This matches the registry's hypomethylated end lengths of 14.58 distal and 2.56 proximal, so ap is the distance from the distal end
#     and ad = N - ap + 1 is the distance from the proximal end. In the G8 map, array_pos 1 is the distal end.
#   Ruler: ru_v11.b, the hard call. A unit's value is the mean over the reads covering it.
import io, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    The roots of the data this script reads are listed in the
    "Running this code" section of README.md."""
    import os as _o, sys as _s
    _x = _o.environ.get(_v)
    if not _x:
        _s.exit(
            "environment variable %s is not set; it is required by %s"
            % (_v, _o.path.basename(__file__))
            + "\n(a module stops the script that imports it)."
            + "\nThe deposited copy has no default; see README.md"
            + " section \"Running this code\"."
        )
    return _x


HERE = os.path.dirname(os.path.abspath(__file__))
S = _need("RDNA_V11_ROOT")
FAMS = ["1kGP", "GIAB", "HPRC"]
EDGE = 20            # how many units from the end to plot


def fam(s):
    # `GIAB_SOMATIC_ONT` is excluded.
    #    It is the tumour sample HG008-T, whose name begins with GIAB.
    #    237 of the 465 units under that prefix, 51%, come from that one individual.
    # The other outputs exclude that sample as well (`arrayprof_noHG008T_…py`).
    #    Without it: 3 individuals (hg002, hg003, hg004), 4 slots, 228 units.
    s = str(s)
    if s.startswith("GIAB_SOMATIC"):
        return "그 밖"
    return ("HPRC" if s.startswith("HPRC") else "GIAB" if s.startswith("GIAB")
            else "1kGP" if s.startswith("1kGP") else "그 밖")


#  The paths are built with os.path.join, so the root needs no trailing slash.
R = pd.read_parquet(os.path.join(S, "02_tables", "ru_v11.parquet"),
                    columns=["sample", "slot", "data_source", "ap", "N", "b"])
R["fam"] = R.data_source.map(fam)
R = R[R.fam.isin(FAMS)].copy()
R["hap"] = R["sample"].astype(str) + "|" + R.slot.astype(str)

# HPRC is restricted to registry H2 (qc_eligible = 1).
#   The registry covers HPRC_ONT_submissions only, so there is no tier to apply to 1kGP or GIAB and those two are left unfiltered.
#   Filtering by name would be wrong: 14 of the 18 1kGP entries share a name with an H2 key, because the same individual appears in both sources.
_TI = pd.read_csv(os.path.join(S, "02_tables", "TIERS_H0H3_SLOTS_2026-10-05.tsv"),
                  sep="	", dtype=str)
_H2 = set(_TI[_TI.tier == "H2"]["sample"] + "|" + _TI[_TI.tier == "H2"]["slot"])
_before = R[R.fam == "HPRC"].hap.nunique()
R = R[(R.fam != "HPRC") | (R.hap.isin(_H2))].copy()
print("[H2 gate] HPRC haplotypes %d -> %d (registry H2 %d); 1kGP and GIAB have no tier, kept as is"
      % (_before, R[R.fam == "HPRC"].hap.nunique(), len(_H2)))
U = R.groupby(["fam", "sample", "hap", "ap", "N"], as_index=False).b.mean()
U["ad"] = U.N - U.ap + 1
U["rel"] = (U.ap - 1) / (U.N - 1).replace(0, np.nan)      # 0 = distal end, 1 = proximal end
U.to_csv(os.path.join(HERE, "FIGX25_UNIT.tsv"), sep="\t", index=False, float_format="%.5f")
print("유닛 %d · " % len(U) +
      " · ".join("%s 하플로 %d/개체 %d"
                 % (f, U[U.fam == f].hap.nunique(), U[U.fam == f]["sample"].nunique())
                 for f in FAMS))
print()

# ---------- 1) is it bimodal ----------
G = np.linspace(0.01, 0.99, 4901)


def peaks(y):
    """local maxima (positions higher than both neighbours), ordered by height"""
    i = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    return sorted(i, key=lambda k: -y[k])


def shape(v, bw=0.03):
    y = gaussian_kde(v, bw_method=bw)(G)
    pk = peaks(y)
    big = [k for k in pk if y[k] >= 0.05 * y[pk[0]]]      # only peaks at least 5% as high as the tallest are counted
    lo = G[int(np.argmax(np.where(G < 0.4, y, -1)))]
    hi = G[int(np.argmax(np.where(G > 0.6, y, -1)))]
    m = (G > lo) & (G < hi)
    j = int(np.argmin(np.where(m, y, 1e18)))
    return dict(봉우리수=len(big), 저봉=round(float(lo), 3), 고봉=round(float(hi), 3),
                골=round(float(G[j]), 4),
                골봉비=round(float(y[j] / max(y[int(np.argmax(np.where(G < 0.4, y, -1)))],
                                            y[int(np.argmax(np.where(G > 0.6, y, -1)))])), 4),
                y=y)


rows, curves = [], []
rng = np.random.default_rng(20261008)
for f in FAMS:
    d = U[U.fam == f]
    v = d.b.dropna().values
    s = shape(v)
    # split-half: the haplotypes are halved to test whether the trough reproduces (50 replicates)
    haps = d.hap.unique()
    by = {h: x.b.dropna().values for h, x in d.groupby("hap")}
    sh = []
    for _ in range(50):
        p = rng.permutation(haps)
        a = np.concatenate([by[h] for h in p[:len(p) // 2]])
        b2 = np.concatenate([by[h] for h in p[len(p) // 2:]])
        try:
            sh.append((shape(a)["골"], shape(b2)["골"]))
        except Exception:
            pass
    sh = np.array(sh)
    diff = np.abs(sh[:, 0] - sh[:, 1]) if len(sh) else np.array([np.nan])
    rows.append(dict(자료원=f, 하플로타입=len(haps), 개체=d["sample"].nunique(), 유닛=len(v),
                     전체평균=round(float(v.mean()), 4),
                     봉우리수=s["봉우리수"], 저봉=s["저봉"], 고봉=s["고봉"], 골=s["골"],
                     골봉비=s["골봉비"],
                     양끝몫=round(float(((v <= 0.2) | (v >= 0.8)).mean()), 4),
                     가운데몫=round(float(((v > 0.2) & (v < 0.8)).mean()), 4),
                     splithalf_골차_중앙=round(float(np.median(diff)), 4),
                     splithalf_골차_90p=round(float(np.percentile(diff, 90)), 4)))
    for x, yy in zip(G[::10], s["y"][::10]):
        curves.append(dict(fam=f, beta=round(float(x), 4), density=round(float(yy), 5)))
B = pd.DataFrame(rows)
#  FIGX25_BIMODAL.tsv is written with English column names (the source values were already
#    English). B keeps its names for the printout below. Reader : figS4_A_unit_beta_three_sources.py
B.rename(columns={"자료원": "source", "하플로타입": "haplotypes", "개체": "individuals", "유닛": "units",
                  "전체평균": "mean_beta", "봉우리수": "n_peaks", "저봉": "low_peak", "고봉": "high_peak",
                  "골": "valley", "골봉비": "valley_to_peak_ratio", "양끝몫": "fraction_at_ends",
                  "가운데몫": "fraction_in_middle", "splithalf_골차_중앙": "splithalf_valley_difference_median",
                  "splithalf_골차_90p": "splithalf_valley_difference_p90"}
         ).to_csv(os.path.join(HERE, "FIGX25_BIMODAL.tsv"), sep="\t", index=False)
pd.DataFrame(curves).to_csv(os.path.join(HERE, "FIGX25_DENSITY.tsv"), sep="\t", index=False)
pd.set_option("display.width", 240)
print("## (1) Is unit β bimodal (per data source)")
print(B.to_string(index=False))
print()

# ---------- 2) the distal-to-proximal position profile ----------
#   The distal end is binned by ap and the proximal end by ad. The proximal end requires knowing where the array ends,
#   so only haplotypes whose largest observed position equals N are used; otherwise the end is unknown.
# obs_max must be grouped by data source as well as hap (= sample|slot). Grouping by hap alone means that
#   when the same individual is in both 1kGP and HPRC, rows from the other source supply that individual's array end.
#   That inflates the 1kGP haplotypes usable at the proximal end to 17 of 18; the correct count is 6 of 18.
#   A haplotype observed all the way to the end cannot grow in count with distance, so a count that rises is a sign of haplotypes missing their proximal terminal unit.
obs = U.groupby(["fam", "hap"]).ap.max().rename("obs_max").reset_index()
U2 = U.merge(obs, on=["fam", "hap"])
U2["_ok"] = U2.obs_max == U2.N

print("   haplotypes whose proximal end is usable: " +
      " · ".join("%s %d / %d" % (f, U2[(U2.fam == f) & U2._ok].hap.nunique(),
                                 U2[U2.fam == f].hap.nunique()) for f in FAMS))

prof = []
for f in FAMS:
    d = U2[U2.fam == f]
    for k in range(1, EDGE + 1):
        v = d[d.ap == k].b
        prof.append(dict(fam=f, 쪽="원위", 거리=k, n_unit=len(v), n_hap=d[d.ap == k].hap.nunique(),
                         mean=float(v.mean()) if len(v) else np.nan,
                         median=float(v.median()) if len(v) else np.nan))
    dp = d[d._ok]
    for k in range(1, EDGE + 1):
        v = dp[dp.ad == k].b
        prof.append(dict(fam=f, 쪽="근위", 거리=k, n_unit=len(v), n_hap=dp[dp.ad == k].hap.nunique(),
                         mean=float(v.mean()) if len(v) else np.nan,
                         median=float(v.median()) if len(v) else np.nan))
    mid = d[(d.ap > EDGE) & (d.ad > EDGE)].b
    prof.append(dict(fam=f, 쪽="가운데", 거리=0, n_unit=len(mid),
                     n_hap=d[(d.ap > EDGE) & (d.ad > EDGE)].hap.nunique(),
                     mean=float(mid.mean()), median=float(mid.median())))
PR = pd.DataFrame(prof)
PR.to_csv(os.path.join(HERE, "FIGX25_PROFILE.tsv"), sep="\t", index=False, float_format="%.5f")
print()
print("## (2) Mean β by position counted from the distal end (end 1-10)")
print(PR[(PR.쪽 == "원위") & (PR.거리 <= 10)]
      .pivot_table(index="거리", columns="fam", values="mean").round(3).to_string())
print()
print("## Mean β by position counted from the proximal end (end 1-10)")
print(PR[(PR.쪽 == "근위") & (PR.거리 <= 10)]
      .pivot_table(index="거리", columns="fam", values="mean").round(3).to_string())
print()
print("## Interior")
print(PR[PR.쪽 == "가운데"][["fam", "n_unit", "n_hap", "mean"]].round(3).to_string(index=False))

# Do the profiles resemble each other: a correlation per pair of data sources.
print()
print("## (3) Are the profiles similar (correlation between data sources of mean β at end 1-20)")
for side in ("원위", "근위"):
    w = PR[(PR.쪽 == side) & (PR.거리 <= EDGE)].pivot_table(index="거리", columns="fam",
                                                          values="mean")
    out = []
    for a in range(len(FAMS)):
        for b2 in range(a + 1, len(FAMS)):
            x, y = FAMS[a], FAMS[b2]
            m = w[[x, y]].dropna()
            out.append("%s-%s rho %+.3f (n %d)" % (x, y, m[x].corr(m[y], method="spearman"), len(m)))
    print("   %-4s : %s" % (side, " · ".join(out)))
print()
print("[done] FIGX25_UNIT, FIGX25_BIMODAL, FIGX25_DENSITY, FIGX25_PROFILE")
