# -*- coding: utf-8 -*-
# FigX25 data preparation, per data source: 1) is the unit beta bimodal, 2) do the distal-to-proximal position profiles resemble each other.
# A copy of the same name is also in the repository: analysis/producers_v11_2026-10-11/figx25_shape_prep.py.
#   The name and the logic are the same; the only differences are the comments and how the paths are supplied
#   (DATA_ROOT and TABLE_DIR against _need, string concatenation against os.path.join).
#   The md5 of the output body, with the # header removed, is the same for both copies,
#     `FIGX25_BIMODAL` `38ff2f78943a` · `FIGX25_DENSITY` `4012e473da5a` ·
#     FIGX25_PROFILE at 2151a2cc4d17.
#   The md5 of the whole file always differs, because the .tsv header carries a date.
#   Run with the same TABLE_DIR, each copy overwrites the other's table, the names being identical.
#     So any test run must use its own TABLE_DIR.
#   Figure 3 and Supplementary Figure 3 read the table this script writes; the values are the same either way.
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
# The data roots come from environment variables.
#   They have no default: the script stops with a message if one is not set.
DATA_ROOT = os.environ["DATA_ROOT"]

# Where the tables are read and written. The tables are not deposited next to the code, so these paths come from the environment.
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
S = DATA_ROOT + "/"
FAMS = ["1kGP", "GIAB", "HPRC"]
EDGE = 20            # how many units from the end to plot


def fam(s):
    s = str(s)
    # The tumour sample is not included in the GIAB curve.
    #   GIAB_SOMATIC_ONT is the single individual HG008-T, and accounted for 237 of the 465 GIAB units (51%).
    #   It remains in supplementary panel F, normal against tumour, where that is its purpose.
    if s.startswith("GIAB_SOMATIC"):
        return "그 밖"
    return ("HPRC" if s.startswith("HPRC") else "GIAB" if s.startswith("GIAB")
            else "1kGP" if s.startswith("1kGP") else "그 밖")


R = pd.read_parquet(S + "02_tables/ru_v11.parquet",
                    columns=["sample", "slot", "data_source", "ap", "N", "b"])
R["fam"] = R.data_source.map(fam)
R = R[R.fam.isin(FAMS)].copy()
R["hap"] = R["sample"].astype(str) + "|" + R.slot.astype(str)

# HPRC is restricted to registry H2 (qc_eligible = 1).
#   The registry covers HPRC_ONT_submissions only, so there is no tier to apply to 1kGP or GIAB and those two are left unfiltered.
#   Filtering by name would be wrong: 14 of the 18 1kGP entries share a name with an H2 key, because the same individual appears in both sources.
_TI = pd.read_csv(S + "02_tables/TIERS_H0H3_SLOTS_2026-10-05.tsv", sep="	", dtype=str)
_H2 = set(_TI[_TI.tier == "H2"]["sample"] + "|" + _TI[_TI.tier == "H2"]["slot"])
_before = R[R.fam == "HPRC"].hap.nunique()
R = R[(R.fam != "HPRC") | (R.hap.isin(_H2))].copy()
print("[H2 관문] HPRC 하플로 %d -> %d (대장 H2 %d) · 1kGP·GIAB 는 층 없어 그대로"
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
                     interior_share=round(float(((v > 0.2) & (v < 0.8)).mean()), 4),
                     splithalf_골차_중앙=round(float(np.median(diff)), 4),
                     splithalf_골차_90p=round(float(np.percentile(diff, 90)), 4)))
    for x, yy in zip(G[::10], s["y"][::10]):
        curves.append(dict(fam=f, beta=round(float(x), 4), density=round(float(yy), 5)))
B = pd.DataFrame(rows)
B.to_csv(os.path.join(HERE, "FIGX25_BIMODAL.tsv"), sep="\t", index=False)
pd.DataFrame(curves).to_csv(os.path.join(HERE, "FIGX25_DENSITY.tsv"), sep="\t", index=False)
pd.set_option("display.width", 240)
print("## ⑴ 유닛 β 가 이봉인가 (자료원마다)")
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

print("   근위 끝을 쓸 수 있는 하플로타입 : " +
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
print("## ⑵ 원위 끝에서 센 자리별 평균 β (끝 1-10)")
print(PR[(PR.쪽 == "원위") & (PR.거리 <= 10)]
      .pivot_table(index="거리", columns="fam", values="mean").round(3).to_string())
print()
print("## 근위 끝에서 센 자리별 평균 β (끝 1-10)")
print(PR[(PR.쪽 == "근위") & (PR.거리 <= 10)]
      .pivot_table(index="거리", columns="fam", values="mean").round(3).to_string())
print()
print("## 가운데")
print(PR[PR.쪽 == "가운데"][["fam", "n_unit", "n_hap", "mean"]].round(3).to_string(index=False))

# Do the profiles resemble each other: a correlation per pair of data sources.
print()
print("## ⑶ 윤곽이 닮았나 (끝 1-20 의 평균 β 를 자료원끼리 상관)")
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
print("[끝] FIGX25_UNIT · FIGX25_BIMODAL · FIGX25_DENSITY · FIGX25_PROFILE")
