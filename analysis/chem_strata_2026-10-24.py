# -*- coding: utf-8 -*-
#  Re-measurement stratified by R9 / R10
#
# Why : the basecaller version shifts beta by 0.13, and
#    the 10 `1kGP` individuals of `S4 A` are of mixed chemistry.
#
#  How: each panel's own ruler is used unchanged and only the stratum is split out. No new ruler is made.
#    (1) `S4 A`: the recipe of `figx25_shape_prep.py` (H2 gate; unit beta = mean of b; `shape()`)
#    (2) the valley of `Fig 4 a and b`: the gate of `valley_hardcall.py`
#       (HPRC_ONT_submissions; n_cpg >= 20; H2 slots; ao/jx deduplicated with jx preferred; hard-call beta)
# The HPRC release 2 ONT sequencing index is read from RDNA_HPRC_R2_ONT_INDEX.
# Set RDNA_HPRC_R2_ONT_INDEX to a local copy; without it the script falls back to a file of
# that name beside itself.
import collections
import csv
import gzip
import io
import os

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def _need(k):
    """The deposited copy keeps no default: it stops rather than read the wrong place."""
    v = os.environ.get(k)
    if not v:
        raise SystemExit("environment variable %s is not set (the "
                         "deposited copy has no default; see README.md,"
                         " section 'Environment variables')" % k)
    return v


SP = os.path.dirname(os.path.abspath(__file__))
NL, TAB = chr(10), chr(9)
V11 = _need("RDNA_V11_READVAR_DIR")
W = _need("RDNA_METH_BLOCKS_DIR")
L1P = os.path.join(V11, "02_tables", "L1_READ_UNIT.parquet")

#  --------- chemistry map (both taken from the measured files)
hp = {}
for r in csv.DictReader(io.open(os.environ.get("RDNA_HPRC_R2_ONT_INDEX",
                                     os.path.join(SP, "hprc_ont_r2_index_2026-10-24.csv")),
                                encoding="utf8", errors="replace")):
    c = (r.get("sequencing_chemistry") or "").strip()
    if c:
        hp.setdefault(r["sample_id"], set()).add(c)
HPCHEM = {k: (sorted(v)[0] if len(v) == 1 else "mixed") for k, v in hp.items()}
kg = {}
for r in csv.DictReader(io.open(os.path.join(V11, "03_figures", "paper_Fig", "00_figure_ledger",
                                             "source_evidence",
                                             "1kGP_LRSC_500_ONT_Metadata_2026-10-09.tsv"),
                                encoding="utf8", errors="replace"), delimiter=TAB):
    kg[r["NHGRI_ID"]] = (r.get("ONT_pore") or "?").strip()
print("chemistry map: HPRC index %d individuals, 1kGP metadata table %d individuals" % (len(HPCHEM), len(kg)))
print()

#  --------- (1) `S4 A` - the figx25 recipe unchanged, with the stratum added
G = np.linspace(0.01, 0.99, 4901)


def peaks(y):
    i = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    return sorted(i, key=lambda k: -y[k])


def shape(v, bw=0.03):
    y = gaussian_kde(v, bw_method=bw)(G)
    pk = peaks(y)
    big = [k for k in pk if y[k] >= 0.05 * y[pk[0]]]
    lo = G[int(np.argmax(np.where(G < 0.4, y, -1)))]
    hi = G[int(np.argmax(np.where(G > 0.6, y, -1)))]
    m = (G > lo) & (G < hi)
    j = int(np.argmin(np.where(m, y, 1e18)))
    ylo = y[int(np.argmax(np.where(G < 0.4, y, -1)))]
    yhi = y[int(np.argmax(np.where(G > 0.6, y, -1)))]
    return dict(봉우리수=len(big), 저봉=round(float(lo), 3), 고봉=round(float(hi), 3),
                골=round(float(G[j]), 4), 골봉비=round(float(y[j] / max(ylo, yhi)), 4))


def famof(s):
    s = str(s)
    return ("HPRC" if s.startswith("HPRC") else "GIAB" if s.startswith("GIAB")
            else "1kGP" if s.startswith("1kGP") else "그 밖")


R = pd.read_parquet(os.path.join(V11, "02_tables", "ru_v11.parquet"),
                    columns=["sample", "slot", "data_source", "ap", "N", "b"])
R["fam"] = R.data_source.map(famof)
R = R[R.fam.isin(["1kGP", "GIAB", "HPRC"])].copy()
R["hap"] = R["sample"].astype(str) + "|" + R.slot.astype(str)
_TI = pd.read_csv(os.path.join(V11, "02_tables", "TIERS_H0H3_SLOTS_2026-10-05.tsv"),
                  sep=TAB, dtype=str)
_H2 = set(_TI[_TI.tier == "H2"]["sample"] + "|" + _TI[_TI.tier == "H2"]["slot"])
R = R[(R.fam != "HPRC") | (R.hap.isin(_H2))].copy()
U = R.groupby(["fam", "sample", "hap", "ap", "N", "data_source"], as_index=False).b.mean()


def stratum(row):
    f, s, d = row.fam, str(row["sample"]), str(row.data_source)
    if f == "1kGP":
        return kg.get(s, "unknown")
    if f == "GIAB":
        return "tumour(HG008-T)" if d == "GIAB_SOMATIC_ONT" else "normal(hg002-4)"
    return HPCHEM.get(s, "unknown (outside the index)")


U["층"] = U.apply(stratum, axis=1)
print("## (1) `SuppFig S4 A`: data source × chemistry/sample stratum  (same measure as figx25)")
print()
print("| data source | stratum | haplotypes | individuals | units | mean β | number of peaks | low peak | high peak | trough | trough/peak ratio |")
print("|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
rows = []
for f in ["1kGP", "GIAB", "HPRC"]:
    sub = U[U.fam == f]
    groups = [("all", sub)] + [(k, sub[sub["층"] == k]) for k in sorted(sub["층"].unique())]
    for k, g in groups:
        if len(g) < 60:
            print("| %s | %s | %d | %d | %d | %.4f | too few units to measure | | | | |"
                  % (f, k, g.hap.nunique(), g["sample"].nunique(), len(g),
                     g.b.mean() if len(g) else float("nan")))
            continue
        sh = shape(g.b.values)
        print("| %s | %s | %d | %d | %s | %.4f | %d | %.3f | %.3f | %.4f | %.4f |"
              % (f, k, g.hap.nunique(), g["sample"].nunique(), format(len(g), ","),
                 g.b.mean(), sh["봉우리수"], sh["저봉"], sh["고봉"], sh["골"], sh["골봉비"]))
        rows.append([f, k, g.hap.nunique(), g["sample"].nunique(), len(g),
                     round(float(g.b.mean()), 4), sh["봉우리수"], sh["저봉"], sh["고봉"],
                     sh["골"], sh["골봉비"]])
print()

#  --------- (2) the valley of `Fig 4 a and b` - the valley_hardcall gate unchanged, split by stratum
h2 = pd.read_parquet(L1P, columns=["sample", "slot", "gate"])
H2 = set(map(tuple, h2.loc[h2.gate == 1, ["sample", "slot"]].drop_duplicates()
             .astype(str).values))
best = {}
RU = os.path.join(W, "out_v11_fix4", "RU.tsv.gz")
for i, L in enumerate(gzip.open(RU, "rt", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split(TAB)
    if f[10] != "HPRC_ONT_submissions" or int(f[7]) < 20:
        continue
    if (f[0], f[1]) not in H2:
        continue
    k = (f[0], f[1], f[2], f[3])
    if k not in best or (f[6] == "jx" and best[k][0] != "jx"):
        best[k] = (f[6], float(f[9]), f[0])
print("## (2) trough of `Fig 4 a, b`: same gate as valley_hardcall, split by chemistry stratum only")
print()
by = collections.defaultdict(list)
for v in best.values():
    by[HPCHEM.get(v[2], "unknown (outside the index)")].append(v[1])
    by["all"].append(v[1])


def valley(b, bw=0.002, smooth=5, lo=0.20, hi=0.80):
    b = np.asarray(b)
    n = int(round(1.0 / bw))
    hist, edges = np.histogram(b, bins=n, range=(0.0, 1.0))
    dens = hist / hist.sum()
    mid = (edges[:-1] + edges[1:]) / 2.0
    d = np.convolve(dens, np.ones(smooth) / smooth, mode="same")
    m = (mid >= lo) & (mid <= hi)
    j = int(np.argmin(d[m]))
    peak_lo = mid[mid < 0.4][int(np.argmax(d[mid < 0.4]))]
    peak_hi = mid[mid > 0.6][int(np.argmax(d[mid > 0.6]))]
    return mid[m][j], peak_lo, peak_hi, float(d[m][j] / max(d[mid < 0.4].max(),
                                                            d[mid > 0.6].max()))


print("| stratum | individuals | read-units | trough | low peak | high peak | trough/peak | difference from 0.508 |")
print("|---|--:|--:|--:|--:|--:|--:|--:|")
for k in ["all"] + sorted(x for x in by if x != "all"):
    vals = by[k]
    ns = len({v[2] for v in best.values() if HPCHEM.get(v[2], "unknown (outside the index)") == k}) \
        if k != "all" else len({v[2] for v in best.values()})
    if len(vals) < 500:
        print("| %s | %d | %s | too few to measure | | | | |" % (k, ns, format(len(vals), ",")))
        continue
    v, pl, ph, r = valley(vals)
    print("| %s | %d | %s | **%.4f** | %.3f | %.3f | %.4f | %+.4f |"
          % (k, ns, format(len(vals), ","), v, pl, ph, r, v - 0.508))
    rows.append(["H2-trough", k, "", ns, len(vals), round(float(v), 4), "", round(float(pl), 3),
                 round(float(ph), 3), "", round(r, 4)])

P = os.path.join(SP, "chem_strata_2026-10-24.tsv")
out = [TAB.join(["block", "stratum", "haplotypes", "individuals", "units_or_readunits",
                 "mean_or_valley", "n_peaks", "peak_low", "peak_high", "valley", "ratio"])]
for r in rows:
    out.append(TAB.join(str(x) for x in r))
io.open(P, "w", encoding="utf8", newline=NL).write(NL.join(out) + NL)
print()
print("wrote: %s" % P)
