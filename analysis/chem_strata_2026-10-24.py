# -*- coding: utf-8 -*-
#  R9/R10 층별 재측정
#
# Why : the basecaller version shifts beta by 0.13,
#    기록에 *"R9/R10 층을 나눈 재측정 미실시 · Fig 4 a·b 에 걸린다"* 가 남아 있다.
# The 10 `1kGP` individuals of `S4 A` turned out to be of mixed chemistry.
#
#  어떻게 — 각 패널이 자기 자로 잰 것을 그대로 쓰고 층만 가른다. 새 자를 만들지 않는다.
#    ⑴ `S4 A` : `figx25_shape_prep.py` 의 그 recipe (H2 관문 · 유닛 β = b 평균 · `shape()`)
#    ⑵ `Fig 4 a·b` 의 골 : `valley_hardcall.py` 의 그 관문
#       (HPRC_ONT_submissions · n_cpg>=20 · H2 슬롯 · ao/jx 중복제거 jx 우선 · 하드콜 β)
#  읽기만 한다. 쓰는 것은 모래상자뿐.
# The HPRC release 2 ONT sequencing index this script reads is not part of this archive.
# Set RDNA_HPRC_R2_ONT_INDEX to a local copy; without it the script falls back to a file of
# that name beside itself. The README says where the index comes from and gives its md5.
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

#  ───────── 화학 지도 (둘 다 실측 파일에서)
hp = {}
for r in csv.DictReader(io.open(os.environ.get("RDNA_HPRC_R2_ONT_INDEX",
                                     os.path.join(SP, "hprc_ont_r2_index_2026-10-24.csv")),
                                encoding="utf8", errors="replace")):
    c = (r.get("sequencing_chemistry") or "").strip()
    if c:
        hp.setdefault(r["sample_id"], set()).add(c)
HPCHEM = {k: (sorted(v)[0] if len(v) == 1 else "혼재") for k, v in hp.items()}
kg = {}
for r in csv.DictReader(io.open(os.path.join(V11, "03_figures", "paper_Fig", "00_figure_ledger",
                                             "source_evidence",
                                             "1kGP_LRSC_500_ONT_Metadata_2026-10-09.tsv"),
                                encoding="utf8", errors="replace"), delimiter=TAB):
    kg[r["NHGRI_ID"]] = (r.get("ONT_pore") or "?").strip()
print("화학 지도 : HPRC 색인 %d 개체 · 1kGP 메타표 %d 개체" % (len(HPCHEM), len(kg)))
print()

#  ───────── ⑴ `S4 A` — figx25 의 recipe 그대로, 층만 더한다
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
        return kg.get(s, "모름")
    if f == "GIAB":
        return "종양(HG008-T)" if d == "GIAB_SOMATIC_ONT" else "정상(hg002-4)"
    return HPCHEM.get(s, "모름(색인밖)")


U["층"] = U.apply(stratum, axis=1)
print("## ⑴ `SuppFig S4 A` — 자료원 × 화학/시료 층  (figx25 의 자 그대로)")
print()
print("| 자료원 | 층 | 하플로 | 개체 | 유닛 | 평균 β | 봉우리수 | 저봉 | 고봉 | 골 | 골봉비 |")
print("|---|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|")
rows = []
for f in ["1kGP", "GIAB", "HPRC"]:
    sub = U[U.fam == f]
    groups = [("전체", sub)] + [(k, sub[sub["층"] == k]) for k in sorted(sub["층"].unique())]
    for k, g in groups:
        if len(g) < 60:
            print("| %s | %s | %d | %d | %d | %.4f | 유닛 적어 못 잼 | | | | |"
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

#  ───────── ⑵ `Fig 4 a·b` 의 골 — valley_hardcall 의 관문 그대로, 층만 가른다
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
print("## ⑵ `Fig 4 a·b` 의 골 — valley_hardcall 의 관문 그대로, 화학 층만 가름")
print()
by = collections.defaultdict(list)
for v in best.values():
    by[HPCHEM.get(v[2], "모름(색인밖)")].append(v[1])
    by["전체"].append(v[1])


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


print("| 층 | 개체 | read-유닛 | 골 | 저봉 | 고봉 | 골/봉 | 0.508 과 차 |")
print("|---|--:|--:|--:|--:|--:|--:|--:|")
for k in ["전체"] + sorted(x for x in by if x != "전체"):
    vals = by[k]
    ns = len({v[2] for v in best.values() if HPCHEM.get(v[2], "모름(색인밖)") == k}) \
        if k != "전체" else len({v[2] for v in best.values()})
    if len(vals) < 500:
        print("| %s | %d | %s | 적어 못 잼 | | | | |" % (k, ns, format(len(vals), ",")))
        continue
    v, pl, ph, r = valley(vals)
    print("| %s | %d | %s | **%.4f** | %.3f | %.3f | %.4f | %+.4f |"
          % (k, ns, format(len(vals), ","), v, pl, ph, r, v - 0.508))
    rows.append(["H2-골", k, "", ns, len(vals), round(float(v), 4), "", round(float(pl), 3),
                 round(float(ph), 3), "", round(r, 4)])

P = os.path.join(SP, "chem_strata_2026-10-24.tsv")
out = [TAB.join(["block", "stratum", "haplotypes", "individuals", "units_or_readunits",
                 "mean_or_valley", "n_peaks", "peak_low", "peak_high", "valley", "ratio"])]
for r in rows:
    out.append(TAB.join(str(x) for x in r))
io.open(P, "w", encoding="utf8", newline=NL).write(NL.join(out) + NL)
print()
print("기록 : %s" % P)
