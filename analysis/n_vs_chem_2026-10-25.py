# -*- coding: utf-8 -*-
#  물음 : `1kGP` 곡선이 흔들리는 것이 화학 혼합 때문인가, 유닛 수가 적어서인가.
#  제가 어제 화학 탓으로 적었는데 그 측정은 둘을 가르지 못합니다.
#    단서 : R9 만(719 유닛)이 봉우리 23 개로 합친 것(1,265 유닛 · 16 개)보다 더 많습니다.
# If chemistry drove it, splitting should make it cleaner, and it does the opposite, so `n` is the suspect.
#  자 : HPRC(화학이 섞였지만 47,483 유닛)에서 1,265 유닛씩 50 번 뽑아
#        `figx25` 의 같은 `shape()` 로 봉우리 수·골을 센다.
#        HPRC 는 전수로는 봉우리 3 · 골 0.5068 로 단단하다 — 작게 뽑아도 그러면 n 탓이 아니다.
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


V11 = _need("RDNA_V11_READVAR_DIR")
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
    return len(big), float(G[j])


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
                  sep=chr(9), dtype=str)
_H2 = set(_TI[_TI.tier == "H2"]["sample"] + "|" + _TI[_TI.tier == "H2"]["slot"])
R = R[(R.fam != "HPRC") | (R.hap.isin(_H2))].copy()
U = R.groupby(["fam", "sample", "hap", "ap", "N"], as_index=False).b.mean()
H = U[U.fam == "HPRC"].b.values
print("HPRC 유닛 %s" % format(len(H), ","))
np.random.seed(20261025)
print()
print("## HPRC 를 작게 뽑아 보면 — `figx25` 의 같은 자")
print()
print("| 뽑은 유닛 수 | 봉우리수 중앙 | 봉우리수 범위 | 골 중앙 | 골 범위 |")
print("|---|--:|--:|--:|--:|")
for n in (1265, 719, 546, 465, 228, 5000, 47483):
    if n > len(H):
        continue
    pk, vl = [], []
    reps = 1 if n == len(H) else 50
    for _ in range(reps):
        s = H if n == len(H) else np.random.choice(H, n, replace=False)
        a, b = shape(s)
        pk.append(a)
        vl.append(b)
    print("| %s%s | %d | %d - %d | %.4f | %.3f - %.3f |"
          % (format(n, ","), " (전수)" if n == len(H) else "",
             int(np.median(pk)), min(pk), max(pk),
             float(np.median(vl)), min(vl), max(vl)))
print()
print("견줄 값 — 1kGP 합친 것 : 봉우리 16 · 골 0.1860 (1,265 유닛)")
print("           1kGP R9      : 봉우리 23 · 골 0.3924 (  719 유닛)")
print("           1kGP R10     : 봉우리  6 · 골 0.5192 (  546 유닛)")
print("           GIAB 합친 것 : 봉우리  7 · 골 0.5810 (  465 유닛)")
