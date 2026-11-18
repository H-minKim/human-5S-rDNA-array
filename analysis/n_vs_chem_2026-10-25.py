# -*- coding: utf-8 -*-
#  The question: is the 1kGP curve unstable because the chemistries are mixed, or because there are few units?
#  Splitting by chemistry alone cannot separate the two.
#    The clue: R9 alone (719 units) gives 23 peaks, more than the pooled set (1,265 units, 16 peaks).
# If chemistry drove it, splitting should make it cleaner, and it does the opposite, so `n` is the suspect.
#  The ruler: draw 1,265 units 50 times from HPRC, which also mixes chemistries but holds 47,483 units,
#        and count the peaks and the trough with the same shape() used in figx25.
#        On the full set HPRC is stable at 3 peaks and a trough of 0.5068; if a small draw stays stable, n is not the cause.
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
print("HPRC units %s" % format(len(H), ","))
np.random.seed(20261025)
print()
print("## HPRC subsampled to small sizes: same measure as `figx25`")
print()
print("| units sampled | median number of peaks | range of number of peaks | median trough | trough range |")
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
          % (format(n, ","), " (full set)" if n == len(H) else "",
             int(np.median(pk)), min(pk), max(pk),
             float(np.median(vl)), min(vl), max(vl)))
print()
print("compare to 1kGP pooled  : peaks 16, trough 0.1860 (1,265 units)")
print("           1kGP R9      : peaks 23, trough 0.3924 (  719 units)")
print("           1kGP R10     : peaks  6, trough 0.5192 (  546 units)")
print("           GIAB pooled  : peaks  7, trough 0.5810 (  465 units)")
