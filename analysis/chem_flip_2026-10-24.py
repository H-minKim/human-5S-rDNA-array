# -*- coding: utf-8 -*-
#  The per-stratum valley came out at 0.4470 on R10 and 0.4950 on R9. So the question is:
#    how many calls move if one cut (0.508) is applied to both chemistries.
#  Ruler: a deep valley means few units move when the cut shifts a little, so that fraction is counted.
#    The fraction of read-units inside [min(valleys), 0.508] = those that could flip wherever in between the cut is placed
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

hp = {}
for r in csv.DictReader(io.open(os.environ.get("RDNA_HPRC_R2_ONT_INDEX",
                                     os.path.join(SP, "hprc_ont_r2_index_2026-10-24.csv")),
                                encoding="utf8", errors="replace")):
    c = (r.get("sequencing_chemistry") or "").strip()
    if c:
        hp.setdefault(r["sample_id"], set()).add(c)
CH = {k: (sorted(v)[0] if len(v) == 1 else "혼재") for k, v in hp.items()}

h2 = pd.read_parquet(os.path.join(V11, "02_tables", "L1_READ_UNIT.parquet"),
                     columns=["sample", "slot", "gate"])
H2 = set(map(tuple, h2.loc[h2.gate == 1, ["sample", "slot"]].drop_duplicates()
             .astype(str).values))
best = {}
for i, L in enumerate(gzip.open(os.path.join(W, "out_v11_fix4", "RU.tsv.gz"),
                                "rt", encoding="utf8")):
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

by = collections.defaultdict(list)
for v in best.values():
    by[CH.get(v[2], "모름(색인밖)")].append(v[1])
    by["전체"].append(v[1])

LOW, CUT = 0.4470, 0.508          # The R10 valley and the published cut
print("## 컷을 0.4470(R10 골) ~ 0.508(발행 컷) 사이로 옮길 때 뒤집힐 수 있는 read-유닛의 몫")
print()
print("| 층 | read-유닛 | [0.4470, 0.508) 안 | 그 몫 | 0.508 아래 몫 | 0.4470 아래 몫 | 차 |")
print("|---|--:|--:|--:|--:|--:|--:|")
rows = []
for k in ["전체", "R941", "R1041", "모름(색인밖)"]:
    if k not in by:
        continue
    b = np.asarray(by[k])
    band = int(((b >= LOW) & (b < CUT)).sum())
    f_cut = float((b < CUT).mean())
    f_low = float((b < LOW).mean())
    print("| %s | %s | %s | **%.3f%%** | %.3f%% | %.3f%% | %+.3f%%p |"
          % (k, format(len(b), ","), format(band, ","), 100.0 * band / len(b),
             100.0 * f_cut, 100.0 * f_low, 100.0 * (f_cut - f_low)))
    rows.append([k, len(b), band, round(100.0 * band / len(b), 3),
                 round(100.0 * f_cut, 3), round(100.0 * f_low, 3)])
P = os.path.join(SP, "chem_flip_2026-10-24.tsv")
out = [TAB.join(["stratum", "read_units", "in_band_0.4470_0.508", "pct_in_band",
                 "pct_below_0.508", "pct_below_0.4470"])]
for r in rows:
    out.append(TAB.join(str(x) for x in r))
io.open(P, "w", encoding="utf8", newline=NL).write(NL.join(out) + NL)
print()
print("기록 : %s" % P)
