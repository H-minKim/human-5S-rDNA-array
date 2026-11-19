# -*- coding: utf-8 -*-
# A per-array summary of the three Arabidopsis contexts (CG, CHG, CHH), written to tables/FIGDATA_SUPP4_ATH_CONTEXT.tsv.
#
#   Three numbers from this table (CG 0.844, CHG 0.423, CHH 0.146) appear in the manuscript text and legend.
#
#   Ruler: the mean beta is call-weighted, that is sum(n_5mC) / sum(n_calls) over every read x unit in the array.
#     It is a different number from the arithmetic mean of the read-unit betas, which gives CG 0.841.
#     The unit of aggregation differs, so the two must not be mixed; the manuscript uses the value this script produces, 0.844.
#
#   Order: by chromosome number, then by the array's start coordinate. This is the order of the data, not the layout of the figure,
#     and it appears only in the three result tables PERARRAY_META, CHROM_NAMES and UNITBETA_PERARRAY_v13.
#
#   Inputs: tables/ATH_READUNIT_CONTEXT.tsv.gz and tables/ATH_CONTEXT_PERARRAY.tsv, both from x13_ath_context.py.
#           tables/UNITBETA_PERARRAY_v13.tsv.gz (`x9D`), tables/PERARRAY_META.tsv (`x34`), tables/CHROM_NAMES.tsv (`x9Q`)
#   Output: tables/FIGDATA_SUPP4_ATH_CONTEXT.tsv
import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'Environment variables')" % k)
    return v


import os
import sys

import numpy as np
import pandas as pd

NEW = _need("RDNA_OUT_DIR")
sys.path.insert(0, NEW + "/code")
import arch        # noqa: E402
import figorder    # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
REASON = "애기장대 세 문맥 array 요약 — 그리는 글에서 분석 글로 옮김"
CTX = ["CG", "CHG", "CHH"]

CX = pd.read_csv(NEW + "/tables/ATH_READUNIT_CONTEXT.tsv.gz", sep="\t")
CXM = pd.read_csv(NEW + "/tables/ATH_CONTEXT_PERARRAY.tsv", sep="\t")
CXM = CXM[CXM.part == "whole"]
G = pd.read_csv(NEW + "/tables/UNITBETA_PERARRAY_v13.tsv.gz", sep="\t",
                usecols=["array_id", "묶음"]).drop_duplicates()
GRP = dict(zip(G.array_id, G.묶음))
ORD = figorder.panel_order(NEW)

A_IDS = [a for a in CXM.array_id.unique()]
A_IDS.sort(key=lambda a: ORD.get(GRP.get(a, ""), (999, 0.0)))

rows = []
for aid in A_IDS:
    for ctx in CTX:
        d = CX[(CX.array_id == aid) & (CX.context == ctx)]
        d = d[d.n_calls > 0]
        if not len(d):
            continue
        rows.append(dict(array_id=aid, 문맥=ctx, read유닛=len(d),
                         자리수중앙=float(np.median(d.n_calls.values)),
                         평균β=round(float(d.n_5mC.sum() / d.n_calls.sum()), 4)))
T = pd.DataFrame(rows)

TMP = os.environ.get("TEMP", ".") + "/x14"
os.makedirs(TMP, exist_ok=True)
T.to_csv(TMP + "/FIGDATA_SUPP4_ATH_CONTEXT.tsv", sep="\t", index=False, encoding="utf8",
         lineterminator="\n")
arch.put(TMP + "/FIGDATA_SUPP4_ATH_CONTEXT.tsv", "tables/FIGDATA_SUPP4_ATH_CONTEXT.tsv", REASON)

mu = T.groupby("문맥").평균β.mean()
print("arrays %d, rows %d" % (T.array_id.nunique(), len(T)))
print("array mean β (call-weighted): %s"
      % " · ".join("%s %.3f" % (k, mu[k]) for k in CTX if k in mu))
print("   These three numbers appear in the manuscript text and the Supp 7-1C legend.")
