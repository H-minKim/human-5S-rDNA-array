# -*- coding: utf-8 -*-
# 애기장대 세 문맥(CG·CHG·CHH)의 array 별 요약 — `tables/FIGDATA_SUPP4_ATH_CONTEXT.tsv`
#
#   왜 이 글이 생겼나
#     이 표의 세 수(CG 0.844 · CHG 0.423 · CHH 0.146)가 원고 본문과 legend 에 실린다.
#     그런데 지금까지 그 표를 쓰는 유일한 코드가 그림을 그리는 글(`xb2_supp.py`)이었다.
#     논문에 들어가는 수는 결과를 만드는 코드 가 내야 한다 — 그리는 코드는 논문에 넣지 않는다.
#     ⇒ 계산을 그대로 이리로 옮겼다. 값은 한 자리도 바뀌지 않는다(바이트 대조로 확인한다).
#
#   자(ruler) : `평균β` 는 호출 가중 이다 — array 안 모든 read×유닛에 대해 `Σn_5mC / Σn_calls`.
#     read 유닛 β 의 산술평균(`ATH_CONTEXT_PERARRAY` 의 CG 0.841)과 다른 수 다.
#     집계 단위가 다르므로 섞어 쓰면 틀린다 — 원고가 쓰는 값은 이 글이 내는 쪽(0.844) 이다.
#
#   차례 : 염색체 번호 → array 시작좌표(`figorder.panel_order`). 그림 배치가 아니라 자료 차례 이고,
#     `PERARRAY_META`·`CHROM_NAMES`·`UNITBETA_PERARRAY_v13` 세 결과표에서만 나온다.
#
#   입력  : tables/ATH_READUNIT_CONTEXT.tsv.gz · tables/ATH_CONTEXT_PERARRAY.tsv (둘 다 `x13_ath_context.py`)
#           tables/UNITBETA_PERARRAY_v13.tsv.gz (`x9D`) · tables/PERARRAY_META.tsv (`x34`) · tables/CHROM_NAMES.tsv (`x9Q`)
#   산출  : tables/FIGDATA_SUPP4_ATH_CONTEXT.tsv
import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
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
CXM = CXM[CXM.part == "전체"]
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
print("array %d 개 · 줄 %d" % (T.array_id.nunique(), len(T)))
print("array 평균 β (호출가중) : %s"
      % " · ".join("%s %.3f" % (k, mu[k]) for k in CTX if k in mu))
print("   이 세 수가 원고 본문과 Supp 7-1C legend 에 실린다.")
