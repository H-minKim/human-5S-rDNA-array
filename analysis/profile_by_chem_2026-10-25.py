# -*- coding: utf-8 -*-
#  사용자 물음 (2026-10-25) : `Fig 4 a-c` 가 전부 H2(화학 섞임)인데 어떻게 하나.
#  이미 잰 것 : 밀도(b)의 두 봉우리·골·컷 띠 — 층마다 같았다.
#  아직 안 잰 것 : `a` 는 **위치별 프로파일**이다. 모양 자체를 층별로 재야 한다.
#  자 : `L1_READ_UNIT.parquet` 의 하드콜 β(`b_call`)를 `ap`(근위 1 부터) 띠로 묶어
#        층별 평균·저메틸 몫을 견준다. gate==1(H2)만.
# The HPRC release 2 ONT sequencing index this script reads is not part of this archive.
# Set RDNA_HPRC_R2_ONT_INDEX to a local copy; without it the script falls back to a file of
# that name beside itself. The README says where the index comes from and gives its md5.
import collections
import csv
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

hp = {}
for r in csv.DictReader(io.open(os.environ.get("RDNA_HPRC_R2_ONT_INDEX",
                                     os.path.join(SP, "hprc_ont_r2_index_2026-10-24.csv")),
                                encoding="utf8", errors="replace")):
    c = (r.get("sequencing_chemistry") or "").strip()
    if c:
        hp.setdefault(r["sample_id"], set()).add(c)
CH = {k: (sorted(v)[0] if len(v) == 1 else "혼재") for k, v in hp.items()}

cols = ["sample", "slot", "gate", "ap", "ad", "N", "b_call"]
D = pd.read_parquet(os.path.join(V11, "02_tables", "L1_READ_UNIT.parquet"), columns=cols)
D = D[D.gate == 1].copy()
D["층"] = D["sample"].map(lambda s: CH.get(str(s), "모름(색인밖)"))
print("H2 read-유닛 %s · 개체 %d" % (format(len(D), ","), D["sample"].nunique()))


def band(ap, ad):
    if ap <= 3:
        return "근위 1-3"
    if ap <= 10:
        return "근위 4-10"
    if ad <= 3:
        return "원위 1-3"
    if ad <= 10:
        return "원위 4-10"
    return "내부"


D["띠"] = [band(a, d) for a, d in zip(D.ap.values, D.ad.values)]
ORD = ["근위 1-3", "근위 4-10", "내부", "원위 4-10", "원위 1-3"]
print()
print("## `Fig 4 a` 의 위치별 프로파일 — 화학 층별 (하드콜 β 평균 · 저메틸 몫 β<0.508)")
print()
hdr = "| 띠 |" + "".join(" %s 평균 | %s 저메틸몫 |" % (k, k)
                       for k in ["전체", "R941", "R1041", "모름"])
print(hdr)
print("|---" * (1 + 8) + "|")
rows = []
for b in ORD:
    g = D[D["띠"] == b]
    cells = []
    for k in ["전체", "R941", "R1041", "모름(색인밖)"]:
        s = g if k == "전체" else g[g["층"] == k]
        if len(s) < 100:
            cells += ["-", "-"]
            continue
        cells += ["%.4f" % s.b_call.mean(), "%.3f" % (s.b_call < 0.508).mean()]
    print("| %s | %s |" % (b, " | ".join(cells)))
    rows.append([b] + cells)

print()
print("## 층별 근위-원위 차 (프로파일의 기울기 — 이것이 패널 a 의 주장이다)")
print()
print("| 층 | 개체 | read-유닛 | 근위 1-3 평균 | 원위 1-3 평균 | 차 | 근위 저메틸몫 | 원위 저메틸몫 |")
print("|---|--:|--:|--:|--:|--:|--:|--:|")
for k in ["전체", "R941", "R1041", "모름(색인밖)"]:
    s = D if k == "전체" else D[D["층"] == k]
    p = s[s["띠"] == "근위 1-3"]
    q = s[s["띠"] == "원위 1-3"]
    if len(p) < 100 or len(q) < 100:
        continue
    print("| %s | %d | %s | %.4f | %.4f | **%+.4f** | %.3f | %.3f |"
          % (k, s["sample"].nunique(), format(len(s), ","), p.b_call.mean(),
             q.b_call.mean(), p.b_call.mean() - q.b_call.mean(),
             (p.b_call < 0.508).mean(), (q.b_call < 0.508).mean()))
P = os.path.join(SP, "profile_by_chem_2026-10-25.tsv")
out = [TAB.join(["band", "all_mean", "all_lowfrac", "R941_mean", "R941_lowfrac",
                 "R1041_mean", "R1041_lowfrac", "unknown_mean", "unknown_lowfrac"])]
for r in rows:
    out.append(TAB.join(str(x) for x in r))
io.open(P, "w", encoding="utf8", newline=NL).write(NL.join(out) + NL)
print()
print("기록 : %s" % P)
