# -*- coding: utf-8 -*-
# The two measurements differ : `haplotype effect` 0.508 (54,472 units) against 0.495 here (778,798 read-units)
#   ⇒ 조리법을 그쪽과 똑같이 맞춰 본다 : 칸 0.005 · 이동평균 5 칸 · 탐색 0.20-0.80
#   ⇒ 그리고 집계 단위(read-유닛 ↔ 유닛)를 바꿔 가며 같은 조리법으로 잰다
#      (기억 : "집계 범위를 바꿔 보기 전에 출처불명 금지")
import collections
import gzip
import os

import numpy as np
import pandas as pd

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    이 글이 읽는 자료의 뿌리는 README.md 의
    "Running this code" 절에 적혀 있다."""
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


W = os.path.join(_need("RDNA_WORK_ROOT"))
RU = os.path.join(W, "out_v11_fix4", "RU.tsv.gz")
L1 = _need("RDNA_HUMAN_L1_PARQUET")

h2 = pd.read_parquet(L1, columns=["sample", "slot", "gate"])
H2 = set(map(tuple, h2.loc[h2.gate == 1, ["sample", "slot"]].drop_duplicates().astype(str).values))

best = {}
for i, L in enumerate(gzip.open(RU, "rt", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    if f[10] != "HPRC_ONT_submissions" or int(f[7]) < 20:
        continue
    if (f[0], f[1]) not in H2:
        continue
    k = (f[0], f[1], f[2], f[3])
    if k not in best or (f[6] == "jx" and best[k][0] != "jx"):
        best[k] = (f[6], float(f[9]), f[0], f[1], f[3])

readunit = np.array([v[1] for v in best.values()])
# 유닛 수준 : 같은 (sample, slot, 유닛번호) 의 read 들을 하나로
g = collections.defaultdict(list)
for v in best.values():
    g[(v[2], v[3], v[4])].append(v[1])
unit_mean = np.array([float(np.mean(x)) for x in g.values()])
unit_med = np.array([float(np.median(x)) for x in g.values()])
print("read-유닛 %s · 유닛(자리) %s" % (format(len(readunit), ","), format(len(unit_mean), ",")))
print()


def valley(b, bw=0.005, sm=5, lo=0.20, hi=0.80):
    n = int(round(1.0 / bw))
    hist, edges = np.histogram(b, bins=n, range=(0.0, 1.0))
    dens = hist / hist.sum()
    mid = (edges[:-1] + edges[1:]) / 2.0
    d = np.convolve(dens, np.ones(sm) / sm, mode="same") if sm > 1 else dens
    m = (mid >= lo) & (mid <= hi)
    j = int(np.argmin(d[m]))
    v = float(mid[m][j])
    dmin = float(d[m][j])
    peak = float(d.max())
    flat = mid[m][d[m] <= dmin * 1.10]
    return v, dmin, peak, (float(flat.min()), float(flat.max()))


print("그쪽 조리법 그대로 : 칸 0.005 · 이동평균 5 · 탐색 0.20-0.80")
print("%-18s %10s %10s %12s %s" % ("집계 단위", "골", "0.508 차", "골/최고", "평평한 구간"))
for nm, b in (("read-유닛", readunit), ("유닛 평균", unit_mean), ("유닛 중앙", unit_med)):
    v, dmin, peak, fl = valley(b)
    print("%-18s %10.4f %10.4f %12.4f %.3f - %.3f"
          % (nm, v, v - 0.508, dmin / peak, fl[0], fl[1]))
print()
print("평활 칸수를 바꿔 보면 (칸 0.005 · read-유닛)")
for sm in (1, 3, 5, 7, 9, 11):
    v, dmin, peak, fl = valley(readunit, sm=sm)
    print("   평활 %2d 칸 : 골 %.4f (0.508 차 %+.4f) · 골/최고 %.4f" % (sm, v, v - 0.508, dmin / peak))
print()
print("평활 칸수를 바꿔 보면 (칸 0.005 · 유닛 평균)")
for sm in (1, 3, 5, 7, 9, 11):
    v, dmin, peak, fl = valley(unit_mean, sm=sm)
    print("   평활 %2d 칸 : 골 %.4f (0.508 차 %+.4f) · 골/최고 %.4f" % (sm, v, v - 0.508, dmin / peak))
print()
print("컷을 옮기면 꺼진 몫이 얼마나 (유닛 평균 기준)")
for cut in (0.45, 0.495, 0.508, 0.538):
    print("   컷 %.3f : read-유닛 %.2f%% · 유닛 평균 %.2f%%"
          % (cut, 100 * (readunit <= cut).mean(), 100 * (unit_mean <= cut).mean()))
