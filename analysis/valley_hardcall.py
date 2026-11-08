# -*- coding: utf-8 -*-
# 하드콜 β 분포의 골을 다시 잰다
#   바탕은 `fig40_claim_A_B.py` 패널 a 와 똑같다 :
#     RU.tsv.gz · HPRC_ONT_submissions · n_cpg>=20 · H2 슬롯 · ao/jx 중복제거(jx 우선) · 하드콜 β
#   다른 것은 칸 폭뿐 : 0.02 (그림) → 0.005 · 0.002 · 0.001
#   읽기만 한다. 아무것도 쓰지 않는다(표 하나만 모래상자에).
import collections
import gzip
import io
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
print("H2 하플로타입(슬롯) %d" % len(H2))

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
        best[k] = (f[6], float(f[9]))
b = np.array([v[1] for v in best.values()])
print("read-유닛 %s" % format(len(b), ","))
print()


def valley(bw, lo=0.20, hi=0.80, smooth=0):
    n = int(round(1.0 / bw))
    hist, edges = np.histogram(b, bins=n, range=(0.0, 1.0))
    dens = hist / hist.sum()
    mid = (edges[:-1] + edges[1:]) / 2.0
    d = dens
    if smooth:                      # 이동평균 (칸 수 smooth)
        k = np.ones(smooth) / smooth
        d = np.convolve(dens, k, mode="same")
    m = (mid >= lo) & (mid <= hi)
    j = np.argmin(d[m])
    v = mid[m][j]
    dmin = d[m][j]
    flat = mid[m][d[m] <= dmin * 1.10]
    return v, dmin, (flat.min(), flat.max()), hist, mid, dens


print("%-10s %-10s %-12s %-22s %s" % ("칸 폭", "골 위치", "그 칸 몫", "골 평평한 구간(+10%)", "0.508 과 차"))
rows = []
for bw, sm in ((0.02, 0), (0.005, 0), (0.005, 3), (0.002, 0), (0.002, 5), (0.001, 0), (0.001, 11)):
    v, dmin, fl, hist, mid, dens = valley(bw, smooth=sm)
    tag = "%.3f%s" % (bw, (" sm%d" % sm) if sm else "")
    rows.append((tag, v, dmin, fl))
    print("%-10s %-10.4f %-12.5f %-22s %+.4f"
          % (tag, v, dmin, "%.3f - %.3f" % fl, v - 0.508))
print()

# 0.508 이 골 바닥에서 얼마나 떨어져 있나 (칸 0.002, 평활 5)
v, dmin, fl, hist, mid, dens = valley(0.002, smooth=5)
k = np.ones(5) / 5
d = np.convolve(dens, k, mode="same")
i508 = int(0.508 / 0.002)
print("칸 0.002(평활 5) 기준 : 골 바닥 %.4f (몫 %.5f) · 0.508 자리의 몫 %.5f · 비 %.3f"
      % (v, dmin, d[i508], d[i508] / dmin if dmin else float("nan")))
print()
# 컷을 0.508 로 둘 때와 실측 골로 둘 때 꺼진 몫이 얼마나 달라지나
for cut in (0.508, round(v, 4)):
    print("   컷 %.4f → 꺼진 read-유닛 몫 %.4f%%" % (cut, 100.0 * (b <= cut).mean()))
print()
# 부트스트랩 : 골 위치가 얼마나 흔들리나
rng = np.random.default_rng(20261007)
vs = []
for _ in range(200):
    s = rng.choice(b, size=len(b), replace=True)
    hist, edges = np.histogram(s, bins=500, range=(0.0, 1.0))
    dens = hist / hist.sum()
    mid = (edges[:-1] + edges[1:]) / 2.0
    d = np.convolve(dens, np.ones(5) / 5, mode="same")
    m = (mid >= 0.20) & (mid <= 0.80)
    vs.append(mid[m][np.argmin(d[m])])
vs = np.array(vs)
print("부트스트랩 200 회 (칸 0.002 · 평활 5) : 골 중앙 %.4f · 2.5-97.5%% %.4f - %.4f"
      % (np.median(vs), np.percentile(vs, 2.5), np.percentile(vs, 97.5)))

SP = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(SP, "VALLEY_HARDCALL_H2.tsv")
with io.open(P, "w", encoding="utf8", newline=chr(10)) as fh:
    fh.write(chr(9).join(["bin", "valley", "dens_min", "flat_lo", "flat_hi", "diff_vs_0.508"]) + chr(10))
    for tag, v2, dm, fl2 in rows:
        fh.write(chr(9).join([tag, "%.4f" % v2, "%.6f" % dm, "%.3f" % fl2[0],
                              "%.3f" % fl2[1], "%+.4f" % (v2 - 0.508)]) + chr(10))
print()
print("적었다 :", P)
