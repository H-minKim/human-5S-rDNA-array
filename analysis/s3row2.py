# -*- coding: utf-8 -*-
"""표S3 2 번 줄을 양분으로 (2026-10-14).

옛 값 : 중메틸 0.582 / 고메틸 0.570  ← 삼분위 시절 · 정의 불명 (재현 못 함)
새로 낼 때는 정의를 표에 적는다.

재는 것 : TFIIIB 구역 창(-28..-10) 안에서,
          그 구간에 뉴클레오솜이 없는 분자의 m6A/A
대안     : "footprint 는 사실 작은 뉴클레오솜이다"
예측     : 그렇다면 뉴클레오솜 없는 분자의 m6A/A 도 낮아야 한다
"""
import os
import numpy as np

DD = r"$RDNA_FHMM_DIR\data"
z = np.load(os.path.join(DD, "nucfree_flankcore_span.npz"))
nA, nM, nmol = z["nA"], z["nM"], z["nmol"]
rel = np.arange(int(z["lo"]), int(z["hi"]) + 1)
SITE = [("TFIIIB", 0, (-28, -10)), ("ICR", 1, (57, 68))]

print("## 구역 창 안 m6A/A — 뉴클레오솜 없는 분자 ↔ 있는 분자")
print("%-8s %-6s %-12s %10s %10s %9s"
      % ("자리", "층", "무리", "A", "m6A", "m6A/A"))
OUT = {}
for nm, j, (a, b) in SITE:
    m = (rel >= a) & (rel <= b)
    for si, sl in ((0, "저메틸"), (1, "고메틸")):
        for fr, fn in ((1, "뉴클 없음"), (0, "뉴클 있음")):
            A, M = nA[j, si, fr][m].sum(), nM[j, si, fr][m].sum()
            r = M / A if A else float("nan")
            if fr == 1:
                OUT[(nm, sl)] = (r, int(A), int(nmol[j, si, 1]))
            print("%-8s %-6s %-12s %10s %10s %9.4f"
                  % (nm, sl, fn, format(int(A), ","), format(int(M), ","), r))
    print()

print("## 표S3 2 번 줄에 넣을 값 (TFIIIB · 양분)")
for sl in ("저메틸", "고메틸"):
    r, A, n = OUT[("TFIIIB", sl)]
    print("   %-6s %.4f   (A %s · 분자 %s)" % (sl, r, format(A, ","), format(n, ",")))
print()
print("   판정 : 뉴클레오솜 있는 분자보다 훨씬 높으면 '작은 뉴클레오솜' 설명은 기각")
