# -*- coding: utf-8 -*-
# read 단위 메틸 — 같은 read 안의 CpG 는 독립이 아니다.
#   CpG 수로 낸 이항 구간은 표본크기를 read 수가 아니라 CpG 수로 세어
#   불확실성을 크게 과소평가한다(실측: Potok 알 chr18 array 를 덮은 read 가 17 개뿐인데
#   CpG 82 개로 ±0.107 이 나왔다).
#   그래서 read 마다 CpG 메틸 비율을 하나 내고, read 를 표본단위로 삼아 평균과 구간을 낸다.
#
#   인자 : bam 라벨 구역목록.tsv 출력.tsv  [최소 read당 CpG 수(기본 3)]
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import math
import os
import sys
import pysam

TAB = chr(9)
NL = chr(10)
BAM, LAB, REG, OUT = sys.argv[1:5]
MINC = int(sys.argv[5]) if len(sys.argv) > 5 else 3

b = pysam.AlignmentFile(BAM, "rb")
REFS = set(b.references)
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline=NL)
if new:
    o.write(TAB.join(["sample", "region", "n_read", "read_mean",
                      "read_sd", "read_se", "n_cpg_total"]) + NL)

for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    nm, ch, s, e = f[0], f[1], int(f[2]), int(f[3])
    if ch not in REFS:
        o.write(TAB.join([LAB, nm] + ["-"] * 5) + NL)
        continue
    vals = []
    tot = 0
    for r in b.fetch(ch, max(0, s - 1), e):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        try:
            xm = r.get_tag("XM")
        except KeyError:
            continue
        cg = me = 0
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp is None or qp is None or not (s <= rp + 1 <= e):
                continue
            c = xm[qp] if qp < len(xm) else "."
            if c in ("z", "Z"):
                cg += 1
                if c == "Z":
                    me += 1
        if cg >= MINC:
            vals.append(me / float(cg))
            tot += cg
    n = len(vals)
    if n == 0:
        o.write(TAB.join([LAB, nm, "0"] + ["-"] * 4) + NL)
        continue
    mu = sum(vals) / n
    if n > 1:
        sd = math.sqrt(sum((v - mu) ** 2 for v in vals) / (n - 1))
        se = sd / math.sqrt(n)
    else:
        sd = se = float("nan")
    o.write(TAB.join([LAB, nm, str(n), "%.4f" % mu,
                      "%.4f" % sd, "%.4f" % se, str(tot)]) + NL)
o.close()
b.close()
print("%s 완료" % LAB)
