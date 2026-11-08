# -*- coding: utf-8 -*-
# array 를 코딩부 / 스페이서 / 바깥(flank) 로 나눠 CpG 메틸을 낸다.
#
#   구간이 잘게 쪼개져 있으므로(chr18 코딩부 65 조각) 구간마다 fetch 하지 않고
#     array+flank 를 한 번만 훑으면서 CpG 자리마다 어느 부류인지 판정한다.
#
#   구간은 read 묶음 부트스트랩으로 낸다 — 같은 read 안의 CpG 는 독립이 아니다.
#     read 를 단위로 복원추출해 비율을 다시 계산하고 2.5/97.5 백분위를 쓴다.
#     코딩부는 CpG 가 몇 개뿐일 수 있어 이항 구간을 쓰면 크게 과소평가된다.
#
#   인자 : bam 라벨 구역조각.tsv 출력.tsv [부트스트랩 횟수(기본 2000)]
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import random
import sys
import collections
import pysam

TAB = chr(9)
NL = chr(10)
BAM, LAB, REG, OUT = sys.argv[1:5]
NB = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
random.seed(0)

# 구역 조각 → 부류별 구간 목록
parts = collections.defaultdict(list)          # (base, kind) -> [(chrom,s,e)]
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    base, _, kind = f[0].rpartition("__")
    parts[(base, kind)].append((f[1], int(f[2]), int(f[3])))

b = pysam.AlignmentFile(BAM, "rb")
REFS = set(b.references)
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline=NL)
if new:
    o.write(TAB.join(["sample", "locus", "part", "n_read", "n_cpg", "n_meth",
                      "frac", "lo95", "hi95", "bp"]) + NL)


def boot(per_read):
    """read 를 단위로 복원추출 — 같은 read 안 CpG 는 같이 움직인다"""
    n = len(per_read)
    if n < 2:
        return ("-", "-")
    out = []
    for _ in range(NB):
        c = m = 0
        for _ in range(n):
            a, d = per_read[random.randrange(n)]
            c += a
            m += d
        if c:
            out.append(m / float(c))
    if not out:
        return ("-", "-")
    out.sort()
    return ("%.4f" % out[int(0.025 * len(out))], "%.4f" % out[int(0.975 * len(out))])


for (base, kind) in sorted(parts):
    iv = parts[(base, kind)]
    ch = iv[0][0]
    if ch not in REFS:
        o.write(TAB.join([LAB, base, kind] + ["-"] * 7) + NL)
        continue
    iv_s = sorted((s, e) for _, s, e in iv)
    bp = sum(e - s + 1 for s, e in iv_s)
    # 빠른 판정을 위해 자리 집합 대신 정렬된 구간에 이분 탐색
    import bisect
    starts = [s for s, e in iv_s]
    ends = [e for s, e in iv_s]

    def inside(p):
        i = bisect.bisect_right(starts, p) - 1
        return i >= 0 and p <= ends[i]

    lo, hi = iv_s[0][0], iv_s[-1][1]
    per_read = []
    tot_c = tot_m = 0
    for r in b.fetch(ch, max(0, lo - 1), hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        try:
            xm = r.get_tag("XM")
        except KeyError:
            continue
        c = m = 0
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp is None or qp is None:
                continue
            if not inside(rp + 1):
                continue
            ch2 = xm[qp] if qp < len(xm) else "."
            if ch2 in ("z", "Z"):
                c += 1
                if ch2 == "Z":
                    m += 1
        if c:
            per_read.append((c, m))
            tot_c += c
            tot_m += m
    fr = "%.4f" % (tot_m / float(tot_c)) if tot_c else "-"
    l, h = boot(per_read)
    o.write(TAB.join([LAB, base, kind, str(len(per_read)), str(tot_c), str(tot_m),
                      fr, l, h, str(bp)]) + NL)
o.close()
b.close()
print("%s 완료" % LAB)
