# -*- coding: utf-8 -*-
# 하플로타입 중복 array 에서 read 를 배정할 수 있는가
#
#   금화조 chr2 는 mat/pat array 서열이 바이트 단위로 같다.
#   → array 안쪽만 덮은 read 는 두 곳에 똑같이 맞아 원리적으로 배정 불가(mapq 0).
#   → 그러나 flank 까지 문 read 는 flank 서열이 다르므로 배정이 가능해야 한다.
#   이 스크립트는 read 를 두 무리로 갈라 mapq 분포를 비교해 그 예측을 검사한다.
#     A 무리 : array 안쪽만 (flank 물림 < ANCHOR)
#     B 무리 : 한쪽 이상 flank 를 ANCHOR bp 이상 물었다
#   또 mat/pat 에 붙은 read 수가 반반인지(= 정렬기가 임의로 가른 흔적인지) 센다.
#   인자 : bam chrom array시작 array끝 라벨 [ANCHOR=1000]
# MIT License · Copyright (c) 2026 Hungmin Kim
import sys
import collections
import pysam

BAM, CHROM = sys.argv[1], sys.argv[2]
A0, A1 = int(sys.argv[3]), int(sys.argv[4])
LAB = sys.argv[5]
ANCHOR = int(sys.argv[6]) if len(sys.argv) > 6 else 1000
b = pysam.AlignmentFile(BAM, "rb")
grp = collections.defaultdict(lambda: collections.Counter())
mq = collections.defaultdict(list)
for r in b.fetch(CHROM, max(0, A0 - 1 - 50000), A1 + 50000):
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue
    s, e = r.reference_start + 1, r.reference_end
    if e < A0 or s > A1:
        continue
    left = max(0, A0 - s)      # array 왼쪽 바깥으로 문 길이
    right = max(0, e - A1)
    k = "B_flank앵커" if (left >= ANCHOR or right >= ANCHOR) else "A_array내부"
    grp[k]["n"] += 1
    grp[k]["mq0"] += 1 if r.mapping_quality == 0 else 0
    grp[k]["mq60"] += 1 if r.mapping_quality >= 60 else 0
    grp[k]["mq1_59"] += 1 if 0 < r.mapping_quality < 60 else 0
    mq[k].append(r.mapping_quality)
b.close()
print("%-26s %-12s %6s %7s %8s %8s %9s" % (LAB, "무리", "read", "mapq0", "mapq1-59", "mapq>=60", "mapq중앙"))
for k in sorted(grp):
    g = grp[k]
    v = sorted(mq[k])
    print("%-26s %-12s %6d %6.1f%% %7.1f%% %7.1f%% %9d"
          % ("", k, g["n"], 100.0 * g["mq0"] / g["n"], 100.0 * g["mq1_59"] / g["n"],
             100.0 * g["mq60"] / g["n"], v[len(v) // 2]))
