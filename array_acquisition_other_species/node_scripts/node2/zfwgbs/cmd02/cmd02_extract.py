# -*- coding: utf-8 -*-
# ONT 배치 BAM(MM/ML) → read × array × 부류(코딩/스페이서) CpG 호출 표
#   CpG 전용 모델(금화조 dorado 5mCG_5hmCG · 제브라 guppy 5mC_CG)만 — 문맥 판정 없이 5mC('m') 호출을 CpG 로 본다.
#   메틸 판정 : p(5mC) = ML/255 ≥ 0.5 (배치 경로 modkit 기준과 같은 문턱)
#   부류 : CpG 의 C 위치가 cmsearch 유전자 구간 안이면 코딩, array 안 나머지는 스페이서
#   보조/이차 정렬 제외 · mapq 로 거르지 않음(탠덤 안은 mapq0 정상) — mapq 는 칸으로 남김
#   인자 : bam 시료 REGIONS.tsv 출력.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import bisect
import io
import os
import sys
import pysam

T = "\t"
BAM, LAB, REG, OUT = sys.argv[1:5]
regs = []
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip("\n").split(T)
    if f[0] == "array_id":
        continue
    gi = []
    for x in f[7].split(","):
        if x:
            a, b = x.split("-")
            gi.append((int(a), int(b)))
    regs.append((f[0], f[1], int(f[2]), int(f[3]), gi))
b = pysam.AlignmentFile(BAM, "rb")
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline="\n")
if new:
    o.write(T.join(["sample", "array_id", "read_id", "mapq", "part", "n_calls", "n_meth"]) + "\n")
tot = 0
for aid, ch, a0, a1, gi in regs:
    gs = [x for x, _ in gi]
    ge = [y for _, y in gi]

    def coding(p):
        i = bisect.bisect_right(gs, p) - 1
        return i >= 0 and p <= ge[i]

    for r in b.fetch(ch, a0 - 1, a1):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        mb = r.modified_bases
        if not mb:
            continue
        q2r = dict(r.get_aligned_pairs(matches_only=True))
        cnt = {"코딩부": [0, 0], "스페이서": [0, 0]}
        for key, calls in mb.items():
            if key[2] != "m":
                continue
            rev = key[1] == 1
            for qp, ql in calls:
                rp = q2r.get(qp)
                if rp is None:
                    continue
                cpos = rp if not rev else rp - 1      # CpG 의 C (1 기반은 +1)
                p1 = cpos + 1
                if p1 < a0 or p1 > a1:
                    continue
                part = "코딩부" if coding(p1) else "스페이서"
                cnt[part][0] += 1
                if ql / 255.0 >= 0.5:
                    cnt[part][1] += 1
        for part, (c, m) in cnt.items():
            if c:
                o.write(T.join([LAB, aid, r.query_name, str(r.mapping_quality), part, str(c), str(m)]) + "\n")
                tot += 1
o.close()
b.close()
print("%s 행 %d" % (LAB, tot))
