# -*- coding: utf-8 -*-
# PBAT 단일말단 read 전처리 — cutadapt 대신(노드에 없다)
#   ① 5' 9 bp 자름 : PBAT 는 무작위 프라이머로 합성해 앞쪽이 참조와 안 맞는다
#   ② 3' 어댑터 자름 : Illumina 공통 `AGATCGGAAGAGC`(어댑터는 바이설파이트 뒤에 붙어 변환 안 됨)
#      전체 일치가 없으면 read 끝에 걸친 앞부분(≥5 bp)도 잘라 낸다
#   ③ 30 bp 미만 버림
#   인자 : in.fastq.gz out.fastq.gz
# MIT License · Copyright (c) 2026 Hungmin Kim
import gzip
import sys

ADP = "AGATCGGAAGAGC"
CLIP5 = 9
MINL = 30
src, dst = sys.argv[1], sys.argv[2]
n = kept = adp = 0
with gzip.open(src, "rt") as f, gzip.open(dst, "wt") as o:
    while True:
        h = f.readline()
        if not h:
            break
        s = f.readline().rstrip("\n")
        f.readline()
        q = f.readline().rstrip("\n")
        n += 1
        i = s.find(ADP)
        if i < 0:
            for k in range(len(ADP) - 1, 4, -1):
                if s.endswith(ADP[:k]):
                    i = len(s) - k
                    break
        if i >= 0:
            adp += 1
            s, q = s[:i], q[:i]
        s, q = s[CLIP5:], q[CLIP5:]
        if len(s) < MINL:
            continue
        kept += 1
        o.write(h + s + "\n+\n" + q + "\n")
sys.stderr.write("%s read %d · 어댑터 %d(%.1f%%) · 남김 %d(%.1f%%)\n" % (
    src.split("/")[-1], n, adp, 100.0 * adp / max(n, 1), kept, 100.0 * kept / max(n, 1)))
