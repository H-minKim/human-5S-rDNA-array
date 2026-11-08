# -*- coding: utf-8 -*-
# 군 BAM 마다 read 하나당 (좌위 갈래, 부류, CpG 수, 메틸 수) 를 뽑는다 — 부트스트랩·군간검정은 로컬에서
#   좌위 갈래 : chr4(6 array 묶음) / chr18
#   부류 : 코딩부 / 스페이서 / 바깥
#   인자 : 구역조각.tsv 출력.tsv bam:군 [bam:군 ...]
# MIT License · Copyright (c) 2026 Hungmin Kim
import bisect
import collections
import io
import sys
import pysam

TAB = chr(9)
NL = chr(10)
REG, OUT = sys.argv[1], sys.argv[2]
parts = collections.defaultdict(list)   # (갈래, 부류) -> [(chrom,s,e)]
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    base, _, kind = f[0].rpartition("__")
    br = "chr18" if base.startswith("chr18") else "chr4"
    k = kind
    parts[(br, k)].append((f[1], int(f[2]), int(f[3])))

o = io.open(OUT, "w", encoding="utf8", newline=NL)
o.write(TAB.join(["group", "branch", "part", "read", "n_cpg", "n_meth"]) + NL)
for arg in sys.argv[3:]:
    bam, grp = arg.split(":")
    b = pysam.AlignmentFile(bam, "rb")
    for (br, k), iv in sorted(parts.items()):
        bych = collections.defaultdict(list)
        for ch, s, e in iv:
            bych[ch].append((s, e))
        for ch, ivs in bych.items():
            ivs.sort()
            st = [s for s, e in ivs]
            en = [e for s, e in ivs]
            # 구간이 흩어져 있으므로 구간마다 fetch 하되 read 는 한 번만 센다
            seen = {}
            for s, e in ivs:
                for r in b.fetch(ch, s - 1, e):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                        continue
                    key = (r.query_name, r.reference_start)
                    if key in seen:
                        continue
                    try:
                        xm = r.get_tag("XM")
                    except KeyError:
                        continue
                    c = m = 0
                    for qp, rp in r.get_aligned_pairs(matches_only=True):
                        if rp is None or qp is None:
                            continue
                        p = rp + 1
                        i = bisect.bisect_right(st, p) - 1
                        if i < 0 or p > en[i]:
                            continue
                        z = xm[qp] if qp < len(xm) else "."
                        if z in ("z", "Z"):
                            c += 1
                            m += (z == "Z")
                    seen[key] = (c, m)
            for (qn, rs), (c, m) in seen.items():
                if c:
                    o.write(TAB.join([grp, br, k, "%s:%d" % (qn, rs), str(c), str(m)]) + NL)
    b.close()
o.close()
print("완료")
