# -*- coding: utf-8 -*-
# bismark BAM 에서 구역별 CpG 메틸을 낸다 (array 와 그 flank 비교용)
#
#   bismark 는 MM/ML 이 아니라 XM 태그에 메틸 호출을 담는다.
#     Z = 메틸된 CpG · z = 안 된 CpG · X/x = CHG · H/h = CHH · U/u = 미지 · . = C 아님
#   그래서 ONT 용 region_meth.py 와 별도로 둔다.
#
#   mapq 로 거르지 않는다(탠덤 array 는 mapq 0 이 정상).
#     대신 mapq>=10 만 쓴 값도 같이 적어 편향을 볼 수 있게 한다.
#   중복(PCR) 표시된 read 는 뺀다.
#
#   인자 : bam 라벨 구역목록.tsv 출력.tsv
#     구역목록 : 이름 <tab> chrom <tab> start <tab> end
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import pysam

TAB = chr(9)
NL = chr(10)
BAM, LAB, REG, OUT = sys.argv[1:5]
b = pysam.AlignmentFile(BAM, "rb")
REFS = set(b.references)
new = not os.path.exists(OUT)
o = io.open(OUT, "a", encoding="utf8", newline=NL)
if new:
    o.write(TAB.join(["sample", "region", "chrom", "start", "end", "reads",
                      "n_cpg", "n_meth", "frac_CpG", "n_cpg_mq10", "frac_CpG_mq10",
                      "n_chg", "frac_CHG", "n_chh", "frac_CHH"]) + NL)
for L in io.open(REG, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) < 4:
        continue
    nm, ch, s, e = f[0], f[1], int(f[2]), int(f[3])
    if ch not in REFS:
        o.write(TAB.join([LAB, nm, ch, str(s), str(e)] + ["-"] * 10) + NL)
        continue
    nr = 0
    cg = me = 0
    cgq = meq = 0
    chg = chgm = 0
    chh = chhm = 0
    for r in b.fetch(ch, max(0, s - 1), e):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
            continue
        try:
            xm = r.get_tag("XM")
        except KeyError:
            continue
        nr += 1
        q = r.mapping_quality
        for qp, rp in r.get_aligned_pairs(matches_only=True):
            if rp is None or qp is None or not (s <= rp + 1 <= e):
                continue
            c = xm[qp] if qp < len(xm) else "."
            if c in ("z", "Z"):
                cg += 1
                if c == "Z":
                    me += 1
                if q >= 10:
                    cgq += 1
                    if c == "Z":
                        meq += 1
            elif c in ("x", "X"):
                chg += 1
                if c == "X":
                    chgm += 1
            elif c in ("h", "H"):
                chh += 1
                if c == "H":
                    chhm += 1
    o.write(TAB.join([LAB, nm, ch, str(s), str(e), str(nr), str(cg), str(me),
                      "%.4f" % (me / float(cg)) if cg else "-",
                      str(cgq), "%.4f" % (meq / float(cgq)) if cgq else "-",
                      str(chg), "%.4f" % (chgm / float(chg)) if chg else "-",
                      str(chh), "%.4f" % (chhm / float(chh)) if chh else "-"]) + NL)
o.close()
b.close()
print("%s 완료" % LAB)
