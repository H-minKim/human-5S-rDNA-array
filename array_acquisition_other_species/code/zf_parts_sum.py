# -*- coding: utf-8 -*-
# 부류별 표를 만든다(코딩부 / 스페이서 / 바깥) — 결과를 파일로 쓴다(콘솔 인코딩 회피)
import io
import collections
import sys

TAB = chr(9)
NL = chr(10)
src = sys.argv[1]
dst = sys.argv[2]
LOC = sys.argv[3] if len(sys.argv) > 3 else "chr18_somatic"

d = collections.defaultdict(dict)
for L in io.open(src, encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if f[0] == "sample":
        continue
    d[(f[0], f[1])][f[2]] = f[3:]

order = ["알_Jiang", "알_Potok", "16세포_1.5hpf", "1k세포_3hpf_MBT", "sphere_4hpf",
         "Germring_6hpf", "배아24hpf", "정자_Jiang", "정자_Potok", "정소",
         "근육성체_Potok", "장성체", "심장성체"]

o = io.open(dst, "w", encoding="utf8", newline=NL)
o.write(TAB.join(["sample", "locus", "coding_frac", "coding_lo", "coding_hi",
                  "coding_read", "coding_cpg", "spacer_frac", "spacer_lo", "spacer_hi",
                  "spacer_read", "spacer_cpg", "flank_frac", "flank_read", "flank_cpg",
                  "coding_minus_spacer", "coding_minus_flank", "spacer_minus_flank"]) + NL)

loci = sorted(set(k[1] for k in d)) if LOC == "ALL" else [LOC]
for lo_ in loci:
    for s in order:
        k = (s, lo_)
        if k not in d:
            continue
        c = d[k].get("코딩부")
        p = d[k].get("스페이서")
        fc = fm = fr = 0
        for q in ("바깥_flankL", "바깥_flankR"):
            v = d[k].get(q)
            if v and v[1] != "-":
                fr += int(v[0])
                fc += int(v[1])
                fm += int(v[2])
        ff = fm / float(fc) if fc else None

        def g(v, i):
            return v[i] if v and v[i] != "-" else "-"
        cf = float(c[3]) if c and c[3] != "-" else None
        pf = float(p[3]) if p and p[3] != "-" else None
        row = [s, lo_,
               g(c, 3), g(c, 4), g(c, 5), g(c, 0), g(c, 1),
               g(p, 3), g(p, 4), g(p, 5), g(p, 0), g(p, 1),
               ("%.4f" % ff) if ff is not None else "-", str(fr), str(fc),
               ("%+.4f" % (cf - pf)) if (cf is not None and pf is not None) else "-",
               ("%+.4f" % (cf - ff)) if (cf is not None and ff is not None) else "-",
               ("%+.4f" % (pf - ff)) if (pf is not None and ff is not None) else "-"]
        o.write(TAB.join(row) + NL)
o.close()
print("wrote %s" % dst)
