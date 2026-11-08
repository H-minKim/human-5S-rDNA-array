# -*- coding: utf-8 -*-
# 시료×좌위갈래×부류(코딩부/스페이서/flankL/flankR) 메틸 — read 묶음 부트스트랩(독립, 정렬 전 차이)
import io
import collections
import numpy as np

rng = np.random.default_rng(0)
NB = 2000
TAB = chr(9)
d = collections.defaultdict(list)
for L in io.open("ADULT_PERREAD.tsv", encoding="utf8"):
    f = L.rstrip("\n").split(TAB)
    if f[0] == "group":
        continue
    d[(f[0], f[1], f[2])].append((int(f[4]), int(f[5])))


def boot(v):
    a = np.array(v, dtype=float)
    n = len(a)
    idx = rng.integers(0, n, size=(NB, n))
    c = a[idx, 0].sum(1)
    m = a[idx, 1].sum(1)
    return m / c          # 정렬하지 않는다


def pt(v):
    a = np.array(v, dtype=float)
    return a[:, 1].sum() / a[:, 0].sum()


groups = ["알_Jiang", "알_Potok", "정자_Jiang", "정자_Potok", "16세포", "1k세포", "sphere", "Germring",
          "배아24hpf", "F_GSC_35to70", "M_GSC_35to70", "early_GSC_1to28", "F_Control_35to70",
          "M_Control_35to70", "early_Control_1to28", "근육_Potok", "장", "심장", "정소"]
parts = ["코딩부", "스페이서", "바깥_flankL", "바깥_flankR"]
out = io.open("ADULT_STATS.tsv", "w", encoding="utf8", newline="\n")
out.write(TAB.join(["sample", "branch", "part", "n_read", "n_cpg", "frac", "lo95", "hi95"]) + "\n")
diffo = io.open("ADULT_DIFFS.tsv", "w", encoding="utf8", newline="\n")
diffo.write(TAB.join(["sample", "branch", "contrast", "diff", "lo95", "hi95"]) + "\n")
for br in ("chr4", "chr18"):
    print("\n#### %s   (값 [95%%] r=read)" % br)
    print("시료 | 코딩부 | 스페이서 | flankL | flankR | 코딩−flank(L+R) | 스페이서−flank")
    for g in groups:
        cells = []
        B = {}
        for p in parts:
            v = d.get((g, br, p), [])
            if len(v) >= 3:
                b = boot(v)
                B[p] = b
                lo, hi = np.percentile(b, [2.5, 97.5])
                cells.append("%.3f [%.3f,%.3f] r%d" % (pt(v), lo, hi, len(v)))
                out.write(TAB.join([g, br, p, str(len(v)), str(sum(x[0] for x in v)), "%.4f" % pt(v), "%.4f" % lo, "%.4f" % hi]) + "\n")
            else:
                cells.append("- r%d" % len(v))
        fl = d.get((g, br, "바깥_flankL"), []) + d.get((g, br, "바깥_flankR"), [])
        extra = []
        if len(fl) >= 3:
            bf = boot(fl)
            for p in ("코딩부", "스페이서"):
                if p in B:
                    dd = B[p] - bf
                    lo, hi = np.percentile(dd, [2.5, 97.5])
                    val = pt(d[(g, br, p)]) - pt(fl)
                    extra.append("%+.3f [%+.3f,%+.3f]" % (val, lo, hi))
                    diffo.write(TAB.join([g, br, p + "-flank", "%+.4f" % val, "%+.4f" % lo, "%+.4f" % hi]) + "\n")
                else:
                    extra.append("-")
        else:
            extra = ["-", "-"]
        print("%s | %s | %s" % (g, " | ".join(cells), " | ".join(extra)))
out.close()
diffo.close()
