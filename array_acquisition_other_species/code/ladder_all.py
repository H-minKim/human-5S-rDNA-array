# -*- coding: utf-8 -*-
# 77 array 전부 : 자기-도트플롯 옵셋 사다리로 반복 주기를 직접 판정한다
#   진짜 탠덤 = 옵셋 히스토그램이 p · 2p · 3p ... 로 서는 사다리 + 옵셋 p 의 정렬량이 조각 길이에 육박
#   (사람 실측 : 옵셋 2240 에 147,862 bp / 조각 150,000 · HSP 66 개
#    붉은털원숭이 : 최상위 옵셋도 7,714 bp 뿐 · HSP 2,898 개 = 사다리 없음)
#   유닛·cmsearch·공개 서열 전부 안 쓰고 array 서열만 본다
#   산출 : LADDER_ALL.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import glob
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser("~/lit5s")
BIN = 20
MINALN = 200
MINP = 130
w = io.open(os.path.join(H, "LADDER_ALL.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "slice_len", "n_hsp", "ladder_period", "aln_at_p", "aln_frac_of_slice",
                  "n_rungs", "rung_offsets", "verdict", "top_offsets"]) + NL)
print("%-50s %7s %6s %7s %8s %5s %5s  %s" % ("array_id", "slice", "HSP", "주기p", "p정렬량", "비율", "단수", "판정"))


def ladder(off, slen):
    if not off:
        return 0, 0, 0, []
    # 후보 : 정렬량 상위 옵셋 가운데 가장 작은 것부터
    cand = [k for k, v in off.most_common(12) if k >= MINP]
    best = (0, 0, 0, [])
    for p in sorted(cand):
        tol = max(BIN, int(0.03 * p))
        rungs, tot = [], 0
        for k in range(1, 13):
            d = p * k
            if d > slen * 0.8:
                break
            hit = [(o, off[o]) for o in off if abs(o - d) <= tol]
            if not hit:
                break
            o, v = max(hit, key=lambda x: x[1])
            rungs.append(o)
            tot += v
        if len(rungs) >= 3 and rungs[0] and off[rungs[0]] > best[1]:
            best = (p, off[rungs[0]], len(rungs), rungs)
    if best[0]:
        return best
    p, v = off.most_common(1)[0]
    return p, v, 1, [p]


for sfa in sorted(glob.glob(os.path.join(H, "slices", "*.arr.fa"))):
    aid = os.path.basename(sfa).replace(".arr.fa", "")
    hd = io.open(sfa, encoding="utf8").readline()
    kv = dict(x.split("=", 1) for x in hd.split() if "=" in x)
    slen = int(kv["slice_len"])
    r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "16", "-dust", "no", "-evalue", "1e-20",
                        "-query", sfa, "-subject", sfa, "-outfmt", "6 qstart qend sstart send length"],
                       capture_output=True, text=True)
    off = collections.Counter()
    n = 0
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        qs, qe, ss, se, aln = int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4])
        if se < ss or aln < MINALN:
            continue
        d = ss - qs
        if d <= 0:
            continue
        n += 1
        off[int(round(d / float(BIN))) * BIN] += aln
    p, ap, nr, rungs = ladder(off, slen)
    frac = ap / float(slen) if slen else 0.0
    if p == 0:
        v = "반복 없음"
    elif nr >= 3 and frac >= 0.50:
        v = "탠덤 확인 (주기 %d · 사다리 %d 단)" % (p, nr)
    elif nr >= 3 and frac >= 0.15:
        v = "탠덤이나 일부만 (주기 %d · 비율 %.2f)" % (p, frac)
    elif nr >= 3:
        v = "약한 사다리 (주기 %d · 비율 %.2f)" % (p, frac)
    else:
        v = "사다리 없음 = 일정한 유닛 주기 없음"
    print("%-50s %7d %6d %7d %8d %5.2f %5d  %s" % (aid, slen, n, p, ap, frac, nr, v))
    w.write(TAB.join([aid, str(slen), str(n), str(p), str(ap), "%.3f" % frac, str(nr),
                      ",".join(str(x) for x in rungs), v,
                      ",".join("%d:%d" % (k, vv) for k, vv in off.most_common(6))]) + NL)
    sys.stdout.flush()
w.close()
print("표 : ~/lit5s/LADDER_ALL.tsv")
