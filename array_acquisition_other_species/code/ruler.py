# -*- coding: utf-8 -*-
# 공개 전장 반복단위를 자로 써서 실측 유닛을 잰다
#   자 : X12811 (인간 5S 반복단위 2,231 bp · 전장) · AY130738 (애기장대 495 bp) · AF213516 (제브라피시 175 bp)
#   재는 대상 : 같은 계통 실측 대표 유닛
#   n_HSP 1 + 덮힘 ~1.0 = 유닛 하나 / HSP 여러개 규칙적 = 여러 유닛 덩어리 / 덮힘 낮음 = 자 밖 서열이 큼
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import glob
import subprocess

NL = chr(10)
TAB = chr(9)


def rdfa(p):
    d, k = {}, None
    for L in io.open(p, encoding="utf8"):
        L = L.strip()
        if not L:
            continue
        if L[0] == ">":
            k = L[1:].split()[0]
            d[k] = []
        else:
            d[k].append(L.upper())
    return dict((a, "".join(b)) for a, b in d.items())


U = {}
for p in sorted(glob.glob("units/*.unit.fa")):
    U[os.path.basename(p).replace(".unit.fa", "")] = list(rdfa(p).values())[0]
ALL = rdfa("ALL_lit.fa")
RULERS = [("X12811.1", ("human", "rhesus", "chimp", "bonobo", "gorilla", "orang", "siamang")),
          ("AY130738.1", ("arabidopsis",)),
          ("AF213516.1", ("zebrafish",)),
          ("D43630.1", ("mouse",)),
          ("AF419701.3", ("zfinch",))]
for racc, tags in RULERS:
    if racc not in ALL:
        print("자 없음 " + racc)
        continue
    io.open("ruler.fa", "w", encoding="utf8", newline=NL).write(">" + racc + NL + ALL[racc] + NL)
    print(NL + "=== 자 : %s (%d bp)  -> %s" % (racc, len(ALL[racc]), "/".join(tags)))
    print("%-52s %6s %4s %6s %7s %s" % ("array_id", "unit", "nHSP", "cov", "pid_med", "HSP 위치(유닛 좌표)"))
    keep = [a for a in U if any(t in a for t in tags)]
    for a in sorted(keep, key=lambda z: (z.split("__")[0], -len(U[z]))):
        io.open("one.fa", "w", encoding="utf8", newline=NL).write(">" + a + NL + U[a] + NL)
        r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "11", "-dust", "no", "-evalue", "1e-5",
                            "-query", "ruler.fa", "-subject", "one.fa",
                            "-outfmt", "6 pident length sstart send bitscore"], capture_output=True, text=True)
        hs = []
        for L in r.stdout.splitlines():
            f = L.split(TAB)
            hs.append((float(f[0]), int(f[1]), int(f[2]), int(f[3]), float(f[4])))
        hs.sort(key=lambda x: -x[-1])
        taken, pids = [], []
        for pid, aln, ss, se, bs in hs:
            x, y = min(ss, se), max(ss, se)
            if any(not (y < p or x > q) for p, q in taken):
                continue
            taken.append((x, y))
            pids.append(pid)
        if not taken:
            print("%-52s %6d %4s %6s %7s NOHIT" % (a, len(U[a]), "0", "0.000", "-"))
            continue
        cov = sum(y - x + 1 for x, y in taken)
        pids.sort()
        sp = sorted(taken)
        print("%-52s %6d %4d %6.3f %7.1f %s" % (a, len(U[a]), len(taken), cov / float(len(U[a])),
                                                pids[len(pids) // 2],
                                                ",".join("%d-%d" % (x, y) for x, y in sp[:5]) +
                                                (" ... 총 %d" % len(sp) if len(sp) > 5 else "")))
