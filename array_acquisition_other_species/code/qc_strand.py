# -*- coding: utf-8 -*-
# 양가닥 변환 검사(단일말단용) — C→T 가닥(C 적음) · G→A 가닥(G 적음) 을 read 마다 판정.
#   PBAT 는 read 가 두 가닥에서 섞여 나오므로 두 서명의 합으로 본다.
import gzip
import sys

for path in sys.argv[1:]:
    n = ct = ga = 0
    C = G = T = A = tot = 0
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 != 1:
                continue
            s = line.strip()
            L = len(s)
            if L < 30:
                continue
            c, g, t, a = s.count("C"), s.count("G"), s.count("T"), s.count("A")
            C += c; G += g; T += t; A += a; tot += L
            n += 1
            if 100.0 * c / L <= 5:
                ct += 1
            elif 100.0 * g / L <= 5:
                ga += 1
            if n >= 200000:
                break
    print("%s\tread %d\tC %.1f%%\tG %.1f%%\tT %.1f%%\tA %.1f%%\tC→T형 %.1f%%\tG→A형 %.1f%%\t합 %.1f%%" % (
        path.split("/")[-1], n, 100.0 * C / tot, 100.0 * G / tot, 100.0 * T / tot, 100.0 * A / tot,
        100.0 * ct / n, 100.0 * ga / n, 100.0 * (ct + ga) / n))
