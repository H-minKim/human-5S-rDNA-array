# -*- coding: utf-8 -*-
# A both-strand conversion check for single-end data: each read is called as a C-to-T strand, low in C, or a G-to-A strand, low in G.
#   PBAT produces reads from both strands mixed together, so the two signatures are taken together.
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
    print("%s\tread %d\tC %.1f%%\tG %.1f%%\tT %.1f%%\tA %.1f%%\tC->T type %.1f%%\tG->A type %.1f%%\tsum %.1f%%" % (
        path.split("/")[-1], n, 100.0 * C / tot, 100.0 * G / tot, 100.0 * T / tot, 100.0 * A / tot,
        100.0 * ct / n, 100.0 * ga / n, 100.0 * (ct + ga) / n))
