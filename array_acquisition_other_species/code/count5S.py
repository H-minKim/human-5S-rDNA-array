# -*- coding: utf-8 -*-
# Stream the first part of a public RNA-seq run and count the maternal-type (chr4) and somatic-type (chr18) 5S transcripts.
#   The two types differ at 9 positions, 22 to 53. Only the k-mers unique to one type (k = 15) are kept and matched against the read on both strands.
#   A read carrying unique k-mers of both types is counted separately as ambiguous.
#   Arguments: the maximum number of reads, then run:label:url [run:label:url ...]
# MIT License. Copyright (c) 2026 Hyungmin Kim
import gzip
import io
import sys
import urllib.request

K = 15
seqs = {}
name = None
for L in io.open(__file__.replace("count5S.py", "zf5S_types.fa")):
    L = L.strip()
    if L.startswith(">"):
        name = L[1:]
    elif L:
        seqs[name] = L
RC = str.maketrans("ACGT", "TGCA")


def kmers(s):
    out = set()
    for t in (s, s.translate(RC)[::-1]):
        for i in range(len(t) - K + 1):
            out.add(t[i:i + K])
    return out


km = {n: kmers(s) for n, s in seqs.items()}
shared = km["maternal_chr4"] & km["somatic_chr18"]
M = km["maternal_chr4"] - shared
S = km["somatic_chr18"] - shared
ANY = km["maternal_chr4"] | km["somatic_chr18"]
NMAX = int(sys.argv[1])
print("sample\treads\t5S_any\tmaternal\tsomatic\tambiguous\tmaternal_frac")
for arg in sys.argv[2:]:
    run, lab, url = arg.split("|")
    n = anyc = m = s = amb = 0
    try:
        resp = urllib.request.urlopen("https://" + url, timeout=120)
        f = gzip.GzipFile(fileobj=resp)
        for i, line in enumerate(f):
            if i % 4 != 1:
                continue
            r = line.decode().strip()
            n += 1
            hit = False
            hm = hs = False
            for j in range(0, len(r) - K + 1):
                w = r[j:j + K]
                if w in ANY:
                    hit = True
                    if w in M:
                        hm = True
                    elif w in S:
                        hs = True
            if hit:
                anyc += 1
                if hm and hs:
                    amb += 1
                elif hm:
                    m += 1
                elif hs:
                    s += 1
            if n >= NMAX:
                break
        resp.close()
    except Exception as e:
        print("%s\terror %s" % (lab, e))
        continue
    fr = "%.3f" % (m / float(m + s)) if (m + s) else "-"
    print("%s\t%d\t%d\t%d\t%d\t%d\t%s" % (lab, n, anyc, m, s, amb, fr))
    sys.stdout.flush()
