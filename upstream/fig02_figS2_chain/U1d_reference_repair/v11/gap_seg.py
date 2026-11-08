# -*- coding: utf-8 -*-
# Step 1: extract only the anchor segments on each node (2026-09-12), the reads are split across nodes and must be gathered
#   anchor = start of the intact unit before the gap to end of the intact unit after the gap
#   Even if the anchor ends are not exactly aligned, the nearest aligned position within +-100bp is used
#   (requiring exactness discards reads in bulk, measured 15 -> 2)
#   Output: <output folder>/seg_<individual>_<slot>_<front unit>.<node>.fa  (header line has read id and actual anchor offset)
import io
import os
import re
import sys
import subprocess
import collections

H = "/home/ubuntu"
NL = chr(10)
TAB = chr(9)
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
NODE = subprocess.run(["hostname"], capture_output=True, text=True).stdout.strip()
SLACK = 100

U = collections.defaultdict(list)
f = io.open(os.path.join(H, "g8canonv8", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] == "UNIT":
        U[(r[0], r[1])].append((int(r[ix["array_pos"]]), int(r[ix["unit_start"]]), int(r[ix["unit_end"]])))
SITES = []
for k, v in U.items():
    v.sort()
    for i in range(len(v) - 1):
        g = v[i + 1][1] - v[i][2] - 1
        if g != 0:
            SITES.append((k[0], k[1], v[i][0], v[i][1], v[i + 1][2], g))
SITES.sort()

tot = collections.Counter()
for (smp, sl, ap, a0, a1, g) in SITES:
    recs = []
    for d, off in ((os.path.join(H, "placed_bam"), 5000), (os.path.join(H, "placed_bam_g8"), 0)):
        if not os.path.isdir(d):
            continue
        for fn in os.listdir(d):
            if not (fn.startswith(smp + ".") and fn.endswith(".bam")):
                continue
            b = os.path.join(d, fn)
            src = fn.split(".")[2]
            reg = "%s_h%s:%d-%d" % (smp, sl, max(1, a0 + off), a1 + off)
            out = subprocess.run(["samtools", "view", b, reg], capture_output=True, text=True).stdout
            for L in out.splitlines():
                r = L.split(TAB)
                if len(r) < 10 or (int(r[1]) & 0x900) or r[9] == "*":
                    continue
                pos = int(r[3]) - off
                q = 0
                rp = pos - 1
                m = {}
                for gg in CIG.finditer(r[5]):
                    k2, op = int(gg.group(1)), gg.group(2)
                    if op == "S":
                        q += k2
                    elif op == "H":
                        pass
                    elif op in "M=X":
                        for t in range(k2):
                            m[rp + t] = q + t
                        rp += k2
                        q += k2
                    elif op == "I":
                        q += k2
                    elif op in "DN":
                        rp += k2
                if not m:
                    continue
                lo, hi = min(m), max(m)
                if lo > a0 - 1 + SLACK or hi < a1 - 1 - SLACK:
                    continue                      # Anchor not covered enough
                p0 = next((p for p in range(a0 - 1, a0 - 1 + SLACK) if p in m), None)
                p1 = next((p for p in range(a1 - 1, a1 - 1 - SLACK, -1) if p in m), None)
                if p0 is None or p1 is None:
                    continue
                q0, q1 = m[p0], m[p1]
                if q1 - q0 < 500:
                    continue
                recs.append((r[0], src, p0 - (a0 - 1), (a1 - 1) - p1, r[9][q0:q1 + 1]))
    if not recs:
        continue
    p = os.path.join(OUT, "seg_%s_%s_%d.%s.fa" % (smp, sl, ap, NODE))
    with io.open(p, "w", newline=NL) as gfh:
        for rid, src, d0, d1, sq in recs:
            gfh.write(">%s|%s|%d|%d%s%s%s" % (rid, src, d0, d1, NL, sq, NL))
    tot[(smp, sl, ap)] = len(recs)
for k, v in sorted(tot.items()):
    print("  %s h%s 유닛%d : 조각 %d" % (k[0], k[1], k[2], v))
print("자리 %d · 조각 합 %d" % (len(tot), sum(tot.values())))
