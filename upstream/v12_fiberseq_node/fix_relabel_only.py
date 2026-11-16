# -*- coding: utf-8 -*-
# Undo the relabelling, by the same rule as build_ao_final.py.
#   Why: ao_g8cut moves a read that appears to disagree with the stored placement onto the stored index, as relabeled_d.
# The number is the real position; forcing the stored value is wrong (measured over 66 individuals, 48 against 19).
#        The move is undone as in build_ao_final.py:
#            r[7] = r[7] - d   (array_pos)
#            r[8] = r[8] - d, on the G6 index.
#   Also applied: the coordinate-frame correction for the two repaired slots, by the same rule as fix_shift2.
#   The original is not overwritten.
#   Arguments: the source unit table from ao_g8cut, AO_G8_READS.tsv, and the output.
import io
import re
import sys

NL = chr(10)
TAB = chr(9)
INP, READF, OUTP = sys.argv[1], sys.argv[2], sys.argv[3]
SHIFT = {}

RE = re.compile(r"relabeled_d=([+-]\d+)")
D = {}
for L in io.open(READF, encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) < 11:
        continue
    m = RE.search(r[10])
    if m:
        D[(r[0], r[1])] = int(m.group(1))
print("relabel 된 read %d" % len(D))

n = rel = sh = 0
O = io.open(OUTP, "w", encoding="utf8", newline=NL)
for L in io.open(INP, encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    n += 1
    if len(r) >= 15:
        d = D.get((r[0], r[2]))
        if d and r[7].lstrip("-").isdigit() and r[8].lstrip("-").isdigit():
            r[7] = str(int(r[7]) - d)
            r[8] = str(int(r[8]) - d)
            rel += 1
        k = (r[0], r[3])
        if k in SHIFT and r[7].lstrip("-").isdigit() and int(r[7]) > SHIFT[k]:
            r[7] = str(int(r[7]) - 1)
            sh += 1
    O.write(TAB.join(r) + NL)
O.close()
print("%s -> %s · 전체 %d · relabel 되돌린 유닛 %d · 수리틀 보정 유닛 %d" % (INP, OUTP, n, rel, sh))
