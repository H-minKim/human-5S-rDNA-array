# -*- coding: utf-8 -*-
# Check the indices after the re-cut: find, across every slot, where the read unit index disagrees with the map.
#   The method: per slot, the read unit's variant tokens are compared with the map's array_pos + k for k from -3 to 3, by mean Jaccard.
#          If all is well k = 0 scores highest; a higher score at another k means the index is shifted by that much.
#   Conditions: slots with at least 20 read units only, reported where the best k is not 0 and the difference is at least 0.02.
#   Arguments: the unit table, the map, the output (for example aog8v10/AO_G8_UNITS.tsv, g8canonv10/..., SHIFT_v10.tsv).
import io
import sys
import socket
import collections

NL = chr(10)
TAB = chr(9)
UNITF, MAPF, OUTF = sys.argv[1], sys.argv[2], sys.argv[3]
MINN = 20
MARGIN = 0.02

f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
M = {}
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] != "UNIT":
        continue
    v = r[ix["varlist"]]
    M[(r[0], r[1], int(r[ix["array_pos"]]))] = set() if v == "-" else set(v.split(","))

R = collections.defaultdict(list)
for L in io.open(UNITF, encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) < 15 or r[13] not in ("OK", "OK_FILLED"):
        continue
    R[(r[0], r[3])].append((int(r[7]), set() if r[14] == "-" else set(r[14].split(","))))

OUT = io.open(OUTF, "w", encoding="utf8", newline=NL)
OUT.write(TAB.join(["sample", "slot", "n_units", "best_k", "j_best", "j_0", "gap"]) + NL)
bad = 0
tot = 0
for k2 in sorted(R):
    rows = R[k2]
    if len(rows) < MINN:
        continue
    tot += 1
    sc = {}
    for k in range(-3, 4):
        v = []
        for pos, a in rows:
            b = M.get((k2[0], k2[1], pos + k))
            if b is None:
                continue
            u = a | b
            v.append(1.0 if not u else 1.0 * len(a & b) / len(u))
        sc[k] = sum(v) / len(v) if v else 0.0
    bk = max(sc, key=lambda x: sc[x])
    if bk != 0 and sc[bk] - sc[0] >= MARGIN:
        bad += 1
        OUT.write(TAB.join([k2[0], k2[1], str(len(rows)), str(bk),
                            "%.3f" % sc[bk], "%.3f" % sc[0], "%.3f" % (sc[bk] - sc[0])]) + NL)
OUT.close()
print("%s  검사 슬롯 %d · 번호 어긋난 슬롯 %d  -> %s" % (socket.gethostname(), tot, bad, OUTF))
