# -*- coding: utf-8 -*-
# 재절단 번호 검산 — 전 슬롯에서 read 유닛 번호가 지도와 어긋난 곳을 찾는다 (2026-09-13)
#   방법 : 슬롯마다 read 유닛 변이 토큰을 지도의 array_pos + k (k=-3..3) 와 비교(자카드 평균).
#          제대로면 k=0 이 최고여야 한다. 다른 k 가 더 높으면 그만큼 번호가 밀린 것이다.
#   조건 : read 유닛 20 개 이상인 슬롯만 · 최고 k 가 0 이 아니고 차이가 0.02 이상일 때 보고
#   인자 : 유닛표 · 지도 · 출력   (예: aog8v10/AO_G8_UNITS.tsv  g8canonv10/...  SHIFT_v10.tsv)
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
