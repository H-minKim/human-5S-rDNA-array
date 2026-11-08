# -*- coding: utf-8 -*-
# v11-fix1: G6 map coordinate correction for ARRAY_ONLY (edition 2, 2026-09-17)
#
#   Defect of edition 1 (mkg6map_fix1.py): G6 units overlapping the surgery window were removed entirely (HG00133 12, HG03270 9).
#     For reads whose stored number range ended inside the window, the window (stored G6 range +-250bp) got shorter and slipped by 5-6 units in the tandem
#     (HG03270 read 5288a70a: v11 matched the fingerprint perfectly but edition 1 was -5, 5s-ao-window-slide-trap).
#   Edition 2: remove only the units that were actually lost. Lost units = read unit DP events (old map numbers)
#         HG00133 h2: old G8 units 122-125 (+ gap piece between 124/125)   HG03270 h1: old G8 unit 43
#     For each G6 unit, check which unit (or gap) of the old G8 map (g8canonv10) its center coordinate falls in
#       - before the first lost unit            -> coordinate unchanged
#       - after the last lost unit              -> coordinate + delta (old -> new)
#       - inside a lost unit or a gap between   -> removed
#     delta = splice_win_bp - asm_win_bp.  Check: start of old G8 (last lost unit +1) + delta = start of new map (first lost number)
#   Arguments: original G6 map, old G8 map, new G8 map, repair log table, output, lost units ("sample:slot:first:last,...")
import io, sys, bisect, collections
NL = chr(10); TAB = chr(9)
G6F, OLD8, NEW8, PATCHF, OUTF, DELS = sys.argv[1:7]
DEL = {}
for x in DELS.split(","):
    s, sl, a, b = x.split(":"); DEL[(s, sl)] = (int(a), int(b))
P = {}
f = io.open(PATCHF, encoding="utf8"); h = f.readline().rstrip(NL).split(TAB); ix = dict((n, i) for i, n in enumerate(h))
for L in f:
    r = L.rstrip(NL).split(TAB)
    P[(r[ix["sample"]], r[ix["slot"]])] = int(r[ix["splice_win_bp"]]) - int(r[ix["asm_win_bp"]])
def units(p):
    U = collections.defaultdict(dict)
    f = io.open(p, encoding="utf8"); ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    for L in f:
        r = L.rstrip(NL).split(TAB)
        if (r[0], r[1]) in DEL and r[ix["seg_type"]] == "UNIT":
            U[(r[0], r[1])][int(r[ix["array_pos"]])] = (int(r[ix["unit_start"]]), int(r[ix["unit_end"]]))
    return U
U8o, U8n = units(OLD8), units(NEW8)
for k, (a, b) in DEL.items():
    chk = U8o[k][b + 1][0] + P[k]; new_first = U8n[k][a][0]
    print("  검산 %s h%s : 옛 유닛 %d 시작 %d + delta %d = %d · 새 지도 유닛 %d 시작 %d · 차 %d bp" % (k[0], k[1], b + 1, U8o[k][b + 1][0], P[k], chk, a, new_first, chk - new_first))
fm = io.open(G6F, encoding="utf8"); head = fm.readline()
ci = dict((n, i) for i, n in enumerate(head.rstrip(NL).split(TAB)))
O = io.open(OUTF, "w", encoding="utf8", newline=NL); O.write(head if head.endswith(NL) else head + NL)
stat = collections.defaultdict(collections.Counter)
for L in fm:
    r = L.rstrip(NL).split(TAB); k = (r[0], r[1])
    if k not in DEL:
        O.write(TAB.join(r) + NL); continue
    try:
        s = int(r[ci["unit_start"]]); e = int(r[ci["unit_end"]])
    except ValueError:
        O.write(TAB.join(r) + NL); stat[k]["숫자아님"] += 1; continue
    a, b = DEL[k]; d = P[k]
    lo = U8o[k][a][0]; hi = U8o[k][b][1]          # lost span (old coordinates, including gaps in between)
    mid = (s + e) // 2
    if mid < lo:
        O.write(TAB.join(r) + NL); stat[k]["앞 그대로"] += 1
    elif mid > hi:
        r[ci["unit_start"]] = str(s + d); r[ci["unit_end"]] = str(e + d)
        O.write(TAB.join(r) + NL); stat[k]["뒤 옮김"] += 1; stat[k]["_last"] = int(r[ci["unit_end"]])
    else:
        stat[k]["뺌(없어진 유닛)"] += 1
O.close()
for k in sorted(stat):
    c = stat[k]
    print("  %s h%s : 앞 %d · 옮김 %d · 뺌 %d · 마지막 행 끝 %d" % (k[0], k[1], c["앞 그대로"], c["뒤 옮김"], c["뺌(없어진 유닛)"], c["_last"]))
