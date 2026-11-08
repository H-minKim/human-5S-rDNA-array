# -*- coding: utf-8 -*-
# junction-only map: coordinates of only the 7 operated slots are returned to the old canonical coordinates (2026-09-13)
#
#   Why: junction reads are not re-aligned. The stored alignments are on the old canonical.
#        But the coordinates of the new map are on the new canonical, so they disagree in the 7 operated slots
#        (measured: 72.5% of HG02559_2 junction positions went OK -> NO_FULL).
#   How: relative to the operation window [a0, a1]
#        - unit before the window (unit_end < a0)       -> coordinates unchanged
#        - unit after the window (unit_start > a1_new)  -> subtract delta from the coordinates (new->old)
#        - unit overlapping the window                  -> removed from the map (the old canonical has no matching sequence)
#        The number (array_pos) keeps the new numbering. Only coordinates are moved to the old frame.
#        delta = splice_win_bp - asm_win_bp  (new window length - old window length)
#   The original map is not touched. A new file is written.
#   Arguments: new map, CANON_V10_PATCH.tsv, output map
import io
import sys

NL = chr(10)
TAB = chr(9)
MAPF, PATCHF, OUTF = sys.argv[1], sys.argv[2], sys.argv[3]

P = {}
f = io.open(PATCHF, encoding="utf8")
h = f.readline().rstrip(NL).split(TAB)
ix = dict((n, i) for i, n in enumerate(h))
for L in f:
    r = L.rstrip(NL).split(TAB)
    a0 = int(r[ix["anchor_start"]])
    oldw = int(r[ix["asm_win_bp"]])
    neww = int(r[ix["splice_win_bp"]])
    delta = neww - oldw                       # new - old
    P[(r[ix["sample"]], r[ix["slot"]])] = (a0, a0 + neww - 1, delta)   # window on the new canonical

fm = io.open(MAPF, encoding="utf8")
head = fm.readline()
cols = head.rstrip(NL).split(TAB)
ci = dict((n, i) for i, n in enumerate(cols))
O = io.open(OUTF, "w", encoding="utf8", newline=NL)
O.write(head if head.endswith(NL) else head + NL)
kept = dropped = moved = 0
for L in fm:
    r = L.rstrip(NL).split(TAB)
    k = (r[0], r[1])
    if k not in P:
        O.write(TAB.join(r) + NL)
        kept += 1
        continue
    a0, a1n, delta = P[k]
    try:
        s = int(r[ci["unit_start"]])
        e = int(r[ci["unit_end"]])
    except ValueError:
        O.write(TAB.join(r) + NL)
        kept += 1
        continue
    if e < a0:                                # before the window: unchanged
        O.write(TAB.join(r) + NL)
        kept += 1
    elif s > a1n:                             # after the window: to old coordinates
        r[ci["unit_start"]] = str(s - delta)
        r[ci["unit_end"]] = str(e - delta)
        O.write(TAB.join(r) + NL)
        moved += 1
    else:                                     # overlaps the window: no match on the old canonical
        dropped += 1
O.close()
print("%s -> %s" % (MAPF, OUTF))
print("  그대로 %d · 좌표 되돌림 %d · 뺀 행(수술 창과 겹침) %d" % (kept, moved, dropped))
