# -*- coding: utf-8 -*-
"""
unitdef_v11 - unit definition helpers (built **anew** from the v11 canonical reference)

Why rebuild
  The old `fig_common.py` reads `canon.fa` from the old analysis folder.
  The version has changed, so pulling files from the old folder would contaminate results.
  Here the canonical sequence comes from **v11_io.CAN** (= the central 2,240bp of v11 `03_<map>/g8canonv10/canon_ctx.fa`).

What was carried over (definitions from fig_common.py, verbatim)
  ELEMENTS, element(), iscpg(), _blocks(), parse(), vclass(), zone(), ZONES
  The logic is unchanged down to the character. Only the source of CAN changed.
  (load_map/two_stage/figure helpers of fig_common are not used in v11, so they were not carried over.)
"""
import re
import v11_io as _G

CAN = _G.CAN
LCAN = len(CAN)
assert LCAN == 2240, LCAN
UNIT_BP = 2240
STRATA = [(40, 80, "N40-80"), (81, 110, "N81-110"), (111, 400, "N111+")]
ZONES = ["P1_3", "P4_8", "P9_15", "P16_24", "MID", "D16_24", "D9_15", "D4_8", "D1_3"]
ELEMENTS = [(1, 318, "AluY"), (319, 1178, "NTS1"), (1179, 1551, "GArich"),
            (1552, 2121, "NTS2"), (2122, 2240, "coding")]

def element(p):
    if p is None: return ""
    for a, b, n in ELEMENTS:
        if a <= p <= b: return n
    return ""

def iscpg(p):
    if p is None or p < 1 or p > LCAN: return False
    c = CAN[p - 1]
    return (c == "C" and p < LCAN and CAN[p] == "G") or (c == "G" and p >= 2 and CAN[p - 2] == "C")

def _blocks(seq):
    iv = []
    for m in re.finditer(r"(.)\1{3,}", seq): iv.append((m.start() + 1, m.end()))
    for k in (2, 3, 4, 5, 6):
        for m in re.finditer(r"((.{%d}))\1{2,}" % k, seq): iv.append((m.start() + 1, m.end()))
    iv.sort(); mg = []
    for a, b in iv:
        if mg and a <= mg[-1][1] + 1: mg[-1][1] = max(mg[-1][1], b)
        else: mg.append([a, b])
    return [tuple(x) for x in mg]

BLK = _blocks(CAN)
BID = {}
for i, (a, b) in enumerate(BLK, 1):
    for p in range(a, b + 1): BID[p] = i
# Must equal the repeat blocks counted separately by the gate (same rule)
assert BLK == _G.BLOCKS, "반복블록 정의가 관문과 다르다"

def parse(x):
    """variant string ->(class, pos, ref, alt, len, block, in_repeat)"""
    if x.startswith("DB"):
        b, n = x[2:].split(":"); return ("del", None, "", "", int(n), int(b), True)
    if x.startswith("IB"):
        b, n = x[2:].split(":"); return ("ins", None, "", "", int(n), int(b), True)
    if x[0] == "D":
        p, n = x[1:].split(":"); p = int(p); return ("del", p, CAN[p - 1], "", int(n), BID.get(p), False)
    if x[0] == "I":
        q = x[1:].split(":"); p = int(q[0]); return ("ins", p, "", "", int(q[1]), BID.get(p), False)
    if x[0] == "R":
        p = int(x[1:-2]); return ("snv", p, x[-2], x[-1], 1, BID.get(p), True)
    p = int(x[:-2]); return ("snv", p, x[-2], x[-1], 1, BID.get(p), False)

def vclass(x):
    cl, pos, r, a, ln, blk, rep = parse(x)
    return cl + (" in repeat" if (rep or blk) else "")

def zone(p, N):
    d = N - p + 1
    if p <= 3: return "P1_3"
    if p <= 8: return "P4_8"
    if p <= 15: return "P9_15"
    if p <= 24: return "P16_24"
    if d <= 3: return "D1_3"
    if d <= 8: return "D4_8"
    if d <= 15: return "D9_15"
    if d <= 24: return "D16_24"
    return "MID"
