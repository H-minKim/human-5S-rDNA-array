# -*- coding: utf-8 -*-
# Add the missing unit ticks to the hand-made master SVG.
#
# What is missing: on the second and third purple reads on the left, the ticks stop at x 86.4.
#   The bars run to 104.9 and 119.5, leaving 18.5 pt and 33.1 pt bare.
# Cause: after the ticks were drawn, the label overwrite ran without a graphics argument and
#   erased the artwork its box touched, that is those ticks. The log records 21 and 26 ticks drawn
#   while only 17 survived in the output.
# The coordinates are taken from the corrected PDF rather than assumed; the spacing is 3.08 pt.
# The form of the master file is followed exactly: <path class="st53" d="M{x},{y}v2.9"/>, to one decimal place.
import io, os, re, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8")
# -- (repo) paths come from the environment --------------------------------
#@DROP
SRC = os.environ.get("USER_SVG", "USER_Fig1.svg")
OUT = os.environ.get("PATCHED_SVG", "Fig1_tickfix.svg")

NEW = {
    "229.5": [89.5, 92.6, 95.7, 98.7, 101.8],                       # the bar ends at 104.9
    "231.8": [89.5, 92.6, 95.7, 98.7, 101.8, 104.9, 108.0, 111.1, 114.2, 117.2],  # the bar ends at 119.5
}

s0 = io.open(SRC, encoding="utf8", errors="replace").read()
s = s0
added = 0
for y, xs in NEW.items():
    anchor = '<path class="st53" d="M86.4,%sv2.9"/>' % y
    assert s.count(anchor) == 1, (y, s.count(anchor))
    ins = anchor + "".join('<path class="st53" d="M%s,%sv2.9"/>' % (("%.1f" % x), y) for x in xs)
    s = s.replace(anchor, ins)
    added += len(xs)
    print("   y %s : %d 개 넣음 → x %s" % (y, len(xs), xs))

assert len(s) > len(s0)
io.open(OUT, "w", encoding="utf8", newline="\n").write(s)
print()
print("넣은 눈금 : %d 개" % added)
print("%s  %s -> %s" % (OUT, hashlib.md5(s0.encode("utf8")).hexdigest()[:12],
                         hashlib.md5(s.encode("utf8")).hexdigest()[:12]))

# Check: do the ticks reach the end of the bar on every row?
import collections
BAR = {"229.5": 104.9, "231.8": 119.5, "227.0": 87.3, "234.2": 133.2, "236.6": 123.4}
rows = collections.defaultdict(list)
for m in re.finditer(r'<path class="st53" d="M([0-9.]+),([0-9.]+)v[-0-9.]+"/>', s):
    rows[m.group(2)].append(float(m.group(1)))
print()
print("검산 : 줄마다 마지막 눈금 → 막대 끝")
for y in sorted(BAR, key=float):
    xs = sorted(rows.get(y, []))
    gap = BAR[y] - xs[-1] if xs else 999
    print("   y %-6s 눈금 %2d  마지막 %7.2f  막대끝 %7.2f  빈 구간 %5.2f pt  %s"
          % (y, len(xs), xs[-1], BAR[y], gap, "민둥" if gap > 5 else "ok"))
