#!/bin/bash
# The internal composition of the unit: it is cut into 300 bp windows and each window is counted across that species' whole genome.
#   A window occurring hundreds or thousands of times genome-wide is a dispersed repeat element (AluY in human); one occurring only inside the array is specific to the array.
#   The representative unit is an interior unit of units.fa from array_struct.py, in the frame where 1-119 is the 5S gene and the spacer follows.
#   Arguments: species key, and optionally the array name prefix (default arr0).
#   Output: struct/win_<species>.tsv (window start, window end, alignments, alignments outside any array).
# MIT License. Copyright (c) 2026 Hyungmin Kim
set -u
SP="${1:?species key}"
ARR="${2:-arr0}"
X=~/xsp
U=$X/struct/$SP.units.fa
FAS=$X/asm/$SP.fna
W=300
STEP=150
python3 - "$U" "$ARR" "$W" "$STEP" "/tmp/win_$SP.fa" <<'PYEOF'
import io, sys
fa, arr, w, step, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
d, k = {}, None
for L in io.open(fa, encoding="utf8"):
    if L.startswith(">"):
        k = L[1:].strip()
        d[k] = []
    elif k:
        d[k].append(L.strip())
d = dict((a, "".join(b)) for a, b in d.items())
ks = sorted([a for a in d if ("_" + arr + "_") in a])
if not ks:
    sys.exit("no unit for " + arr)
u = d[ks[len(ks) // 2]]
o = io.open(out, "w", encoding="utf8", newline=chr(10))
n = 0
for s in range(0, max(1, len(u) - w + 1), step):
    seq = u[s:s + w]
    if len(seq) < w // 2:
        break
    o.write(">w_%d_%d%s%s%s" % (s + 1, s + len(seq), chr(10), seq, chr(10)))
    n += 1
o.close()
print("unit %d bp, windows %d" % (len(u), n))
PYEOF
minimap2 -x sr -N 2000 -p 0.6 -t 8 "$FAS" "/tmp/win_$SP.fa" 2> /dev/null \
  | awk -F'\t' -v OFS='\t' '{split($1,a,"_"); c[$1]++} END{for(k in c) print k, c[k]}' | sort -t_ -k2,2n > "$X/struct/win_$SP.tsv"
echo "== $SP alignments in the genome per window (window = coordinate within the unit, 1-119 is the 5S gene)"
cat "$X/struct/win_$SP.tsv"
