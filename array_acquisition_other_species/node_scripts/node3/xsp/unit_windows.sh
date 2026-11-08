#!/bin/bash
# 유닛 내부 구성 : 유닛을 300 bp 창으로 잘라 그 종 게놈 전체에 몇 번 나오는지 센다
#   창이 게놈 곳곳에 수백~수천 번 나오면 = 흩어진 반복요소 (사람이면 AluY), array 안에서만 나오면 = array 고유 구간
#   대표 유닛은 array_struct.py 가 뽑은 units.fa 의 가운데 유닛 (틀 : 1-119 가 5S 유전자, 그 뒤가 스페이서)
#   인자 : 종키 · (선택) array 이름 접두사 (기본 arr0)
#   산출 : ~/xsp/struct/win_<species>.tsv  (창시작 · 창끝 · 정렬수 · array밖 정렬수)
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
SP="${1:?종키}"
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
print("유닛 %d bp · 창 %d 개" % (len(u), n))
PYEOF
minimap2 -x sr -N 2000 -p 0.6 -t 8 "$FAS" "/tmp/win_$SP.fa" 2> /dev/null \
  | awk -F'\t' -v OFS='\t' '{split($1,a,"_"); c[$1]++} END{for(k in c) print k, c[k]}' | sort -t_ -k2,2n > "$X/struct/win_$SP.tsv"
echo "== $SP 창별 게놈 내 정렬 수 (창 = 유닛 안 좌표 · 1-119 가 5S 유전자)"
cat "$X/struct/win_$SP.tsv"
