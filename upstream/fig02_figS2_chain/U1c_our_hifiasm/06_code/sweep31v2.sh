#!/bin/bash
# sweep31 v2: data preparation from 4 passes to 1 pass (2026-09-04)
#
#  Why v1 was slow (measured)
#    `zcat | awk | gzip` ran four times: 1 for counting unit content + 3 ratios.
#    Four single-threaded passes over 460~582 Mb compressed files took 1 hour just for preparation
#    (hifiasm had not even started).
#
#  What v2 changed
#    (1) parallel decompression with unpigz (zcat is single-threaded)
#    (2) one pass writes the three sets f=2,3,4 at the same time and also counts total bases
#    (3) no gzip, hifiasm reads plain fastq. Compression is pure waste.
#  The subsampling method (order c%f) is kept identical to v1, changing it would break comparison with earlier results.
#  As soon as two haplotypes come out, stop and delete the remaining pieces.
set -u
SLOT=${1:-0}; NSLOT=${2:-1}; TH=${TH:-8}
D=~/sweep31; mkdir -p "$D"; SUM=$D/sweep31_summary.tsv
PY=""
for c in ~/miniconda3/bin/python python3; do
  p=$(command -v "$c" 2>/dev/null) || continue
  "$p" -c "import sys" 2>/dev/null && { PY="$p"; break; }
done
[ -n "$PY" ] || { echo "FATAL python 없음"; exit 9; }
UNZ=$(command -v unpigz || command -v zcat)

i=-1
for S in $(tr -d '\r' < ~/sweep31_list.txt); do
  i=$((i+1)); [ $((i % NSLOT)) -ne "$SLOT" ] && continue
  grep -q "^$S	" "$SUM" 2>/dev/null && continue
  Q=""
  for p in ~/s5s_ont/hits/$S.hits.fastq.gz ~/giab_hits/$S.hits.fastq.gz ~/s5s_new/hits/$S.hits.fastq.gz; do
    [ -f "$p" ] && { Q=$p; break; }
  done
  [ -n "$Q" ] || { echo "  $S read 없음"; continue; }
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  [ "${FREE:-0}" -lt 40 ] && { echo "  디스크 ${FREE}G 중단"; break; }

  P="$D/p_${S}_$SLOT"
  rm -rf "$P"; mkdir -p "$P"
  T0=$(date +%s)
  # three ratios + total bases in one pass
  "$UNZ" -c "$Q" | awk -v a="$P/f2.fq" -v b="$P/f3.fq" -v d="$P/f4.fq" -v t="$P/tot" '
    NR%4==1{h=$0} NR%4==2{s=$0; n+=length(s)} NR%4==3{p=$0}
    NR%4==0{c++
      if(c%2==0) printf "%s\n%s\n%s\n%s\n",h,s,p,$0 > a
      if(c%3==0) printf "%s\n%s\n%s\n%s\n",h,s,p,$0 > b
      if(c%4==0) printf "%s\n%s\n%s\n%s\n",h,s,p,$0 > d}
    END{printf "%d", n/1000000 > t}'
  TOT=$(cat "$P/tot" 2>/dev/null || echo 0)
  echo "[$(date +%T)] s$SLOT $S 원본 ${TOT} Mb · 준비 $(( $(date +%s)-T0 ))초"

  for F in 2 3 4; do
    M=$(( TOT / F )); [ "$M" -lt 60 ] && continue
    O="$P/f${F}.fq"; [ -s "$O" ] || continue
    W="$D/w_${S}_${F}_$SLOT"; rm -rf "$W"; mkdir -p "$W"
    ( cd "$W" && timeout 7200 ~/hifiasm --ont -o "$S" -t "$TH" "$O" > log 2>&1 )
    RC=$?
    R=$("$PY" ~/hiarr.py "${S}_sw${F}" "$W" "$D" 2>/dev/null)
    echo "$R" >> "$SUM"
    echo "    s$SLOT 1/$F (${M} Mb) rc=$RC : $(echo "$R" | tr '\n' ' | ')"
    rm -rf "$W"
    NH=$(awk -F'\t' -v s="$S" '$2=="OK" && index($1,s)==1{print $5}' "$SUM" \
         | sort -rn | awk '{if(!n||($0<last-1200)){print;last=$0;n++}}' | wc -l)
    [ "${NH:-0}" -ge 2 ] && { echo "    s$SLOT 두 하플로타입 확보 — 나머지 건너뜀"; break; }
  done
  rm -rf "$P"
done
echo "SWEEP31V2_s${SLOT}_DONE"
