#!/bin/bash
# fill60 v2, shortened with parallel runs and an adaptive sweep (2026-09-04)
#
#  v1 measured: median 10.8 min per individual. Remaining assemblies node58 39 runs (7.0h), 32 29 runs, 57 22 runs
#
#  Two ways to cut it
#   (1) run in parallel to match the node cores
#      extraction of 1 individual = 3-4 cores, 2 extractions per node need 6-8 cores left free.
#      32-core node : 2 hifiasm x -t 12 = 24 cores (+extraction 8) = 32  -> 2 in parallel
#      16-core node : 1 hifiasm x -t 8  =  8 cores (+extraction 8) = 16  -> 1 in parallel
#   (2) adaptive sweep: if two haplotypes come out at f=2, skip f=3,4.
#      The goal of mode B is to obtain two haplotypes, so once they are obtained there is no reason to keep running.
#      (subsampling is done because too much data tangles the tandem graph, and all measured recoveries were at f=2-3)
#
#  1kGP is not used. Output goes to ~/fill60.
set -u
MODE=${1:?A 또는 B}
SLOT=${2:-0}          # parallel slot number (0,1), the list is split between slots
NSLOT=${3:-1}
D=~/fill60
mkdir -p "$D"
SUM=$D/fill60_summary.tsv
TH=${TH:-8}
PY=""
for c in ~/miniconda3/bin/python python3; do
  p=$(command -v "$c" 2>/dev/null) || continue
  "$p" -c "import sys" 2>/dev/null && { PY="$p"; break; }
done
[ -n "$PY" ] || { echo "FATAL python 없음"; exit 9; }
LIST=~/fill60_$MODE.txt
[ -s "$LIST" ] || { echo "  목록 없음"; exit 0; }

i=-1
for S in $(tr -d '\r' < "$LIST"); do
  i=$((i+1))
  [ $((i % NSLOT)) -ne "$SLOT" ] && continue
  grep -q "^$S	" "$SUM" 2>/dev/null && continue
  Q=""
  for p in ~/s5s_ont/hits/$S.hits.fastq.gz ~/giab_hits/$S.hits.fastq.gz ~/cancer5s/hits/$S.hits.fastq.gz; do
    [ -f "$p" ] && { Q=$p; break; }
  done
  [ -n "$Q" ] || { echo "  $S read 없음"; continue; }
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  [ "${FREE:-0}" -lt 30 ] && { echo "  디스크 ${FREE}G — 중단"; break; }
  MB=$(( $(stat -c %s "$Q") / 1048576 ))
  echo "[$(date +%T)] s$SLOT $S ($MB MB) mode=$MODE"

  if [ "$MODE" = "A" ]; then
    W="$D/w_${S}_$SLOT"; rm -rf "$W"; mkdir -p "$W"
    ( cd "$W" && timeout 10800 ~/hifiasm --ont -o "$S" -t "$TH" "$Q" > log 2>&1 )
    RC=$?
    R=$("$PY" ~/hiarr.py "$S" "$W" "$D" 2>/dev/null)
    echo "$R" >> "$SUM"
    echo "    s$SLOT rc=$RC : $(echo "$R" | tr '\n' ' | ')"
    rm -rf "$W"
  else
    TOT=$(zcat "$Q" | awk 'NR%4==2{t+=length($0)} END{printf "%d", t/1000000}')
    for F in 2 3 4; do
      M=$(( TOT / F ))
      [ "$M" -lt 60 ] && continue
      O="$D/${S}_f${F}_$SLOT.fq.gz"
      zcat "$Q" | awk -v f="$F" 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{c++; if(c%f==0) print h"\n"s"\n"p"\n"$0}' | gzip > "$O"
      W="$D/w_${S}_${F}_$SLOT"; rm -rf "$W"; mkdir -p "$W"
      ( cd "$W" && timeout 7200 ~/hifiasm --ont -o "$S" -t "$TH" "$O" > log 2>&1 )
      RC=$?
      R=$("$PY" ~/hiarr.py "${S}_f${F}" "$W" "$D" 2>/dev/null)
      echo "$R" >> "$SUM"
      echo "    s$SLOT 1/$F (${M} Mb) rc=$RC : $(echo "$R" | tr '\n' ' | ')"
      rm -rf "$W" "$O"
      # adaptive: stop once 2 or more arrays of different lengths have come out for this individual
      NH=$(awk -F'\t' -v s="$S" '$2=="OK" && ($1==s || index($1,s"_f")==1){print $5}' "$SUM" \
           | sort -rn | awk '{if(!n||($0<last-1200)){print;last=$0;n++}}' | wc -l)
      if [ "${NH:-0}" -ge 2 ]; then
        echo "    s$SLOT 두 하플로타입 확보 — 나머지 비율 건너뜀"
        break
      fi
    done
  fi
done
echo "FILL60V2_${MODE}_s${SLOT}_DONE"
