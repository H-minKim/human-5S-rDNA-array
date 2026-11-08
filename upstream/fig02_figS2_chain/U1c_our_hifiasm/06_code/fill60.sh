#!/bin/bash
# Fill the 60 individuals that have no assembly candidate (2026-09-04)
#
#  Measured causes (60 without candidates among our 264 ONT individuals)
#    no attempt record at all  42  (median read data 408 MB, the data are sufficient)
#    NO_BOTH_PROBE          17  (run, but the array could not be joined into one contig)
#    NO_CTG                  1
#    Note: OOM (rc137) occurs 0 times among these 60 individuals.
#
#  Handling
#    A not attempted 42 : one hifiasm run with all reads
#    B NO_BOTH_PROBE 17 : fraction sweep (f=2..6), because too much data tangles the tandem graph
#      and reducing it was already measured to join the array (HG03460 521Mb failed -> 180Mb succeeded 87 units,
#      HG03521 615Mb failed -> 206Mb succeeded 74 units). The optimal data amount differs per haplotype, so
#      the union over several fractions is taken.
#
#  1kGP (kgp_ont/out) is not used, those data were judged unassemblable because they are 27 times too small.
#  Outputs are kept separately in ~/fill60.
set -u
MODE=${1:?A 또는 B}
D=~/fill60
mkdir -p "$D"
SUM=$D/fill60_summary.tsv
PY=""
for c in ~/miniconda3/bin/python python3; do
  p=$(command -v "$c" 2>/dev/null) || continue
  "$p" -c "import sys" 2>/dev/null && { PY="$p"; break; }
done
[ -n "$PY" ] || { echo "FATAL python 없음"; exit 9; }
LIST=~/fill60_$MODE.txt
[ -s "$LIST" ] || { echo "  목록 없음 $LIST"; exit 0; }

for S in $(tr -d '\r' < "$LIST"); do
  grep -q "^$S	" "$SUM" 2>/dev/null && { echo "  $S 이미 함"; continue; }
  Q=""
  for p in ~/s5s_ont/hits/$S.hits.fastq.gz ~/giab_hits/$S.hits.fastq.gz ~/cancer5s/hits/$S.hits.fastq.gz; do
    [ -f "$p" ] && { Q=$p; break; }
  done
  [ -n "$Q" ] || { echo "  $S read 없음"; continue; }
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  [ "${FREE:-0}" -lt 30 ] && { echo "  디스크 ${FREE}G — 중단"; break; }
  MB=$(( $(stat -c %s "$Q") / 1048576 ))
  echo "[$(date +%T)] $S ($MB MB) mode=$MODE"

  if [ "$MODE" = "A" ]; then
    W="$D/w_$S"; rm -rf "$W"; mkdir -p "$W"
    ( cd "$W" && timeout 10800 ~/hifiasm --ont -o "$S" -t 8 "$Q" > log 2>&1 )
    RC=$?
    R=$("$PY" ~/hiarr.py "$S" "$W" "$D" 2>/dev/null)
    echo "$R" >> "$SUM"
    echo "    rc=$RC : $(echo "$R" | tr '\n' ' | ')"
    rm -rf "$W"
  else
    TOT=$(zcat "$Q" | awk 'NR%4==2{t+=length($0)} END{printf "%d", t/1000000}')
    echo "    원본 ${TOT} Mb — 비율 sweep"
    for F in 2 3 4 5 6; do
      M=$(( TOT / F ))
      [ "$M" -lt 60 ] && continue
      O="$D/${S}_f${F}.fq.gz"
      # Subsample by order (cutting by length keeps only long reads and biases the result)
      zcat "$Q" | awk -v f="$F" 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{c++; if(c%f==0) print h"\n"s"\n"p"\n"$0}' | gzip > "$O"
      W="$D/w_${S}_$F"; rm -rf "$W"; mkdir -p "$W"
      ( cd "$W" && timeout 7200 ~/hifiasm --ont -o "$S" -t 8 "$O" > log 2>&1 )
      RC=$?
      R=$("$PY" ~/hiarr.py "${S}_f${F}" "$W" "$D" 2>/dev/null)
      echo "$R" >> "$SUM"
      echo "    1/$F (${M} Mb) rc=$RC : $(echo "$R" | tr '\n' ' | ')"
      rm -rf "$W" "$O"
    done
  fi
done
echo "FILL60_${MODE}_DONE"
