#!/bin/bash
# Transcribed from the run record: the code ran on a compute node and no copy was saved locally.
# Transcribed from the run record.
# hifiasm 0.25.0-r726 (node ~/hifiasm). Reads are the 5S reads extracted by U1a (~/s5s_ont/hits = HPRC ONT).
# Try several subsampling fractions to obtain both haplotypes
#
#  Basis (measured): too much data tangles the tandem graph.
#    HG03460 521 Mb failed -> 180 Mb succeeded (87, matches tier1)
#    HG03521 615 Mb failed -> 206 Mb succeeded (74, matches tier1)
#    HG03868 459 Mb failed -> 230 Mb succeeded (116 and 61, matching tier2 and tier1 respectively)
#  But the first two yielded only one haplotype. The optimal data amount may differ
#  per haplotype, so several fractions are scanned and the union is taken.
#
#  Subsampling is done by order (c%f), cutting by length keeps only long reads and biases the result.
set -u
cd ~ || exit 9
D=~/resub2
mkdir -p "$D"
OUT=~/resub2_summary.tsv
for S in "$@"; do
  Q=~/s5s_ont/hits/$S.hits.fastq.gz
  [ -f "$Q" ] || { echo "  $S read 없음"; continue; }
  TOT=$(zcat "$Q" | awk 'NR%4==2{t+=length($0)} END{printf "%d", t/1000000}')
  echo "== $S 원본 ${TOT} Mb =="
  for F in 2 3 4 5 6; do
    MB=$(( TOT / F ))
    [ "$MB" -lt 60 ] && continue
    O="$D/${S}_f${F}.fq.gz"
    zcat "$Q" | awk -v f="$F" 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{c++; if(c%f==0) print h"\n"s"\n"p"\n"$0}' | gzip > "$O"
    W="$D/${S}_$F"
    rm -rf "$W"; mkdir -p "$W"
    ( cd "$W" && timeout 5400 ~/hifiasm --ont -o "$S" -t 10 "$O" > log 2>&1 )
    RC=$?
    PY=""
    for c in ~/miniconda3/bin/python python3; do
      p=$(command -v "$c" 2>/dev/null) || continue
      "$p" -c "import sys" 2>/dev/null && { PY="$p"; break; }
    done
    R=$("$PY" ~/hiarr.py "${S}_f${F}" "$W" "$D" 2>/dev/null)
    echo "  1/$F (${MB} Mb) rc=$RC : $(echo "$R" | tr '\n' ' | ')"
    echo "$R" >> "$OUT"
    rm -rf "$W" "$O"
  done
done
echo RESUB2_DONE
