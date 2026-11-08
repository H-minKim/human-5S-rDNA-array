#!/bin/bash
# Transcribed from the run record: the code ran on a compute node and no copy was saved locally.
# Transcribed from the run record.
# hifiasm 0.25.0-r726 (node ~/hifiasm). Reads are the 5S reads extracted by U1a (~/s5s_ont/hits = HPRC ONT).
# 4 individuals where hifiasm failed to join the array into one contig - reassemble with downsampled data
#
#  Why downsample (measured evidence)
#    These four have 459-615 Mb of data, twice the successful case (hg002 289 Mb).
#    And when hg002 was reduced to 200 Mb with <50 kb reads,
#    it actually gave a cleaner 103/68 with no extra contigs.
#    In tandem repeats, excessive coverage can tangle the graph and split contigs.
#
#  Failure pattern (measured)
#    HG03460 has 4 contigs with only one probe attached, HG03868 3, HG03521 2
#    -> the array was assembled but the proximal and distal sides were split into different contigs
#    NA20756 has 0 probes at all -> a different problem
#
#  HiFi is not an alternative: these four have only 59-95 HiFi reads.
set -u
cd ~ || exit 9
D=~/resub
mkdir -p "$D"
for S in "$@"; do
  Q=~/s5s_ont/hits/$S.hits.fastq.gz
  [ -f "$Q" ] || { echo "  $S read 없음"; continue; }
  TOT=$(zcat "$Q" | awk 'NR%4==2{t+=length($0)} END{printf "%d", t/1000000}')
  echo "  $S 원본 ${TOT} Mb"
  # Uniform downsampling to a target of ~250 Mb (pick by read order to avoid length bias)
  for TGT in 250; do
    F=$(( (TOT + TGT - 1) / TGT ))
    [ "$F" -lt 1 ] && F=1
    O="$D/${S}_sub${TGT}.fq.gz"
    zcat "$Q" | awk -v f="$F" 'NR%4==1{h=$0} NR%4==2{s=$0} NR%4==3{p=$0} NR%4==0{c++; if(c%f==0) print h"\n"s"\n"p"\n"$0}' | gzip > "$O"
    N=$(zcat "$O" | awk 'NR%4==2{n++;t+=length($0)} END{printf "read %d · %.0f Mb", n, t/1e6}')
    echo "    1/$F 추림 -> $N"
    W="$D/${S}_$TGT"
    rm -rf "$W"; mkdir -p "$W"
    ( cd "$W" && timeout 7200 ~/hifiasm --ont -o "$S" -t 12 "$O" > log 2>&1 )
    echo "    조립 rc=$?"
    PY=""
    for c in ~/miniconda3/bin/python python3; do
      p=$(command -v "$c" 2>/dev/null) || continue
      "$p" -c "import sys" 2>/dev/null && { PY="$p"; break; }
    done
    "$PY" ~/hiarr.py "${S}_sub${TGT}" "$W" ~/resub
    rm -rf "$W" "$O"
  done
done
echo RESUB_DONE
