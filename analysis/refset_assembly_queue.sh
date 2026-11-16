#!/bin/bash
# Reassemble the 5S array with hifiasm --ont, as a queue on the compute node.
#   Step: assembly.
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md
# Reassemble the 5S array with hifiasm --ont, for 1kGP, GIAB, HG008 and the cancer samples.
#
#  The verified premise, measured:
#    across 5 haplotypes of 3 HPRC individuals, our hifiasm --ont assembly agreed
#    with the official HPRC assembly base for base (200410/88, 117677/51, 126824/55,
#    167029/73, and 185018 against 185024/81). The method itself is sound.
#
#  Disk guard: hifiasm uses 676 MB per individual and the node has only 43 GB, so
#    the intermediates are deleted as soon as each individual finishes, keeping only the array sequence and the summary.
#    Running out of memory or disk corrupts the output silently, so each step is checked as it goes.
#
#  Usage: hifi_queue.sh <list file> <output tag> <threads>
set -u
LIST=$1
TAG=$2
TH=${3:-6}
OUT=~/hiasm
mkdir -p $OUT/arr $OUT/log
SUM=$OUT/summary_$TAG.tsv

find_reads () {
  local s=$1
  for p in ~/kgp_ont/out/$s.array.fastq.gz \
           ~/s5s_ont/hits/$s.hits.fastq.gz \
           ~/giab_hits/$s.hits.fastq.gz \
           ~/cancer5s/hits/$s.hits.fastq.gz ; do
    [ -f "$p" ] && { echo "$p"; return; }
  done
  echo ""
}

for S in $(tr -d '\r' < "$LIST"); do
  # individuals already done are skipped, which makes a restart safe
  grep -q "^$S	" $SUM 2>/dev/null && continue
  Q=$(find_reads "$S")
  if [ -z "$Q" ]; then
    printf "%s\tNO_READS\t-\t-\t-\t-\n" "$S" >> $SUM
    continue
  fi
  # the free disk is checked and the run stops below 5 GB, to prevent silent corruption
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$FREE" -lt 5 ]; then
    echo "디스크 부족 ${FREE}G — 중단" >> $OUT/log/$TAG.log
    break
  fi
  D=$(mktemp -d ~/hitmp_XXXXXX)
  ( cd "$D" && timeout 7200 ~/hifiasm --ont -o "$S" -t "$TH" "$Q" > hifiasm.log 2>&1 )
  RC=$?
  NC=$(ls "$D"/*.p_ctg.gfa 2>/dev/null | wc -l)
  if [ "$RC" -ne 0 ] || [ "$NC" -eq 0 ]; then
    printf "%s\tFAIL_rc%s\t-\t-\t-\t-\n" "$S" "$RC" >> $SUM
    tail -3 "$D"/hifiasm.log >> $OUT/log/$TAG.log 2>/dev/null
    rm -rf "$D"
    continue
  fi
  # array extraction, by the same method as hicheck
  ~/miniconda3/bin/python ~/hiarr.py "$S" "$D" "$OUT/arr" >> $SUM 2>> $OUT/log/$TAG.log
  rm -rf "$D"
done
echo "QUEUE_DONE $TAG"
