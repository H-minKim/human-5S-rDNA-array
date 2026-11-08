#!/bin/bash
# Reassemble the 5S array with hifiasm --ont (1kGP, GIAB, HG008, cancer)
#
#  verified premise (measured 2026-09-02):
#    in 5 haplotypes of 3 HPRC individuals our hifiasm --ont assemblies
#    matched the official HPRC assemblies at bp level (200410/88, 117677/51, 126824/55,
#    167029/73, 185018 vs 185024/81). The method itself is valid.
#
#  disk protection: hifiasm uses 676 MB per individual. The node has only 43 GB, so
#    intermediates are deleted as soon as each individual finishes, keeping only the array sequence and the summary.
#    (there were earlier silent corruptions from OOM/disk shortage, so it always checks as it runs)
#
#  usage: hifi_queue.sh <list file> <output tag> <threads>
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
  # skip individuals already done (restart safe)
  grep -q "^$S	" $SUM 2>/dev/null && continue
  Q=$(find_reads "$S")
  if [ -z "$Q" ]; then
    printf "%s\tNO_READS\t-\t-\t-\t-\n" "$S" >> $SUM
    continue
  fi
  # check free disk, stop below 5 GB (prevents silent corruption)
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
  # array extraction (same method as hicheck)
  ~/miniconda3/bin/python ~/hiarr.py "$S" "$D" "$OUT/arr" >> $SUM 2>> $OUT/log/$TAG.log
  rm -rf "$D"
done
echo "QUEUE_DONE $TAG"
