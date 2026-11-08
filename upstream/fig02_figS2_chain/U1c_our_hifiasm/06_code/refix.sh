#!/bin/bash
# Reassemble with HPRC ONT the individuals that were assembled from low-quality 1kGP reads (2026-09-04)
#
#  cause (measured)
#    the search order of find_reads() in hifi_queue.sh was
#      ~/kgp_ont/out/$s.array.fastq.gz   <- 1kGP (10-24 MB) first
#      ~/s5s_ont/hits/$s.hits.fastq.gz   <- HPRC ONT (250-670 MB) later
#    so individuals present in both cohorts were assembled from 1kGP data.
#    1kGP is data already judged "not assemblable" (27x too little data, only 1 of 12 haplotypes with known truth was correct).
#    impact: 17 of the 607 hifi_queue outputs (12 individuals). 6 of those individuals have no replacement candidate.
#
#  Here the order is reversed and HPRC ONT is used first. 1kGP is not used at all.
#  Outputs go separately to ~/refix. They must not mix with the old 1kGP outputs.
set -u
D=~/refix
mkdir -p "$D"
SUM=$D/refix_summary.tsv
PY=""
for c in ~/miniconda3/bin/python python3; do
  p=$(command -v "$c" 2>/dev/null) || continue
  "$p" -c "import sys" 2>/dev/null && { PY="$p"; break; }
done
[ -n "$PY" ] || { echo "FATAL python 없음"; exit 9; }

for S in "$@"; do
  grep -q "^$S	" "$SUM" 2>/dev/null && { echo "  $S 이미 함"; continue; }
  # 1kGP is not used
  Q=""
  for p in ~/s5s_ont/hits/$S.hits.fastq.gz ~/giab_hits/$S.hits.fastq.gz ~/cancer5s/hits/$S.hits.fastq.gz; do
    [ -f "$p" ] && { Q=$p; break; }
  done
  [ -n "$Q" ] || { echo "  $S HPRC/GIAB read 없음 — 건너뜀"; continue; }
  MB=$(( $(stat -c %s "$Q") / 1048576 ))
  echo "  $S : $Q (${MB} MB)"
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  [ "${FREE:-0}" -lt 25 ] && { echo "    디스크 ${FREE}G — 중단"; break; }
  W="$D/w_$S"
  rm -rf "$W"; mkdir -p "$W"
  ( cd "$W" && timeout 10800 ~/hifiasm --ont -o "$S" -t 8 "$Q" > log 2>&1 )
  echo "    조립 rc=$?"
  R=$("$PY" ~/hiarr.py "$S" "$W" "$D" 2>/dev/null)
  echo "$R" >> "$SUM"
  echo "    $(echo "$R" | tr '\n' ' | ')"
  rm -rf "$W"
done
echo REFIX_DONE
