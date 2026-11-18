#!/bin/bash
# Reassemble with HPRC ONT the individuals that hifi_queue.sh assembles from low-quality 1kGP reads
#
#  HPRC ONT is used first. 1kGP is not used at all.
#  Outputs go separately to ~/refix. They must not mix with the 1kGP outputs.
set -u
D=~/refix
mkdir -p "$D"
SUM=$D/refix_summary.tsv
PY=""
for c in ~/miniconda3/bin/python python3; do
  p=$(command -v "$c" 2>/dev/null) || continue
  "$p" -c "import sys" 2>/dev/null && { PY="$p"; break; }
done
[ -n "$PY" ] || { echo "FATAL python not found"; exit 9; }

for S in "$@"; do
  grep -q "^$S	" "$SUM" 2>/dev/null && { echo "  $S already done"; continue; }
  # 1kGP is not used
  Q=""
  for p in ~/s5s_ont/hits/$S.hits.fastq.gz ~/giab_hits/$S.hits.fastq.gz ~/cancer5s/hits/$S.hits.fastq.gz; do
    [ -f "$p" ] && { Q=$p; break; }
  done
  [ -n "$Q" ] || { echo "  $S no HPRC/GIAB read, skipped"; continue; }
  MB=$(( $(stat -c %s "$Q") / 1048576 ))
  echo "  $S : $Q (${MB} MB)"
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  [ "${FREE:-0}" -lt 25 ] && { echo "    disk ${FREE}G: aborting"; break; }
  W="$D/w_$S"
  rm -rf "$W"; mkdir -p "$W"
  ( cd "$W" && timeout 10800 ~/hifiasm --ont -o "$S" -t 8 "$Q" > log 2>&1 )
  echo "    assembly rc=$?"
  R=$("$PY" ~/hiarr.py "$S" "$W" "$D" 2>/dev/null)
  echo "$R" >> "$SUM"
  echo "    $(echo "$R" | tr '\n' ' | ')"
  rm -rf "$W"
done
echo REFIX_DONE
