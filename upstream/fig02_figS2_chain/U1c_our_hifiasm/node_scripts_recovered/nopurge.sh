#!/bin/bash
# hifiasm 0.25.0-r726 (node ~/hifiasm). Reads are the 5S reads extracted by U1a (~/s5s_ont/hits = HPRC ONT).
# Reassemble the same reads with -l 0 (purge off) to recover the shorter haplotype
#
#  Basis: in tandem arrays the two real haplotypes are 99.62% identical (measured).
#        The shorter one looks like a subsequence of the longer one, so purge_dups removes it as a haplotig.
#        There are 26 individuals for which hifiasm outputs only the longer one correctly and not the shorter one
#        (HG03704 only 136 units, spanning reads show 33 units, truth 33,136).
#
#  All other conditions are identical (--ont, same fastq, same threads).
set -u
TH=${TH:-10}
D=~/nopurge; mkdir -p "$D"; SUM=$D/nopurge_summary.tsv
PY=""; for c in ~/miniconda3/bin/python python3; do p=$(command -v "$c" 2>/dev/null) || continue
  "$p" -c "import sys" 2>/dev/null && { PY="$p"; break; }; done
for S in $(tr -d '\r' < ~/nopurge_list.txt); do
  grep -q "^$S	" "$SUM" 2>/dev/null && continue
  Q=""
  for p in ~/s5s_ont/hits/$S.hits.fastq.gz ~/giab_hits/$S.hits.fastq.gz ~/s5s_new/hits/$S.hits.fastq.gz; do
    [ -f "$p" ] && { Q=$p; break; }; done
  [ -n "$Q" ] || { echo "  $S no reads"; continue; }
  W="$D/w_$S"; rm -rf "$W"; mkdir -p "$W"
  echo "[$(date +%T)] $S : -l 0 (purge off)"
  ( cd "$W" && timeout 10800 ~/hifiasm --ont -l 0 -o "$S" -t "$TH" "$Q" > log 2>&1 )
  RC=$?
  R=$("$PY" ~/hiarr.py "${S}_l0" "$W" "$D" 2>/dev/null)
  echo "$R" >> "$SUM"
  echo "    rc=$RC : $(echo "$R" | tr '\n' ' | ')"
  rm -rf "$W"
done
echo NOPURGE_DONE
