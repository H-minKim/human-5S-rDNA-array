#!/bin/bash
# array 안 유닛 종류 판정 2 판 — 전 종
cd ~/xsp
mkdir -p grp3
for d in mapcur/*/; do
  sp=$(basename "$d")
  fa=asm/${sp}.fna
  [ -f "$fa" ] || { echo "SKIP $sp (fasta 없음)"; continue; }
  echo "=== $sp ==="
  python3 unit_groups3.py "$sp" "$fa" "$d" grp3 50
done
echo DONE_ALL
