#!/bin/bash
cd ~/xsp
mkdir -p utypes
for d in mapcur/*/; do
  sp=$(basename "$d"); fa=asm/${sp}.fna
  [ -f "$fa" ] || { echo "SKIP $sp"; continue; }
  python3 unit_types.py "$sp" "$fa" "$d" utypes 5S_unit_status.tsv 50 20
done
echo DONE_ALL
