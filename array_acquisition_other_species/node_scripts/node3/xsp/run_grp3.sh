#!/bin/bash
# The second-stage call on unit types within an array, over all species.
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
