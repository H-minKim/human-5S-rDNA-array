#!/bin/bash
# The second-stage call on unit types within an array, over all species.
cd ~/xsp
mkdir -p grp2
for d in mapcur/*/; do
  sp=$(basename "$d")
  fa=asm/${sp}.fna
  [ -f "$fa" ] || { echo "SKIP $sp (fasta none)"; continue; }
  echo "=== $sp ==="
  python3 unit_groups2.py "$sp" "$fa" "$d" grp2 50
done
echo DONE_ALL
