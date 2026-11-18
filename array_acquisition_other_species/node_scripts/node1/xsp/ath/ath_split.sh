#!/bin/bash
# Arabidopsis: distribute the pod5 slices holding only the array reads across the nodes, by size, weighted by cores at 16, 28 and 16.
#   Only slices that kept at least one read, where column 5 of members.tsv is above 0.
#   writes : grp_A.txt, grp_W.txt, grp_B.txt (lists of absolute paths; the three nodes share one home directory, so the paths are the same on all of them)
set -u
H=~/xsp/ath
rm -f $H/grp_A.txt $H/grp_W.txt $H/grp_B.txt
awk -F'\t' '$5>0 {print $1}' $H/members.tsv | while read b; do
  f=$H/sub/$b
  [ -s "$f" ] && printf '%s\t%s\n' "$(stat -c %s "$f")" "$f"
done | sort -k1,1nr | awk -F'\t' -v H="$H" 'BEGIN{w["A"]=16; w["W"]=28; w["B"]=16}
  {best=""; for(k in w){r=(ld[k]+$1)/w[k]; if(best=="" || r<br){best=k; br=r}}; ld[best]+=$1; n[best]++; print $2 > (H "/grp_" best ".txt")}
  END{for(k in w) printf "%s : chunks %d, %.2f GB\n", k, n[k], ld[k]/1e9}'
