#!/bin/bash
# Arabidopsis: claim one batch of filtered, unclaimed pod5 chunks, one batch at a time under flock
#   One row of members.tsv means that member has been filtered; column 5 is the number of reads kept, and 0 means it is excluded.
#   Arguments: the target number of reads per batch, and force (1 = take everything left even below the target, once the streaming has finished).
#   Prints the batch name (b001 ...) or NONE. Outputs: batch/<batch>.txt (a list of absolute paths) and batch/claimed.txt.
set -u
H=~/xsp/ath
N="${1:-300}"
FORCE="${2:-0}"
mkdir -p $H/batch
touch $H/batch/claimed.txt
exec 9> $H/batch/.lock
flock 9
# NR==FNR mistakes the second file for the first when the first (claimed.txt) is empty, so the files are distinguished by name instead.
avail=$(awk -F'\t' 'FILENAME==ARGV[1]{c[$1]; next} NF>=5 && $5>0 && !($1 in c){print $1"\t"$5}' $H/batch/claimed.txt $H/members.tsv)
tot=$(printf '%s\n' "$avail" | awk -F'\t' 'NF==2{s+=$2} END{print s+0}')
if [ "$tot" -eq 0 ] || { [ "$tot" -lt "$N" ] && [ "$FORCE" != 1 ]; }; then
  echo NONE
  exit 0
fi
k=$(( $(ls $H/batch/b*.txt 2> /dev/null | wc -l) + 1 ))
id=$(printf 'b%03d' $k)
printf '%s\n' "$avail" | awk -F'\t' -v N="$N" -v H="$H" 'NF==2 && s<N {print H "/sub44/" $1; s+=$2}' > $H/batch/$id.txt
sed 's#.*/##' $H/batch/$id.txt >> $H/batch/claimed.txt
echo $id
