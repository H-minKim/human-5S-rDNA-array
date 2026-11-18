#!/bin/bash
# Arabidopsis: convert slices written by pod5 0.3.47 back into a format dorado 2.1.2 can read.
#   Cause: in a file written by 0.3.47 the channel column of the reads table is uint32, which dorado (and 0.3.44) reject. Written by 0.3.46 or earlier it reads.
#   Route: pod5 (0.3.47) --to_fast5--> fast5 --convert with 0.3.44--> pod5 (0.3.44). In testing dorado read these and the MM tags were intact.
#   Argument: the number of concurrent jobs.
#   Outputs: sub44/<same_name>.pod5 and sub44/converted.tsv (name, read count).
# MIT License. Copyright (c) 2026 Hyungmin Kim
set -u
H=~/xsp/ath
J="${1:-6}"
mkdir -p $H/sub44 $H/conv_log
: > $H/sub44/converted.tsv
one() {
  f="$1"
  b=$(basename "$f")
  H=~/xsp/ath
  [ -s $H/sub44/$b ] && return 0
  t=/tmp/conv_$b.d
  rm -rf $t
  mkdir -p $t/f5
  if ! ~/xsp/venv/bin/pod5 convert to_fast5 "$f" --output $t/f5 > $H/conv_log/$b.log 2>&1; then
    printf '%s\tFAIL_TO_FAST5\n' "$b" >> $H/sub44/converted.tsv
    rm -rf $t
    return 1
  fi
  if ! ~/xsp/venv44/bin/pod5 convert fast5 $t/f5/*.fast5 --output $H/sub44/$b >> $H/conv_log/$b.log 2>&1; then
    printf '%s\tFAIL_FROM_FAST5\n' "$b" >> $H/sub44/converted.tsv
    rm -rf $t
    return 1
  fi
  rm -rf $t
  n=$(~/xsp/venv44/bin/python -c "import pod5,sys; print(len(pod5.Reader(sys.argv[1]).read_ids))" $H/sub44/$b 2> /dev/null)
  printf '%s\t%s\n' "$b" "${n:-NA}" >> $H/sub44/converted.tsv
}
export -f one
ls $H/sub/*.pod5 | xargs -P "$J" -I{} bash -c 'one "$1"' _ {}
echo "conversion done $(date -u +%H:%M) UTC, files $(wc -l < $H/sub44/converted.tsv)"
awk -F'\t' '$2 ~ /^[0-9]+$/ {n++; s+=$2} $2 !~ /^[0-9]+$/ {f++} END{printf "succeeded %d files, read %d, failed %d\n", n, s, f+0}' $H/sub44/converted.tsv
echo "source(0.3.47) read sum $(awk -F'\t' 'NF>=5 && $5>0{s+=$5} END{print s}' $H/members.tsv)"
