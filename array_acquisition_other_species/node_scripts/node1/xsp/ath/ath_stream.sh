#!/bin/bash
# The Arabidopsis ATRun3 pod5 tar (EBI ERR17795792, 164 GB) is streamed rather than stored whole, keeping only the array reads as pod5.
#   pdl.py (parallel range requests, order preserved, 32 connections) -> tar -xz --to-command, one member at a time through ath_member.sh.
#   Argument: MODE = test (the first few GB only, filtered against the full id list to confirm fastq id = pod5 id) or full (ath.ids).
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
H=~/xsp/ath
MODE="${1:?test|full}"
U=https://ftp.sra.ebi.ac.uk/vol1/run/ERR177/ERR17795792/ATRun3_pod5.tar.gz
if [ "$MODE" = test ]; then
  export IDS=$H/all.ids
  : > $H/members.tsv
  python3 ~/xsp/pdl.py $U - 32 16 2> $H/pdl_pod5_test.log | head -c "${2:-3000000000}" | tar -xzf - --to-command="bash $H/ath_member.sh" 2> $H/tar_test.log
  echo "시험 끝 rc=${PIPESTATUS[*]}"
  cat $H/members.tsv
  rm -rf $H/sub
else
  export IDS=$H/ath.ids
  : > $H/members.tsv
  python3 ~/xsp/pdl.py $U - 32 16 2> $H/pdl_pod5.log | tar -xzf - --to-command="bash $H/ath_member.sh" 2> $H/tar_full.log
  echo "전체 끝 rc=${PIPESTATUS[*]} $(date -u +%H:%M) UTC"
  awk -F'\t' '{n++; a+=$4; b+=$5; if($3!=0) f++} END{print "멤버 " n " · read " a " · 남김 " b " · filter 실패 " f+0}' $H/members.tsv
fi
