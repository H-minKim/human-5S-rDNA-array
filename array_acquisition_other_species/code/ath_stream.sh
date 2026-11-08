#!/bin/bash
# 애기장대 ATRun3 pod5 tar (EBI ERR17795792 · 164 GB) 를 디스크에 통째로 두지 않고 흘려 받으며 array read 만 pod5 로 남긴다
#   pdl.py (병렬 범위 받기 · 순서 보존 · 32 연결) -> tar -xz --to-command (멤버 하나씩 ath_member.sh)
#   인자 : MODE = test (앞 GB 만 · 전체 id 로 걸러 fastq id = pod5 id 확인) | full (ath.ids)
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
