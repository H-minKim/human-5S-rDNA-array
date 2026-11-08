#!/bin/bash
# Arabidopsis: claim one batch of filtered, unclaimed pod5 chunks, one batch at a time under flock
#   members.tsv 한 줄 = 멤버 필터 완료 (5 열 = 남긴 read 수 · 0 이면 제외)
#   인자 : 묶음 read 목표 수 · 강제(1 = 목표 미달이어도 남은 것 모두 · 흘려받기가 끝난 뒤)
#   출력 : 묶음 이름 (b001 ...) 또는 NONE   · 산출 : batch/<batch>.txt (절대경로 목록) · batch/claimed.txt
set -u
H=~/xsp/ath
N="${1:-300}"
FORCE="${2:-0}"
mkdir -p $H/batch
touch $H/batch/claimed.txt
exec 9> $H/batch/.lock
flock 9
# NR==FNR 는 첫 파일(claimed.txt)이 비면 둘째 파일을 첫 파일로 착각한다 -> 파일 이름으로 가른다
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
