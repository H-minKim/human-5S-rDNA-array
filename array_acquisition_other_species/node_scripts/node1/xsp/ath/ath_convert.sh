#!/bin/bash
# 애기장대 : pod5 0.3.47 로 쓴 조각을 dorado 2.1.2 가 읽는 형식으로 되돌린다
#   원인 : 0.3.47 이 쓴 파일은 reads 표의 channel 이 uint32 라 dorado(및 0.3.44)가 거부한다. 0.3.46 이하로 쓰면 읽힌다
#   길 : pod5(0.3.47) --to_fast5--> fast5 --0.3.44 convert--> pod5(0.3.44)  [시험에서 dorado 가 읽고 MM 태그까지 정상]
#   인자 : 동시 작업 수
#   산출 : sub44/<same_name>.pod5 · sub44/converted.tsv (이름 · read 수)
# MIT License · Copyright (c) 2026 Hungmin Kim
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
echo "변환 끝 $(date -u +%H:%M) UTC · 파일 $(wc -l < $H/sub44/converted.tsv)"
awk -F'\t' '$2 ~ /^[0-9]+$/ {n++; s+=$2} $2 !~ /^[0-9]+$/ {f++} END{printf "성공 %d 파일 · read %d · 실패 %d\n", n, s, f+0}' $H/sub44/converted.tsv
echo "원본(0.3.47) read 합 $(awk -F'\t' 'NF>=5 && $5>0{s+=$5} END{print s}' $H/members.tsv)"
