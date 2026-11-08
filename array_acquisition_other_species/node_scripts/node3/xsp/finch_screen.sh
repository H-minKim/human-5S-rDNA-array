#!/bin/bash
# 금화조 ONT fastq 에서 5S 유전자 k-mer 가 든 read 만 고르기 (정렬 전 선별)
#   스트리밍 (S3 -> pigz -> 한 줄 묶기 -> grep -F) · 디스크에는 후보 read 만
#   k-mer = v0.3 mat chr2 5S 유전자 (cmsearch hit CM109790.1:97272500-97272617, - 가닥 = 유전자 5'->3') 의 22-mer 4 개 + 역상보
#     array read 는 5S 사본이 수백 개라 22-mer 정확일치가 사실상 반드시 걸린다 (R10 ~Q20). 흩어진 5S 는 뒤 정렬에서 걸러진다
#   네 줄 묶기 구분자 = 제어문자 0x01 (헤더에 탭이 있어도 레코드가 깨지지 않게)
#   인자 : S3 키 목록 파일 · 동시 스트림 수
#   산출 : ~/xsp/zfont/cand/<flowcell>.cand.fq.gz · <flowcell>.done (후보 수 · 파이프 각 단계 종료코드)
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
L="${1:?키목록}"
P="${2:-4}"
export O=~/xsp/zfont/cand
mkdir -p "$O"
K=~/xsp/zfont/k.txt
printf '%s\n' GATCTCGGAAGCTAAGCAGGGT GGGCCCGGTTAGTACTTGGATG GGAGACCGCCTGGGAATACCGG GCCGCCCTGAGAGCGGCCGGTC > "$K.fw"
rev "$K.fw" | tr ACGT TGCA > "$K.rc"
cat "$K.fw" "$K.rc" > "$K"
export K
one() {
  k="$1"
  b=$(basename "$k" .raw.fastq.gz)
  [ -s "$O/$b.done" ] && return 0
  aws s3 cp --no-sign-request --region us-east-1 "s3://genomeark/species/Taeniopygia_guttata/bTaeGut7/genomic_data/ont/$k" - 2> "$O/$b.err" \
    | pigz -dc -p 3 | paste -d $'\001\001\001' - - - - | grep -F -f "$K" | tr '\001' '\n' | gzip -1 > "$O/$b.cand.fq.gz"
  st="${PIPESTATUS[*]}"
  n=$(( $(zcat "$O/$b.cand.fq.gz" | wc -l) / 4 ))
  printf '%s\t%s\t%s\n' "$b" "$n" "$st" > "$O/$b.done"
}
export -f one
xargs -a "$L" -P "$P" -I{} bash -c 'one "$1"' _ {}
echo "선별 끝 $(date -u +%H:%M) UTC"
cat "$O"/*.done
