#!/bin/bash
# 유인원 : GenomeArk ONT BAM (guppy 6.3.8 메틸 모델 재콜링본, MM/ML 있음) 에서 5S k-mer 가 든 read 만 골라낸다
#   재콜링 불필요 — 이미 메틸 태그가 붙어 있다 (README: "re-called using v6.3.8 with a model that included methylation")
#   흘려받기 : samtools view (https) -> grep -F (k-mer 8 개 = 종 5S 유전자 22-mer 4 개 + 역상보)
#   k-mer 는 그 종의 조립본 array 유닛에서 만든다 (struct/<species>.units.fa 첫 유닛의 유전자 119 bp)
#   인자 : GenomeArk 경로(예 Pan_troglodytes/mPanTro3) · k-mer 파일 · 출력폴더 · 동시 스트림 수
#   산출 : <flowcell>.hit.sam (태그 그대로) · <flowcell>.done (read 수 · MM 있는 수 · 파이프 종료코드 · 초)
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
SPP="${1:?GenomeArk 경로}"
KF="${2:?k-mer 파일}"
export O="${3:?출력폴더}"
P="${4:-2}"
mkdir -p "$O"
export KF
aws s3 ls --no-sign-request --region us-east-1 --recursive "s3://genomeark/species/$SPP/genomic_data/ont/" \
  | awk '$4 ~ /[.]bam$/ {print $4}' > "$O/bams.txt"
echo "BAM $(wc -l < "$O/bams.txt") 개"
one() {
  key="$1"
  b=$(basename "$key" .bam)
  [ -s "$O/$b.done" ] && return 0
  s=$(date +%s)
  samtools view -@3 "https://genomeark.s3.amazonaws.com/$key" 2> "$O/$b.err" | grep -F -f "$KF" > "$O/$b.hit.sam"
  st="${PIPESTATUS[*]}"
  n=$(wc -l < "$O/$b.hit.sam")
  m=$(grep -c "MM:Z" "$O/$b.hit.sam" || true)
  printf '%s\t%s\t%s\t%s\t%s\n' "$b" "$n" "$m" "$st" "$(( $(date +%s) - s ))" > "$O/$b.done"
}
export -f one
xargs -a "$O/bams.txt" -P "$P" -I{} bash -c 'one "$1"' _ {}
echo "선별 끝 $(date -u +%H:%M) UTC"
cat "$O"/*.done
