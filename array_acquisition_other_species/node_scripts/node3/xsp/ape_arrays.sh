#!/bin/bash
# 유인원 등 추가 종 : T2T 조립본 내려받기 -> Rfam RF00001 cmsearch -> 5S array 목록 (다른 종과 같은 규칙)
#   규칙은 종에 무관 : hit 사이 10 kb 이내면 한 array · 3 copy 이상 · copy 수 = hit 수 (bp/주기 환산 금지)
#   씨앗만 종마다 다르다 : 여기서는 공변 모형(RF00001)이 그 역할을 한다 (사람 NTS 씨앗은 쓰지 않는다)
#   인자 : "접근번호:이름" 목록 파일 · cmsearch 코어 수
#   산출 : ~/xsp/cms/<name>.tbl · ~/xsp/5S_array_catalog_apes.tsv (append)
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
X=~/xsp
L="${1:?목록}"
CPU="${2:-6}"
OUT=$X/5S_array_catalog_apes.tsv
mkdir -p $X/asm $X/cms $X/log
for line in $(cat "$L"); do
  acc=${line%%:*}
  nm=${line##*:}
  fa=$X/asm/$nm.fna
  if [ ! -s "$fa" ]; then
    curl -s -L -o $X/asm/$nm.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" || { echo "$nm 내려받기 실패"; continue; }
    unzip -o -q $X/asm/$nm.zip -d $X/asm/$nm.d && cat $X/asm/$nm.d/ncbi_dataset/data/$acc/*.fna > "$fa" && rm -rf $X/asm/$nm.d $X/asm/$nm.zip
  fi
  samtools faidx "$fa"
  echo "[$(date -u +%H:%M)] $nm $acc · $(awk '{s+=$2} END{printf "%.0f Mb · 서열 %d", s/1e6, NR}' $fa.fai)"
  if [ ! -s $X/cms/$nm.tbl ]; then
    cmsearch --cut_ga --cpu "$CPU" --noali --tblout $X/cms/$nm.tbl $X/rfam/RF00001.cm "$fa" > $X/log/$nm.cms.log 2>&1
  fi
  echo "[$(date -u +%H:%M)] $nm cmsearch hit $(grep -v -c '^#' $X/cms/$nm.tbl)"
  python3 $X/array_catalog.py "$nm" "$acc" $X/cms/$nm.tbl "$fa" "$OUT"
done
echo "끝 $(date -u +%H:%M) UTC"
column -t "$OUT" | head -40
