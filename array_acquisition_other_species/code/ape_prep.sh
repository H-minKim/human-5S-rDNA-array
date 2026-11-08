#!/bin/bash
# 유인원 정렬 준비 (노드에서) : 조립본 내려받기 + 색인. cmsearch 표·구조표·modkit 은 지휘 스크립트가 옮겨준다
#   인자 : 종키(예 gorilla_mGorGor1) · 접근번호(예 GCA_029281585.3)
set -u
SP="${1:?종키}"
ACC="${2:?접근번호}"
X=~/xsp
mkdir -p $X/asm $X/cms $X/struct $X/ape
cd $X/asm
if [ ! -s $SP.fna ]; then
  curl -s -L -o $SP.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$ACC/download?include_annotation_type=GENOME_FASTA" || exit 1
  unzip -o -q $SP.zip -d $SP.d && cat $SP.d/ncbi_dataset/data/$ACC/*.fna > $SP.fna && rm -rf $SP.d $SP.zip
fi
[ -s $SP.fna.fai ] || samtools faidx $SP.fna
echo "$(hostname) $SP $(awk '{s+=$2} END{printf "%.0f Mb · 서열 %d", s/1e6, NR}' $SP.fna.fai) · 여유 $(df -h ~ | awk 'NR==2{print $4}')"
