#!/bin/bash
# Preparing the ape alignments (on the compute node): download and index the assembly. The cmsearch tables, structure tables and modkit output are moved in by the driver script.
#   Arguments: species key (e.g. gorilla_mGorGor1), accession (e.g. GCA_029281585.3).
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
