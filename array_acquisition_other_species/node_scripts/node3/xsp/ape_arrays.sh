#!/bin/bash
# Additional species including the great apes: download the T2T assembly -> cmsearch with Rfam RF00001 -> the 5S array catalogue, under the same rules as the other species.
#   The rules do not depend on the species: hits within 10 kb form one array, at least 3 copies, and copy number equals the hit count. It is not derived from bp / period.
#   Only the seed differs between species; here the covariance model RF00001 plays that role, and the human NTS seed is not used.
#   Arguments: a list file of "accession:name", and the number of cmsearch cores.
#   Outputs: cms/<name>.tbl and 5S_array_catalog_apes.tsv (appended to).
# MIT License. Copyright (c) 2026 Hyungmin Kim
set -u
X=~/xsp
L="${1:?list}"
CPU="${2:-6}"
OUT=$X/5S_array_catalog_apes.tsv
mkdir -p $X/asm $X/cms $X/log
for line in $(cat "$L"); do
  acc=${line%%:*}
  nm=${line##*:}
  fa=$X/asm/$nm.fna
  if [ ! -s "$fa" ]; then
    curl -s -L -o $X/asm/$nm.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" || { echo "$nm download failed"; continue; }
    unzip -o -q $X/asm/$nm.zip -d $X/asm/$nm.d && cat $X/asm/$nm.d/ncbi_dataset/data/$acc/*.fna > "$fa" && rm -rf $X/asm/$nm.d $X/asm/$nm.zip
  fi
  samtools faidx "$fa"
  echo "[$(date -u +%H:%M)] $nm $acc, $(awk '{s+=$2} END{printf "%.0f Mb, sequences %d", s/1e6, NR}' $fa.fai)"
  if [ ! -s $X/cms/$nm.tbl ]; then
    cmsearch --cut_ga --cpu "$CPU" --noali --tblout $X/cms/$nm.tbl $X/rfam/RF00001.cm "$fa" > $X/log/$nm.cms.log 2>&1
  fi
  echo "[$(date -u +%H:%M)] $nm cmsearch hit $(grep -v -c '^#' $X/cms/$nm.tbl)"
  python3 $X/array_catalog.py "$nm" "$acc" $X/cms/$nm.tbl "$fa" "$OUT"
done
echo "finished $(date -u +%H:%M) UTC"
column -t "$OUT" | head -40
