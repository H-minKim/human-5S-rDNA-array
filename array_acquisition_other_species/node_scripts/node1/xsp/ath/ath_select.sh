#!/bin/bash
# Arabidopsis ATRun3 (PromethION R10.4.1, SUP v5.0.0 basecalls, no methylation tags): align the reads to Col-CC (PAF) and take the ids of reads overlapping a 5S array.
#   Alignment = minimap2 map-ont --secondary=no. A supplementary alignment also counts as overlapping. No MAPQ filter.
#   Finding the signal in the pod5 files requires the fastq id to equal the pod5 id, so the full id list is kept and split reads (pi:Z) are counted.
#   Outputs: all.ids, ath.ids, sup.colcc.paf.gz, select_summary.txt.
# MIT License. Copyright (c) 2026 Hyungmin Kim
set -u
H=~/xsp/ath
A=~/xsp/asm
BED=$H/arabidopsis_ColCC_5S_arrays.bed
TH="${1:-10}"
zcat $H/ATRun3_SUP.fastq.gz | awk 'NR%4==1{print substr($1,2); if($0 ~ /pi:Z:/) p++} END{print p+0 > "/dev/stderr"}' > $H/all.ids 2> $H/pi_count.txt
minimap2 -x map-ont -t "$TH" --secondary=no $A/arabidopsis_ColCC.map-ont.mmi $H/ATRun3_SUP.fastq.gz 2> $H/mm2.log | gzip -1 > $H/sup.colcc.paf.gz
zcat $H/sup.colcc.paf.gz | awk -F'\t' 'NR==FNR{c[NR]=$1; s[NR]=$2; e[NR]=$3; n=NR; next} {for(i=1;i<=n;i++) if($6==c[i] && $8<e[i] && $9>s[i]) {print $1; break}}' $BED - | sort -u > $H/ath.ids
st=$(zcat $H/sup.colcc.paf.gz | awk -F'\t' 'NR==FNR{w[$1]; next} ($1 in w) && !($1 in seen){seen[$1]; n++; s+=$2} END{printf "%d reads · %d bases", n, s}' $H/ath.ids -)
{
  echo "total read $(wc -l < $H/all.ids), split read(pi:Z) $(cat $H/pi_count.txt)"
  echo "5S array read : $st"
  echo "done $(date -u +%H:%M) UTC"
} | tee $H/select_summary.txt
