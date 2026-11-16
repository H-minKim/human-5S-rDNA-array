#!/bin/bash
# Zebra finch: the 5S k-mer candidate reads are aligned to the v0.3 diploid (maternal + paternal) reference (PAF), and the ids of reads overlapping a 5S array at least once are kept per flowcell.
#   Alignment = minimap2 map-ont --secondary=no; a supplementary alignment also counts as overlapping, and there is no MAPQ filter because this is a tandem array.
#   The arrays are zfinch_v03_5S_arrays.bed: cmsearch RF00001 hits grouped within 10 kb with at least 3 copies, 4 maternal and 4 paternal.
#   Outputs: sel/<fc>.ids, sel/<fc>.paf.gz and sel/summary.tsv (flowcell, candidate reads, array reads, total bases in the array reads).
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
Z=~/xsp/zfont
A=~/xsp/asm
BED=$Z/zfinch_v03_5S_arrays.bed
TH="${1:-12}"
mkdir -p $Z/sel
: > $Z/sel/summary.tsv
for f in $Z/cand/*.cand.fq.gz; do
  b=$(basename "$f" .cand.fq.gz)
  minimap2 -x map-ont -t "$TH" --secondary=no $A/zfinch_v03_dip.map-ont.mmi "$f" 2> /dev/null | gzip -1 > $Z/sel/$b.paf.gz
  zcat $Z/sel/$b.paf.gz | awk -F'\t' 'NR==FNR{c[NR]=$1; s[NR]=$2; e[NR]=$3; n=NR; next} {for(i=1;i<=n;i++) if($6==c[i] && $8<e[i] && $9>s[i]) {print $1; break}}' $BED - | sort -u > $Z/sel/$b.ids
  st=$(zcat "$f" | awk 'NR==FNR{w[$1]; next} FNR%4==1{id=substr($1,2)} FNR%4==2{t++; if(id in w){n++; s+=length($0)}} END{printf "%d\t%d\t%d", t, n, s}' $Z/sel/$b.ids -)
  printf '%s\t%s\n' "$b" "$st" >> $Z/sel/summary.tsv
done
echo "선별 정렬 끝 $(date -u +%H:%M) UTC"
cat $Z/sel/summary.tsv
