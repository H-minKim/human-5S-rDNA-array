#!/bin/bash
# 금화조 : 5S k-mer 후보 read -> v0.3 이배체(mat+pat) 참조 정렬(PAF) -> 5S array 에 한 번이라도 걸친 read id (flowcell 별)
#   정렬 = minimap2 map-ont --secondary=no · 보조정렬(supplementary)도 걸침으로 인정 · MAPQ 로 거르지 않는다 (탠덤 array)
#   array = zfinch_v03_5S_arrays.bed (cmsearch RF00001 hit 10 kb 이내 묶음 · 3 copy 이상 · mat 4 + pat 4)
#   산출 : sel/<fc>.ids · sel/<fc>.paf.gz · sel/summary.tsv (flowcell · 후보 read · array read · array read 염기합)
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
