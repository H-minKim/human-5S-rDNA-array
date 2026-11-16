#!/bin/bash
# Arabidopsis, final step: re-called BAM per machine -> merge -> modkit (CG, CHG, CHH) -> read x unit matrix -> position check against the original SUP alignment
#   BAMs produced on the other nodes are moved into the bc folder beforehand; the Col-CC reference md5 427863e6c8e1 was confirmed identical on all three nodes.
#   modkit: all contexts, since this is a plant (--motif CG 0 --motif CHG 0 --motif CHH 0). CHH is not symmetric, so --combine-strands is not used. 5mC and 5hmC.
#   Unit matrix: units_matrix.py. Context is decided from the reference, and because the input is streamed it must be sorted by read_id then ref_position.
#   Position check: the re-basecalled primary alignment against whichever of that read's original SUP alignments (all PAF rows) overlaps it, comparing start points. The unit is about 500 bp.
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
H=~/xsp/ath
A=~/xsp/asm
MK=~/xsp/tools/modkit
F=$H/final
S=ColCC_ATRun3
mkdir -p $F
cut -f1-3 $H/arabidopsis_ColCC_5S_arrays.bed > $F/arrays3.bed
ls $H/bc/*.bam > $F/bams.txt
samtools merge -f -@4 -b $F/bams.txt -o $F/$S.all.bam && samtools index $F/$S.all.bam
samtools view -b -L $F/arrays3.bed -o $F/$S.5S.bam $F/$S.all.bam && samtools index $F/$S.5S.bam
$MK pileup $F/$S.5S.bam $F/$S.5S.allC.bed --ref $A/arabidopsis_ColCC.fna --motif CG 0 --motif CHG 0 --motif CHH 0 --include-bed $F/arrays3.bed --modified-bases 5mC 5hmC -t 8 --log-filepath $F/pileup.log > /dev/null 2>&1
echo "pileup rc=$? · 행 $(wc -l < $F/$S.5S.allC.bed)"
$MK extract full $F/$S.5S.bam $F/$S.5S.reads.tsv --ref $A/arabidopsis_ColCC.fna --include-bed $F/arrays3.bed -t 8 --log-filepath $F/extract.log > /dev/null 2>&1
echo "extract rc=$? · 행 $(wc -l < $F/$S.5S.reads.tsv)"
x=$F/$S.5S.reads.tsv
c=$(head -1 $x | tr '\t' '\n' | grep -n -x ref_position | cut -d: -f1)
{ head -1 $x; tail -n +2 $x | LC_ALL=C sort -t $'\t' -k1,1 -k${c},${c}n -S 4G --parallel 4 -T $F; } > $x.sorted && mv $x.sorted $x
python3 ~/xsp/units_matrix.py --fa $A/arabidopsis_ColCC.fna --bed $H/arabidopsis_ColCC_5S_arrays.bed --tbl ~/xsp/cms/arabidopsis_ColCC.tbl \
  --samples $S --bam "$F/{s}.5S.bam" --meth "$F/{s}.5S.reads.tsv" --out $F/units
echo "== BAM read $(samtools view -c $F/$S.5S.bam) · MM $(samtools view $F/$S.5S.bam | grep -c MM:Z)"
echo "== pileup 맥락·변형별 행 (이름 열)"; cut -f4 $F/$S.5S.allC.bed | sort | uniq -c
# position check
zcat $H/sup.colcc.paf.gz | awk -F'\t' '{print $1"\t"$6"\t"$8"\t"$9}' > $F/orig_aln.tsv
samtools view -F 0x904 $F/$S.all.bam | awk -F'\t' '{c=$6; e=$4-1; while(match(c,/[0-9]+[MDN=X]/)){e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; print $1"\t"$3"\t"$4-1"\t"e}' > $F/new_aln.tsv
awk -F'\t' 'NR==FNR{o[$1]=o[$1] $2 "," $3 "," $4 ";"; next}
  {best=-1; if($1 in o){n=split(o[$1],L,";"); for(i=1;i<n;i++){split(L[i],q,","); if(q[1]==$2 && q[2]<$4 && q[3]>$3){d=q[2]-$3; if(d<0)d=-d; if(best<0||d<best)best=d}}}
   if(!($1 in o)) cls="원정렬없음"; else if(best<0) cls="겹침없음"; else if(best<=50) cls="0-50bp"; else if(best<=500) cls="51-500bp(1유닛안)"; else cls=">500bp";
   print $1"\t"$2"\t"$3"\t"best"\t"cls}' $F/orig_aln.tsv $F/new_aln.tsv > $F/poscmp.tsv
echo "== 위치 대조 (재콜링 primary vs 원 SUP 정렬 중 겹치는 것)"; cut -f5 $F/poscmp.tsv | sort | uniq -c
echo "== 재콜링 조각 (이름표 · rc · read · MM · 초)"; cat $H/bc/*.done
echo "== array 요약"; column -t $F/units/5S_array_summary.tsv
echo "마무리 끝 $(date -u +%H:%M) UTC"
