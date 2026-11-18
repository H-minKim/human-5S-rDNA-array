#!/bin/bash
# Zebra finch, final step: re-called flowcell BAM -> merge per library -> modkit -> read x unit matrix -> position check against the original alignment
#   Library types are not mixed: LIG = ligation (PAQ25082), UL = ultralong (eight R10 5 kHz flowcells).
# The BAM for this share must be moved to ~/xsp/zfont/bc beforehand
#   modkit with the same options as zebrafish and the pilot (--cpg --combine-strands, 5mC and 5hmC); include-bed is BED3.
#   Position check: how far the start of the re-basecalled primary alignment is from whichever of that read's original fastq alignments (PAF, primary and supplementary) overlaps it.
#     The maternal and paternal chr2 arrays are identical in sequence, so paternal CM109751.1 coordinates are shifted into maternal ones for the comparison (-61,764).
# MIT License. Copyright (c) 2026 Hyungmin Kim
set -u
Z=~/xsp/zfont
A=~/xsp/asm
MK=~/xsp/tools/modkit
F=$Z/final
mkdir -p $F
cut -f1-3 $Z/zfinch_v03_5S_arrays.bed > $F/arrays3.bed
ls $Z/bc/bTaeGut7_PAQ25082_*.bam > $F/LIG.bams
ls $Z/bc/*.bam | grep -v PAQ25082 > $F/UL.bams
for g in LIG UL; do
  samtools merge -f -@4 -b $F/$g.bams -o $F/bTaeGut7_$g.all.bam && samtools index $F/bTaeGut7_$g.all.bam
  samtools view -b -L $F/arrays3.bed -o $F/bTaeGut7_$g.5S.bam $F/bTaeGut7_$g.all.bam && samtools index $F/bTaeGut7_$g.5S.bam
  $MK pileup $F/bTaeGut7_$g.5S.bam $F/bTaeGut7_$g.5S.cpg.bed --ref $A/zfinch_v03_dip.fna --cpg --include-bed $F/arrays3.bed --combine-strands --modified-bases 5mC 5hmC -t 8 --log-filepath $F/$g.pileup.log > /dev/null 2>&1
  $MK extract full $F/bTaeGut7_$g.5S.bam $F/bTaeGut7_$g.5S.reads.tsv --ref $A/zfinch_v03_dip.fna --cpg --include-bed $F/arrays3.bed -t 8 --log-filepath $F/$g.extract.log > /dev/null 2>&1
  echo "$g : BAM read $(samtools view -c $F/bTaeGut7_$g.5S.bam), MM $(samtools view $F/bTaeGut7_$g.5S.bam | grep -c MM:Z), pileup rows $(wc -l < $F/bTaeGut7_$g.5S.cpg.bed), extract rows $(wc -l < $F/bTaeGut7_$g.5S.reads.tsv)"
  # units_matrix.py streams its input, so it must be sorted by read_id then ref_position; the column indices are located from the header.
  x=$F/bTaeGut7_$g.5S.reads.tsv
  c=$(head -1 $x | tr '\t' '\n' | grep -n -x ref_position | cut -d: -f1)
  { head -1 $x; tail -n +2 $x | LC_ALL=C sort -t $'\t' -k1,1 -k${c},${c}n -S 2G --parallel 4; } > $x.sorted && mv $x.sorted $x
done
python3 ~/xsp/units_matrix.py --fa $A/zfinch_v03_dip.fna --bed $Z/zfinch_v03_5S_arrays.bed --tbl $Z/zfinch_v03_mat.tbl,$Z/zfinch_v03_pat.tbl \
  --samples bTaeGut7_LIG,bTaeGut7_UL --bam "$F/{s}.5S.bam" --meth "$F/{s}.5S.reads.tsv" --out $F/units
# position check
norm='{c=$2; s=$3; e=$4; if(c=="CM109751.1" && s>96000000 && e<98000000){c="CM109790.1"; s-=61764; e-=61764}; print $1"\t"c"\t"s"\t"e"\t"$5}'
zcat $Z/sel/*.paf.gz | awk -F'\t' '{print $1"\t"$6"\t"$8"\t"$9"\t"$5}' | awk -F'\t' "$norm" > $F/orig_aln.tsv
for b in $Z/bc/*.bam; do samtools view -F 0x904 $b; done | awk -F'\t' '{c=$6; e=$4-1; while(match(c,/[0-9]+[MDN=X]/)){e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; st=(int($2/16)%2)?"-":"+"; print $1"\t"$3"\t"$4-1"\t"e"\t"st}' | awk -F'\t' "$norm" > $F/new_aln.tsv
awk -F'\t' 'NR==FNR{k=$1; o[k]=o[k] $2 "," $3 "," $4 ";"; next}
  {best=-1; if($1 in o){n=split(o[$1],L,";"); for(i=1;i<n;i++){split(L[i],q,","); if(q[1]==$2 && q[2]<$4 && q[3]>$3){d=q[2]-$3; if(d<0)d=-d; if(best<0||d<best)best=d}}}
   if(!($1 in o)) cls="원정렬없음"; else if(best<0) cls="겹침없음"; else if(best<=50) cls="0-50bp"; else if(best<=368) cls="51-368bp(1유닛안)"; else cls=">368bp";
   print $1"\t"$2"\t"$3"\t"best"\t"cls}' $F/orig_aln.tsv $F/new_aln.tsv > $F/poscmp.tsv
echo "== position comparison (rebasecalled primary vs original alignment, overlapping ones)"; cut -f5 $F/poscmp.tsv | sort | uniq -c
echo "== flowcell (fc, rc, read, MM, seconds)"; cat $Z/bc/*.done
echo "== array summary"; column -t $F/units/5S_array_summary.tsv
echo "finalize finished $(date -u +%H:%M) UTC"
