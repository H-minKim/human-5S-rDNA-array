#!/bin/bash
# Zebra finch pilot: the whole path (re-basecall -> align -> modkit -> unit matrix) on one flowcell and N array reads, plus a position check against the original fastq alignment.
#   The check compares the primary alignment of the same read id: original fastq (PAF, sel/<fc>.paf.gz) against the re-basecalled BAM.
#     The maternal and paternal chr2 arrays have identical sequence (same md5), so which one a read lands on is arbitrary; paternal CM109751.1 coordinates are shifted to maternal CM109790.1 (a difference of 61,764 bp, the same for both arrays).
#     The smaller of the two end-point differences (dmin) is recorded, together with whether the shift is a multiple of the 368 bp unit.
#   modkit with the same options as zebrafish (--cpg --combine-strands) plus 5mC and 5hmC; include-bed is BED3 so that no column is read as a strand.
#   Arguments: flowcell, read count, core range.
# MIT License. Copyright (c) 2026 Hyungmin Kim
set -u
Z=~/xsp/zfont
A=~/xsp/asm
MK=~/xsp/tools/modkit
FC="${1:-bTaeGut7_PAW63621_20240814}"
N="${2:-20}"
CORES="${3:-0-15}"
P=$Z/pilot
rm -rf "$P"
mkdir -p "$P/ids"
head -n "$N" $Z/sel/$FC.ids > $P/ids/$FC.ids
echo "$FC" > $P/fc.txt
cut -f1-3 $Z/zfinch_v03_5S_arrays.bed > $P/arrays3.bed
s=$(date +%s)
BC=$P/bc bash ~/xsp/finch_dorado.sh $P/fc.txt "$CORES" $P/ids > $P/dorado_run.log 2>&1
echo "== rebasecalling (fc, rc, read, MM read, seconds)"; cat $P/bc/$FC.done
grep -E "Simplex|Samples/s" $P/bc/$FC.dorado.log | tail -2
echo "rebasecalled BAM primary alignments $(samtools view -c -F 0x904 $P/bc/$FC.bam) / input id $(wc -l < $P/ids/$FC.ids)"
# 1) position check
norm='{c=$2; s=$3; e=$4; if(c=="CM109751.1" && s>96000000 && e<98000000){c="CM109790.1"; s-=61764; e-=61764}; print $1"\t"c"\t"s"\t"e"\t"$5}'
zcat $Z/sel/$FC.paf.gz | awk -F'\t' '$0 ~ /tp:A:P/ {print $1"\t"$6"\t"$8"\t"$9"\t"$5}' | awk -F'\t' "$norm" | sort -k1,1 > $P/orig.tsv
samtools view -F 0x904 $P/bc/$FC.bam | awk -F'\t' '{c=$6; e=$4-1; while(match(c,/[0-9]+[MDN=X]/)){e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; st=(int($2/16)%2)?"-":"+"; print $1"\t"$3"\t"$4-1"\t"e"\t"st}' | awk -F'\t' "$norm" | sort -k1,1 > $P/new.tsv
join -t $'\t' $P/orig.tsv $P/new.tsv | awk -F'\t' '{ds=$3-$7; if(ds<0)ds=-ds; de=$4-$8; if(de<0)de=-de; d=(ds<de)?ds:de; print $1"\t"$2"\t"$6"\t"$5"\t"$9"\t"ds"\t"de"\t"d"\t"sprintf("%.2f", d/368)}' > $P/poscmp.tsv
awk -F'\t' '{n++; if($2==$3)c++; if($4==$5)t++; if($8<=50)a++; else if($8<=400)b++; else x++} END{printf "== position comparison : read %d, same chromosome %d, same strand %d, end difference<=50bp %d, 50-400bp %d, >400bp %d\n", n, c, t, a, b, x}' $P/poscmp.tsv
awk -F'\t' '$8>50' $P/poscmp.tsv | head -10
# ② modkit
$MK pileup $P/bc/$FC.bam $P/$FC.cpg.bed --ref $A/zfinch_v03_dip.fna --cpg --include-bed $P/arrays3.bed --combine-strands --modified-bases 5mC 5hmC -t 4 --log-filepath $P/pileup.log > $P/pileup.out 2>&1
echo "== modkit pileup rc=$?, rows $(wc -l < $P/$FC.cpg.bed 2>/dev/null)"; tail -3 $P/pileup.out
$MK extract full $P/bc/$FC.bam $P/$FC.reads.tsv --ref $A/zfinch_v03_dip.fna --cpg --include-bed $P/arrays3.bed -t 4 --log-filepath $P/extract.log > $P/extract.out 2>&1
echo "== modkit extract rc=$?, rows $(wc -l < $P/$FC.reads.tsv 2>/dev/null)"; tail -3 $P/extract.out
cut -f14 $P/$FC.reads.tsv 2> /dev/null | sort | uniq -c | head
# 3) unit matrix
python3 ~/xsp/units_matrix.py --fa $A/zfinch_v03_dip.fna --bed $Z/zfinch_v03_5S_arrays.bed --tbl $Z/zfinch_v03_mat.tbl,$Z/zfinch_v03_pat.tbl --samples $FC --bam "$P/bc/{s}.bam" --meth "$P/{s}.reads.tsv" --out $P/units
echo "== array summary"; column -t $P/units/5S_array_summary.tsv | head -20
echo "pilot done, $(( $(date +%s) - s )) sec"
