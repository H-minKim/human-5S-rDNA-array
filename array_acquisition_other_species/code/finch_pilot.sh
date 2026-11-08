#!/bin/bash
# 금화조 파일럿 : flowcell 1 개 · array read N 개로 전 과정 (재콜링 -> 정렬 -> modkit -> 유닛행렬) + 원 fastq 정렬과 위치 대조
#   대조 : 같은 read id 의 primary 정렬 — 원 fastq(PAF, sel/<fc>.paf.gz) vs 재콜링(BAM)
#     mat·pat chr2 array 는 서열이 같아(md5 동일) 어느 쪽에 붙을지 임의 -> pat CM109751.1 좌표를 mat CM109790.1 로 옮겨 비교 (차 61,764 bp, 두 array 모두 같은 차)
#     끝 위치 둘 중 가까운 쪽 차이(dmin) · 368 bp(유닛) 배수로 밀렸는지 표시
#   modkit : 제브라피시와 같은 옵션 (--cpg --combine-strands) + 5mC·5hmC · include-bed 는 BED3 (가닥 열 해석 방지)
#   인자 : flowcell · read 수 · 코어 범위
# MIT License · Copyright (c) 2026 Hungmin Kim
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
echo "== 재콜링 (fc · rc · read · MM read · 초)"; cat $P/bc/$FC.done
grep -E "Simplex|Samples/s" $P/bc/$FC.dorado.log | tail -2
echo "재콜링 BAM primary 정렬 $(samtools view -c -F 0x904 $P/bc/$FC.bam) / 입력 id $(wc -l < $P/ids/$FC.ids)"
# ① 위치 대조
norm='{c=$2; s=$3; e=$4; if(c=="CM109751.1" && s>96000000 && e<98000000){c="CM109790.1"; s-=61764; e-=61764}; print $1"\t"c"\t"s"\t"e"\t"$5}'
zcat $Z/sel/$FC.paf.gz | awk -F'\t' '$0 ~ /tp:A:P/ {print $1"\t"$6"\t"$8"\t"$9"\t"$5}' | awk -F'\t' "$norm" | sort -k1,1 > $P/orig.tsv
samtools view -F 0x904 $P/bc/$FC.bam | awk -F'\t' '{c=$6; e=$4-1; while(match(c,/[0-9]+[MDN=X]/)){e+=substr(c,RSTART,RLENGTH-1); c=substr(c,RSTART+RLENGTH)}; st=(int($2/16)%2)?"-":"+"; print $1"\t"$3"\t"$4-1"\t"e"\t"st}' | awk -F'\t' "$norm" | sort -k1,1 > $P/new.tsv
join -t $'\t' $P/orig.tsv $P/new.tsv | awk -F'\t' '{ds=$3-$7; if(ds<0)ds=-ds; de=$4-$8; if(de<0)de=-de; d=(ds<de)?ds:de; print $1"\t"$2"\t"$6"\t"$5"\t"$9"\t"ds"\t"de"\t"d"\t"sprintf("%.2f", d/368)}' > $P/poscmp.tsv
awk -F'\t' '{n++; if($2==$3)c++; if($4==$5)t++; if($8<=50)a++; else if($8<=400)b++; else x++} END{printf "== 위치 대조 : read %d · 같은 염색체 %d · 같은 가닥 %d · 끝 차이<=50bp %d · 50-400bp %d · >400bp %d\n", n, c, t, a, b, x}' $P/poscmp.tsv
awk -F'\t' '$8>50' $P/poscmp.tsv | head -10
# ② modkit
$MK pileup $P/bc/$FC.bam $P/$FC.cpg.bed --ref $A/zfinch_v03_dip.fna --cpg --include-bed $P/arrays3.bed --combine-strands --modified-bases 5mC 5hmC -t 4 --log-filepath $P/pileup.log > $P/pileup.out 2>&1
echo "== modkit pileup rc=$? · 행 $(wc -l < $P/$FC.cpg.bed 2>/dev/null)"; tail -3 $P/pileup.out
$MK extract full $P/bc/$FC.bam $P/$FC.reads.tsv --ref $A/zfinch_v03_dip.fna --cpg --include-bed $P/arrays3.bed -t 4 --log-filepath $P/extract.log > $P/extract.out 2>&1
echo "== modkit extract rc=$? · 행 $(wc -l < $P/$FC.reads.tsv 2>/dev/null)"; tail -3 $P/extract.out
cut -f14 $P/$FC.reads.tsv 2> /dev/null | sort | uniq -c | head
# ③ 유닛 행렬
python3 ~/xsp/units_matrix.py --fa $A/zfinch_v03_dip.fna --bed $Z/zfinch_v03_5S_arrays.bed --tbl $Z/zfinch_v03_mat.tbl,$Z/zfinch_v03_pat.tbl --samples $FC --bam "$P/bc/{s}.bam" --meth "$P/{s}.reads.tsv" --out $P/units
echo "== array 요약"; column -t $P/units/5S_array_summary.tsv | head -20
echo "파일럿 끝 · $(( $(date +%s) - s )) 초"
