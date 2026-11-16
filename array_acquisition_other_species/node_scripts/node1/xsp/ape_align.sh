#!/bin/bash
# Great apes: SAM records selected by 5S k-mer (MM/ML kept) -> fastq -> alignment to the species assembly -> array read BAM -> modkit -> read x unit matrix.
#   The alignment is the same as for zebrafish, zebra finch and Arabidopsis: samtools fastq -T MM,ML | minimap2 -y -ax map-ont --secondary=no, primary only.
#   No MAPQ filter: mapq 0 is normal in a tandem array.
#   The ONT BAM is R9.4.1 guppy 6.3.8 (5mC CpG model), so modkit is run for 5mC CpG.
#   Arguments: species key (e.g. chimp_mPanTro3), the source BAM URL to take the header from, the selection folder, the output folder.
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
SP="${1:?종키}"
HDR="${2:?머리말 URL}"
HITD="${3:?선별폴더}"
OUT="${4:?출력폴더}"
if [ "$HDR" = auto ]; then HDR="https://genomeark.s3.amazonaws.com/$(head -1 $HITD/bams.txt)"; fi
A=~/xsp/asm/$SP.fna
MK=~/xsp/tools/modkit
mkdir -p "$OUT"
awk -F'\t' -v sp="$SP" 'NR>1 && $1==sp {printf "%s\t%d\t%s\t%s_c%s\n", $2, $3-1, $4, sp, $6}' ~/xsp/struct/ALL.struct.tsv > $OUT/$SP.arrays.bed
cut -f1-3 $OUT/$SP.arrays.bed > $OUT/arrays3.bed
echo "array $(wc -l < $OUT/$SP.arrays.bed) 개"
samtools view -H "$HDR" > $OUT/hdr.sam 2> /dev/null
# Some flowcells were basecalled twice (siamang PAI18731_dorado-0.9.1 and _guppy-6.3.8).
#   -> one version per flowcell (dorado > guppy-6.3.8 > guppy-6.3.7), then deduplicated again by read id.
ls $HITD/*.hit.sam 2> /dev/null | awk -F/ '{f=$NF; split(f, a, "_"); fc=a[1]; p=1; if(f ~ /dorado/) p=3; else if(f ~ /6[.]3[.]8/) p=2; if(p > pr[fc]) {pr[fc]=p; pick[fc]=$0}} END{for(k in pick) print pick[k]}' | sort > $OUT/use.txt
echo "flowcell $(ls $HITD/*.hit.sam 2> /dev/null | wc -l) 개 중 쓰는 판본 $(wc -l < $OUT/use.txt) 개"
n=$(cat $(cat $OUT/use.txt) | wc -l)
u=$(cat $(cat $OUT/use.txt) | awk '!s[$1]++' | wc -l)
echo "중복 제거 : 기록 $n -> read $u"
echo "골라낸 기록 $n"
cat $(cat $OUT/use.txt) | awk '!s[$1]++' | cat $OUT/hdr.sam - | samtools fastq -T MM,ML - 2> /dev/null \
  | minimap2 -y -ax map-ont -t 12 --secondary=no "$A" - 2> $OUT/mm2.log \
  | samtools view -b -F 0x900 - | samtools sort -@4 -o $OUT/$SP.all.bam - && samtools index $OUT/$SP.all.bam
samtools view -b -L $OUT/arrays3.bed -o $OUT/$SP.5S.bam $OUT/$SP.all.bam && samtools index $OUT/$SP.5S.bam
echo "정렬 $(samtools view -c $OUT/$SP.all.bam) · array 에 걸린 것 $(samtools view -c $OUT/$SP.5S.bam) · MM $(samtools view $OUT/$SP.5S.bam | grep -c MM:Z)"
$MK pileup $OUT/$SP.5S.bam $OUT/$SP.5S.cpg.bed --ref "$A" --cpg --include-bed $OUT/arrays3.bed --combine-strands --modified-bases 5mC -t 8 --log-filepath $OUT/pileup.log > /dev/null 2>&1
echo "pileup rc=$? 행 $(wc -l < $OUT/$SP.5S.cpg.bed)"
rm -f $OUT/$SP.5S.reads.tsv
$MK extract full $OUT/$SP.5S.bam $OUT/$SP.5S.reads.tsv --ref "$A" --cpg --include-bed $OUT/arrays3.bed -t 8 --log-filepath $OUT/extract.log > /dev/null 2>&1
echo "extract rc=$? 행 $(wc -l < $OUT/$SP.5S.reads.tsv)"
x=$OUT/$SP.5S.reads.tsv
c=$(head -1 $x | tr '\t' '\n' | grep -n -x ref_position | cut -d: -f1)
{ head -1 $x; tail -n +2 $x | LC_ALL=C sort -t $'\t' -k1,1 -k${c},${c}n -S 2G --parallel 4 -T $OUT; } > $x.sorted && mv $x.sorted $x
python3 ~/xsp/units_matrix.py --fa "$A" --bed $OUT/$SP.arrays.bed --tbl ~/xsp/cms/$SP.tbl \
  --samples "$SP" --bam "$OUT/{s}.5S.bam" --meth "$OUT/{s}.5S.reads.tsv" --out $OUT/units
echo "== array 요약"
column -t $OUT/units/5S_array_summary.tsv
echo "끝 $(date -u +%H:%M) UTC"
