#!/bin/bash
# Zebrafish 5S array reads and methylation extraction (ENA PRJEB53713, 6 control individuals, adult brain, R10.4.1, Guppy 6.4.2, 5mC_CG hac).
#   1) Only the GRCz11 5S regions are extracted remotely from the indexed ENA BAMs aligned to GRCz11, so the full 458 GB is never downloaded.
#      Primary only (-F 0x900): secondary records carry no MM/ML.
#   2) Realigned to GRCz12tu, the gapless T2T assembly, keeping MM/ML: samtools fastq -T MM,ML | minimap2 -y -ax map-ont.
#   3) Only reads overlapping a GRCz12tu 5S array are kept, giving one BAM per individual.
#   4) modkit pileup over the array regions for CpG, and modkit extract per read.
#   Arguments: the GRCz11 5S region BED including padding, and the GRCz12tu 5S array BED.
set -u
X=$HOME/xsp
Z=$X/zf
REG11="${1:?GRCz11 영역 BED}"
ARR12="${2:?GRCz12tu array BED}"
REF=$X/asm/zebrafish_GRCz12tu.fna
MK=$X/tools/modkit
mkdir -p $Z/bam $Z/meth $Z/log
for r in "ERR12118112 Control1" "ERR12118113 Control2" "ERR12118114 Control3" "ERR12118115 Control4" "ERR12118116 Control5" "ERR12118117 Control6"; do
  set -- $r
  acc=$1; nm=$2
  U=https://ftp.sra.ebi.ac.uk/vol1/run/${acc:0:6}/$acc/$nm.bam
  [ -s $Z/bam/$nm.g12.5S.bam ] && { echo "$nm 이미 있음"; continue; }
  echo "[$(date -u +%H:%M)] $nm 원격 영역 추출 시작" >> $Z/log/progress.log
  samtools view -b -F 0x900 -M -L "$REG11" -o $Z/bam/$nm.g11.5S.bam "$U" 2> $Z/log/$nm.view.err
  n11=$(samtools view -c $Z/bam/$nm.g11.5S.bam)
  mm11=$(samtools view $Z/bam/$nm.g11.5S.bam | grep -c "MM:Z")
  samtools fastq -T MM,ML $Z/bam/$nm.g11.5S.bam 2>/dev/null \
    | minimap2 -y -ax map-ont -t 6 --secondary=no "$REF" - 2> $Z/log/$nm.mm2.err \
    | samtools view -b -F 0x900 - \
    | samtools sort -@2 -o $Z/bam/$nm.g12.all.bam - && samtools index $Z/bam/$nm.g12.all.bam
  samtools view -b -L "$ARR12" -o $Z/bam/$nm.g12.5S.bam $Z/bam/$nm.g12.all.bam && samtools index $Z/bam/$nm.g12.5S.bam
  n12=$(samtools view -c $Z/bam/$nm.g12.5S.bam)
  $MK pileup $Z/bam/$nm.g12.5S.bam $Z/meth/$nm.5S.cpg.bed --ref "$REF" --cpg --include-bed "$ARR12" --combine-strands -t 4 --log-filepath $Z/log/$nm.pileup.log > /dev/null 2>&1
  $MK extract full $Z/bam/$nm.g12.5S.bam $Z/meth/$nm.5S.reads.tsv --ref "$REF" --cpg --include-bed "$ARR12" -t 4 --log-filepath $Z/log/$nm.extract.log > /dev/null 2>&1
  echo -e "$nm\t$acc\tGRCz11영역read $n11\tMM태그 $mm11\tGRCz12tu 5S array read $n12\tCpG행 $(wc -l < $Z/meth/$nm.5S.cpg.bed 2>/dev/null)" >> $Z/log/summary.tsv
  echo "[$(date -u +%H:%M)] $nm 끝 · GRCz11 영역 read $n11 · 5S array read $n12" >> $Z/log/progress.log
done
touch $Z/log/ZF_DONE
