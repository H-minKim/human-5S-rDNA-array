#!/bin/bash
# 제브라피시 WGBS 정렬 — 난자부터
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p bam logs
for S in SRR14703401 SRR14703403; do
  [ -f bam/${S}_done ] && { echo "SKIP $S"; continue; }
  [ -f sra/${S}_1.fastq.gz ] || { echo "NOFQ $S"; continue; }
  bismark --genome ref --parallel 4 -p 2 --bam --temp_dir /tmp \
    -o bam -1 sra/${S}_1.fastq.gz -2 sra/${S}_2.fastq.gz > logs/${S}.align.log 2>&1
  if ls bam/${S}_1_bismark_bt2_pe.bam >/dev/null 2>&1; then
    samtools sort -@ 6 -o bam/${S}.sorted.bam bam/${S}_1_bismark_bt2_pe.bam \
      && samtools index bam/${S}.sorted.bam \
      && rm -f bam/${S}_1_bismark_bt2_pe.bam \
      && touch bam/${S}_done
  fi
  echo "DONE $S  $(date -u +%H:%M)"
  df -h / | tail -1
done
echo ALIGN_ALL_DONE
