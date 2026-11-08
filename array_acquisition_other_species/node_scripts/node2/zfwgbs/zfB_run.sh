#!/bin/bash
# node2: builds the bisulfite index and aligns the egg WGBS
#   입력 fastq 를 미리 800 만 spot 으로 줄여 놓았으므로 bismark 임시파일도 작다.
# Checks the free disk at every step and stops below 25 GB, to protect other data.
set -u
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p logs bam tmp

guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP — 디스크 여유 ${F}G"; exit 1; fi
}

guard
if [ ! -d ref/Bisulfite_Genome ]; then
  echo "색인 시작 $(date -u +%H:%M)"
  bismark_genome_preparation --parallel 6 --bowtie2 ref > logs/prep.log 2>&1
  echo "색인 끝 $(date -u +%H:%M) exit=$?"
fi
ls ref/Bisulfite_Genome || { echo PREP_FAIL; exit 1; }
df -h / | tail -1

for S in SRR14703401 SRR14703403; do
  guard
  [ -f bam/${S}_done ] && { echo "SKIP $S"; continue; }
  echo "정렬 시작 $S $(date -u +%H:%M)"
  bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp \
    -o bam -1 sra/${S}_1.fastq.gz -2 sra/${S}_2.fastq.gz > logs/${S}.align.log 2>&1
  B=$(ls bam/${S}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o bam/${S}.sorted.bam "$B" && samtools index bam/${S}.sorted.bam && rm -f "$B" && touch bam/${S}_done
    echo "OK $S $(date -u +%H:%M)"
    grep -E 'Mapping efficiency|C methylated in CpG' logs/${S}.align.log | head -3
  else
    echo "FAIL $S"; tail -4 logs/${S}.align.log
  fi
  rm -f tmp/* 2>/dev/null
  df -h / | tail -1
done
echo RUN_DONE
