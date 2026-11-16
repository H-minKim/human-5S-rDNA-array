#!/bin/bash
# Zebrafish 5S: WGBS alignment of the oocyte, embryo and somatic series.
#   The question: in the soma, the brain ONT data, chr18, the somatic type, was lower than its own flank.
#             Does the direction invert to chr4, the maternal type, in the egg?
#   Only runs that passed the screening check are used, at 88% or more converted reads.
#   One sample is downloaded, aligned and its temporary files deleted before moving to the next. The run stops below 25 GB of free disk.
set -u
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p sra bam logs tmp
N=8000000

guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP — 디스크 여유 ${F}G"; exit 1; fi
}

run_one () {
  R=$1; D=$2
  [ -f bam/${R}_done ] && { echo "SKIP $R ($D)"; return; }
  guard
  if [ ! -f sra/${R}_1.fastq.gz ]; then
    fastq-dump -X $N --split-files --gzip -O sra "$R" > logs/${R}.dl.log 2>&1
  fi
  [ -f sra/${R}_1.fastq.gz ] || { echo "DL_FAIL $R"; return; }
  echo "정렬 시작 $R ($D) $(date -u +%H:%M)"
  bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp \
    -o bam -1 sra/${R}_1.fastq.gz -2 sra/${R}_2.fastq.gz > logs/${R}.align.log 2>&1
  B=$(ls bam/${R}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o bam/${R}.sorted.bam "$B" && samtools index bam/${R}.sorted.bam \
      && rm -f "$B" && touch bam/${R}_done
    echo "OK $R $(date -u +%H:%M)  $(grep 'Mapping efficiency' logs/${R}.align.log | head -1)"
  else
    echo "ALIGN_FAIL $R"; tail -3 logs/${R}.align.log
  fi
  rm -f tmp/* 2>/dev/null
  df -h / | tail -1
}

run_one SRR837875    알_Jiang2013
run_one ERR10788129  알_PRJEB59047
run_one SRR11614917  배아24hpf
run_one SRR12797061  장_체세포
run_one SRR12797059  심장_체세포
echo SERIES_DONE
