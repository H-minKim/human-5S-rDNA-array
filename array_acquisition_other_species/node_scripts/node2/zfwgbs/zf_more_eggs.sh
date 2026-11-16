#!/bin/bash
# Strengthening the egg samples: an egg replicate, sperm and early embryos are taken from the same study, PRJNA188516.
#   PRJNA188516 holds an egg replicate, sperm and a 1k-cell embryo from the same study and the same platform
#   as the egg sample SRR837875, which closes the series within one study.
#   Passing the conversion check does not guarantee the alignment will succeed (ERR10788129 mapped at only 2.8%), so the mapping rate must be checked.
#   8M spots at a time, with the temporary files deleted after each, stopping below 25 GB of free disk.
set -u
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p sra bam logs tmp
N=8000000

guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP — 디스크 여유 ${F}G"; exit 1; fi
}

one () {
  R=$1; D=$2
  [ -f bam/${R}_done ] && { echo "SKIP $R ($D)"; return; }
  guard
  [ -f sra/${R}_1.fastq.gz ] || fastq-dump -X $N --split-files --gzip -O sra "$R" > logs/${R}.dl.log 2>&1
  [ -f sra/${R}_1.fastq.gz ] || { echo "DL_FAIL $R"; return; }
  echo "정렬 시작 $R ($D) $(date -u +%H:%M)"
  bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp -o bam \
    -1 sra/${R}_1.fastq.gz -2 sra/${R}_2.fastq.gz > logs/${R}.align.log 2>&1
  B=$(ls bam/${R}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o bam/${R}.sorted.bam "$B" && samtools index bam/${R}.sorted.bam && rm -f "$B" && touch bam/${R}_done
    echo "OK $R  $(grep 'Mapping efficiency' logs/${R}.align.log | head -1)"
  else
    echo "ALIGN_FAIL $R"; tail -3 logs/${R}.align.log
  fi
  rm -f sra/${R}_*.fastq.gz tmp/* 2>/dev/null
  df -h / | tail -1
}

one SRR837876   알rep2_같은연구
one SRR837872   정자_같은연구
one ERR10788127 알_bouncerKO
one SRR837894   1k세포_같은연구
echo MORE_EGGS_DONE
