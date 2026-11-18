#!/bin/bash
# Increase the number of chr18 reads.
#   Potok egg SRR800059: 9 reads cover the chr18 array, so 32M spots gives roughly four times as many.
#   Potok muscle SRR800081: 33 on chr18, needed to resolve the chr4 confounding within that same study.
#
#   Aligned in 8M chunks one at a time, deleting the temporary files each time. bismark --parallel copies the whole
#     input into chunks first, so --upto does not reduce the temporary space: 32M in one run reaches 31 GB.
#   Only the needed range is fetched with fastq-dump -N/-X, so the whole run is never stored.
#   Stops if less than 25 GB of disk is free.
set -u
cd ~/zfwgbs || exit 1
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p sra bam logs tmp
CH=8000000

guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP: free disk ${F}G"; exit 1; fi
}

chunk () {
  R=$1; K=$2
  T=bam/${R}_c${K}.sorted.bam
  [ -f "$T" ] && { echo "SKIP $R chunk$K"; return 0; }
  guard
  N0=$(( (K-1) * CH + 1 )); N1=$(( K * CH ))
  fastq-dump -N $N0 -X $N1 --split-files --gzip -O sra "$R" > logs/${R}_c${K}.dl.log 2>&1
  [ -f sra/${R}_2.fastq.gz ] || { echo "DL_FAIL $R chunk$K"; rm -f sra/${R}_*.fastq.gz; return 1; }
  echo "  chunk$K alignment $(date -u +%H:%M)"
  bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp -o bam \
    -1 sra/${R}_1.fastq.gz -2 sra/${R}_2.fastq.gz > logs/${R}_c${K}.align.log 2>&1
  B=$(ls bam/${R}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o "$T" "$B" && samtools index "$T" && rm -f "$B"
    echo "  chunk$K OK $(grep 'Mapping efficiency' logs/${R}_c${K}.align.log | head -1)"
  else
    echo "  chunk$K ALIGN_FAIL"; tail -3 logs/${R}_c${K}.align.log
  fi
  rm -f sra/${R}_*.fastq.gz tmp/* 2>/dev/null
  df -h / | tail -1
}

deep () {
  R=$1; D=$2
  echo "=== $R ($D) deepening $(date -u +%H:%M)"
  # chunk 1 already exists as bam/${R}.sorted.bam
  for K in 2 3 4; do chunk "$R" "$K"; done
  guard
  L="bam/${R}.sorted.bam"
  for K in 2 3 4; do [ -f bam/${R}_c${K}.sorted.bam ] && L="$L bam/${R}_c${K}.sorted.bam"; done
  samtools merge -f -@ 6 bam/${R}_deep.bam $L && samtools index bam/${R}_deep.bam
  echo "merged $R read=$(samtools view -c -@6 bam/${R}_deep.bam)  chunk sum=$(for f in $L; do samtools view -c -@6 $f; done | paste -sd+ | bc)"
  # Once the merge is verified the chunks are deleted to return the space; the original of chunk 1 is kept.
  rm -f bam/${R}_c*.sorted.bam*
  df -h / | tail -1
}

deep SRR800059 Potok_알
deep SRR800081 Potok_근육성체
echo DEEPEN_DONE
