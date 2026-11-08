#!/bin/bash
# 알 시료 깊이 올리기 — 8M 씩 조각내어 순차 정렬 후 합친다
#   32M 을 한 번에 돌리면 bismark 임시파일이 31 GB 까지 불어 디스크가 98% 가 됐다(중단함).
#   8M 실행은 문제없이 끝났으므로 같은 크기로 4 조각을 하나씩 처리하고 매번 임시를 지운다.
#   이미 끝난 1 조각(bam/SRR837875.sorted.bam)은 그대로 쓰고 2~4 조각만 만든다.
set -u
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p chunk logs tmp

guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP — 디스크 여유 ${F}G"; exit 1; fi
}

L=32000000          # 8,000,000 read = 32,000,000 줄
for K in 2 3 4; do
  guard
  [ -f bam/eggC${K}.sorted.bam ] && { echo "SKIP 조각$K"; continue; }
  SKIP=$(( (K-1) * L ))
  for M in 1 2; do
    zcat sra/deep_SRR837875_${M}.fastq.gz | tail -n +$((SKIP+1)) | head -n $L | pigz -p 6 > chunk/c${K}_${M}.fastq.gz
  done
  echo "조각$K 정렬 시작 $(date -u +%H:%M)  $(ls -la chunk/c${K}_1.fastq.gz | awk '{printf \"%.2f GB\", $5/1e9}')"
  bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp -o bam \
    -1 chunk/c${K}_1.fastq.gz -2 chunk/c${K}_2.fastq.gz > logs/eggC${K}.log 2>&1
  B=$(ls bam/c${K}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o bam/eggC${K}.sorted.bam "$B" && samtools index bam/eggC${K}.sorted.bam && rm -f "$B"
    echo "조각$K OK $(grep 'Mapping efficiency' logs/eggC${K}.log | head -1)"
  else
    echo "조각$K FAIL"; tail -3 logs/eggC${K}.log
  fi
  rm -f chunk/c${K}_*.fastq.gz tmp/* 2>/dev/null
  df -h / | tail -1
done

guard
samtools merge -f -@ 6 bam/egg_deep.bam bam/SRR837875.sorted.bam bam/eggC2.sorted.bam bam/eggC3.sorted.bam bam/eggC4.sorted.bam 2>/dev/null \
  && samtools index bam/egg_deep.bam && echo "합본 $(ls -la bam/egg_deep.bam | awk '{printf \"%.2f GB\", $5/1e9}')"
df -h / | tail -1
echo CHUNKS_DONE
