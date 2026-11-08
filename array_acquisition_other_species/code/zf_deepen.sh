#!/bin/bash
# chr18 read 수를 올린다 — 지금 결론의 유일한 약점이 read 수다.
#   Potok 알   SRR800059 : chr18 array 를 덮은 read 9 개  → 32M 으로 약 4 배
#   Potok 근육 SRR800081 : chr18 33 개 · chr4 교락을 그 연구 안에서 풀기 위해 필요
#
#   8M 씩 조각내어 하나씩 정렬하고 매번 임시 삭제(bismark --parallel 은 입력 전체를
#     먼저 조각 복사하므로 --upto 로는 임시가 안 줄어든다 — 32M 한 번에 = 31GB).
#   fastq-dump -N/-X 로 구간만 받아 전체를 저장하지 않는다.
#   디스크 25GB 미만이면 중단.
set -u
cd ~/zfwgbs || exit 1
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p sra bam logs tmp
CH=8000000

guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP — 디스크 여유 ${F}G"; exit 1; fi
}

chunk () {
  R=$1; K=$2
  T=bam/${R}_c${K}.sorted.bam
  [ -f "$T" ] && { echo "SKIP $R 조각$K"; return 0; }
  guard
  N0=$(( (K-1) * CH + 1 )); N1=$(( K * CH ))
  fastq-dump -N $N0 -X $N1 --split-files --gzip -O sra "$R" > logs/${R}_c${K}.dl.log 2>&1
  [ -f sra/${R}_2.fastq.gz ] || { echo "DL_FAIL $R 조각$K"; rm -f sra/${R}_*.fastq.gz; return 1; }
  echo "  조각$K 정렬 $(date -u +%H:%M)"
  bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp -o bam \
    -1 sra/${R}_1.fastq.gz -2 sra/${R}_2.fastq.gz > logs/${R}_c${K}.align.log 2>&1
  B=$(ls bam/${R}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o "$T" "$B" && samtools index "$T" && rm -f "$B"
    echo "  조각$K OK $(grep 'Mapping efficiency' logs/${R}_c${K}.align.log | head -1)"
  else
    echo "  조각$K ALIGN_FAIL"; tail -3 logs/${R}_c${K}.align.log
  fi
  rm -f sra/${R}_*.fastq.gz tmp/* 2>/dev/null
  df -h / | tail -1
}

deep () {
  R=$1; D=$2
  echo "=== $R ($D) 깊이 올리기 $(date -u +%H:%M)"
  # 1 조각은 이미 bam/${R}.sorted.bam 으로 있다
  for K in 2 3 4; do chunk "$R" "$K"; done
  guard
  L="bam/${R}.sorted.bam"
  for K in 2 3 4; do [ -f bam/${R}_c${K}.sorted.bam ] && L="$L bam/${R}_c${K}.sorted.bam"; done
  samtools merge -f -@ 6 bam/${R}_deep.bam $L && samtools index bam/${R}_deep.bam
  echo "합본 $R read=$(samtools view -c -@6 bam/${R}_deep.bam)  조각합=$(for f in $L; do samtools view -c -@6 $f; done | paste -sd+ | bc)"
  # 합본이 확인되면 조각은 지워 공간을 돌려준다(1 조각 원본은 남긴다)
  rm -f bam/${R}_c*.sorted.bam*
  df -h / | tail -1
}

deep SRR800059 Potok_알
deep SRR800081 Potok_근육성체
echo DEEPEN_DONE
