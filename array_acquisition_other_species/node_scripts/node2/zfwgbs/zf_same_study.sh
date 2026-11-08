#!/bin/bash
# 같은 연구(PRJNA188516) 안에서 계열을 닫는다 — 대조군 보강
#   지금까지 알(SRR837875)만 이 연구 것이고 체세포 대조(장·심장)는 다른 연구였다.
#   이 연구 하나 안에 : 알 · 정자 · 16세포(ZGA 전) · 1k세포(ZGA 무렵) · Germring(6hpf) · 정소
#   가 다 있으므로 연구·플랫폼 차이를 빼고 비교할 수 있다.
#   앞 작업(zf_more_eggs.sh)이 끝난 뒤 이어서 돈다.
set -u
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p sra bam logs tmp
N=8000000

until grep -q MORE_EGGS_DONE more.out 2>/dev/null; do sleep 60; done
echo "앞 작업 끝 — 대조군 시작 $(date -u +%H:%M)"

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

one SRR837878 16세포_ZGA전
one SRR837895 Germring_6hpf
one SRR837899 정소_조직
echo SAME_STUDY_DONE
