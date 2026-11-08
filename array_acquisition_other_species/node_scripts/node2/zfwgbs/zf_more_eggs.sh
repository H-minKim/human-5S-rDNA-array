#!/bin/bash
# 알 시료 보강 — 알 재현본과 정자·초기배아를 같은 연구(PRJNA188516)에서 채운다
#   지금까지 알이 SRR837875 하나뿐이었다(두 번째 알 ERR10788129 는 매핑률 2.8% 로 실패).
#   PRJNA188516 에는 우리가 성공한 알과 같은 연구·같은 플랫폼으로
#   알 rep2 · 정자 · 1k세포 배아가 있다. 이것들을 채우면 계열이 한 연구 안에서 닫힌다.
#   변환 검사는 통과했어도 정렬이 실패할 수 있다(ERR10788129 전례) → 매핑률을 반드시 확인한다.
#   8M spot 씩 하나하나 · 매번 임시 삭제 · 디스크 25 GB 미만이면 중단.
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
