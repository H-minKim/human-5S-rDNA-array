#!/bin/bash
# 제브라피시 5S — 난자(알) → 배아 → 체세포 계열 WGBS 정렬
#   묻는 것 : 체세포(뇌 ONT)에서 chr18(체세포형)이 제 flank 보다 낮았다.
#             알에서는 chr4(모계형)로 방향이 뒤집히는가?
#   선별 검사를 통과한 run 만 쓴다(변환 read >= 88%).
#   한 시료씩 받고·정렬하고·임시 지우고 다음으로 간다. 디스크 25 GB 미만이면 멈춘다.
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
