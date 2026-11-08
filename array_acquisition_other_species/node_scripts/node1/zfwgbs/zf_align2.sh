#!/bin/bash
# 제브라피시 WGBS 정렬 — 작은 발자국판
#   · --upto 로 앞 800 만 read 쌍만 쓴다(array 평균에는 충분: 전장 약 1.6x, array 800kb 에 read 수천)
#   · --parallel 2 로 임시파일을 줄인다
#   · 임시는 zfwgbs/tmp 에 두고, 디스크 25G 미만이면 스스로 멈춘다
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p bam logs tmp
for S in SRR14703401 SRR14703403; do
  [ -f bam/${S}_done ] && { echo "SKIP $S"; continue; }
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$FREE" -lt 25 ]; then echo "STOP $S — 여유 ${FREE}G"; break; fi
  bismark --genome ref --parallel 2 -p 2 --bam --temp_dir tmp --upto 8000000     -o bam -1 sra/${S}_1.fastq.gz -2 sra/${S}_2.fastq.gz > logs/${S}.align2.log 2>&1
  B=$(ls bam/${S}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o bam/${S}.sorted.bam "$B" && samtools index bam/${S}.sorted.bam && rm -f "$B" && touch bam/${S}_done
    echo "OK $S $(date -u +%H:%M) $(ls -la bam/${S}.sorted.bam | awk '{printf \"%.1f GB\", \/1e9}')"
  else
    echo "FAIL $S"; tail -3 logs/${S}.align2.log
  fi
  rm -f tmp/* 2>/dev/null
  df -h / | tail -1
done
echo ALIGN2_DONE
