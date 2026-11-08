#!/bin/bash
# 제브라피시 WGBS 추가 시료 — 하나씩 받고 즉시 압축(디스크 아끼기)
cd ~/zfwgbs/sra
for R in SRR19364705 SRR11614917 SRR12797061; do
  [ -f ${R}.done ] && { echo "SKIP $R"; continue; }
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  [ "$FREE" -lt 25 ] && { echo "STOP $R — 디스크 여유 ${FREE}G"; break; }
  fasterq-dump --split-files --threads 6 -O . $R > ../logs/${R}.dl.log 2>&1
  if ls ${R}_1.fastq >/dev/null 2>&1; then
    pigz -p 8 ${R}_*.fastq && touch ${R}.done
    echo "OK $R  $(ls -la ${R}*.fastq.gz | awk '{s+=$5} END{printf \"%.1f GB\", s/1e9}')"
  else
    echo "FAIL $R"; tail -2 ../logs/${R}.dl.log
  fi
  df -h / | tail -1
done
echo DL2_DONE
