#!/bin/bash
set -u
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
[ "$F" -lt 25 ] && { echo STOP_DISK; exit 1; }
bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp -o bam   -1 sra/deep_SRR837875_1.fastq.gz -2 sra/deep_SRR837875_2.fastq.gz > logs/deep2.align.log 2>&1
B=$(ls bam/deep_SRR837875_1*bismark*pe.bam 2>/dev/null | head -1)
if [ -n "$B" ]; then
  samtools sort -@ 6 -o bam/deep_egg.sorted.bam "$B" && samtools index bam/deep_egg.sorted.bam && rm -f "$B"
  grep 'Mapping efficiency' logs/deep2.align.log | head -1
  ls -la bam/deep_egg.sorted.bam | awk '{printf "BAM %.2f GB\n",$5/1e9}'
else
  echo DEEP2_FAIL; tail -3 logs/deep2.align.log
fi
rm -f tmp/* 2>/dev/null
df -h / | tail -1
echo DEEP2_DONE
