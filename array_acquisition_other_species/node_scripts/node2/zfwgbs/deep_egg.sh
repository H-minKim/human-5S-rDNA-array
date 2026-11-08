#!/bin/bash
set -u
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
[ "$F" -lt 25 ] && { echo STOP_DISK; exit 1; }
R=SRR837875; N=32000000
[ -f sra/deep_${R}_1.fastq.gz ] || fastq-dump -X $N --split-files --gzip -O sra $R > logs/deep.dl.log 2>&1
mv -f sra/${R}_1.fastq.gz sra/old_${R}_1.fastq.gz 2>/dev/null
mv -f sra/${R}_2.fastq.gz sra/old_${R}_2.fastq.gz 2>/dev/null
ls -la sra/${R}_*.fastq.gz 2>/dev/null | awk '{printf "%s %.2f GB\n",$NF,$5/1e9}'
bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp -o bam   -1 sra/${R}_1.fastq.gz -2 sra/${R}_2.fastq.gz > logs/deep_${R}.align.log 2>&1
B=$(ls bam/${R}_1*bismark*pe.bam 2>/dev/null | head -1)
[ -n "$B" ] && samtools sort -@ 6 -o bam/deep_${R}.sorted.bam "$B" && samtools index bam/deep_${R}.sorted.bam && rm -f "$B"
grep 'Mapping efficiency' logs/deep_${R}.align.log | head -1
rm -f tmp/* 2>/dev/null
df -h / | tail -1
echo DEEP_DONE
