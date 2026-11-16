#!/bin/bash
cd ~/zfwgbs
# 1) compress the fastq to free disk; bismark reads .gz
for f in sra/*.fastq; do [ -f "$f" ] && pigz -p 8 "$f" 2>/dev/null || gzip "$f"; done
df -h / | tail -1
# 2) bismark accepts only .fa and .fasta, so the file is renamed
[ -f ref/zebrafish_GRCz12tu.fna ] && mv ref/zebrafish_GRCz12tu.fna ref/zebrafish_GRCz12tu.fa
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
bismark_genome_preparation --parallel 6 --bowtie2 ref > logs/prep2.log 2>&1
echo PREP_EXIT=$?
ls ref/Bisulfite_Genome 2>/dev/null && echo PREP_OK
df -h / | tail -1
