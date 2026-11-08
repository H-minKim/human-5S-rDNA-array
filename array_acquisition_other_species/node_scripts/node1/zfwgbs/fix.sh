#!/bin/bash
cd ~/zfwgbs
# 1) fastq 압축으로 디스크 확보 (bismark 는 .gz 를 읽는다)
for f in sra/*.fastq; do [ -f "$f" ] && pigz -p 8 "$f" 2>/dev/null || gzip "$f"; done
df -h / | tail -1
# 2) bismark 는 .fa/.fasta 만 받는다 → 이름 바꾸기
[ -f ref/zebrafish_GRCz12tu.fna ] && mv ref/zebrafish_GRCz12tu.fna ref/zebrafish_GRCz12tu.fa
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
bismark_genome_preparation --parallel 6 --bowtie2 ref > logs/prep2.log 2>&1
echo PREP_EXIT=$?
ls ref/Bisulfite_Genome 2>/dev/null && echo PREP_OK
df -h / | tail -1
