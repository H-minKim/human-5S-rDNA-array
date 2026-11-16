#!/bin/bash
# Zebrafish WGBS: cut the fastq down first.
#   bismark's --parallel copies the whole input into chunks before applying --upto,
#   so --upto does not reduce the temporary files; measured, they still reach 21 GB. The input itself must be made smaller.
#   8,000,000 read pairs is 32,000,000 lines. That is about 1.6x over the whole genome, which is enough for an array mean.
cd ~/zfwgbs/sra
N=32000000
for S in SRR14703401 SRR14703403; do
  if [ -f sub_${S}_1.fastq.gz ]; then echo "SKIP $S"; continue; fi
  for M in 1 2; do
    zcat ${S}_${M}.fastq.gz | head -n $N | pigz -p 6 > sub_${S}_${M}.fastq.gz
  done
  ls -la sub_${S}_*.fastq.gz | awk '{printf "%s %.2f GB\n", $NF, $5/1e9}'
done
df -h / | tail -1
echo SUB_DONE
