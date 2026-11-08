#!/bin/bash
# 제브라피시 WGBS — fastq 를 먼저 작게 자른다
#   bismark 의 --parallel 은 입력 전체를 먼저 조각으로 복사한 뒤 --upto 를 적용한다.
#   그래서 --upto 로는 임시파일이 줄지 않는다(실측 21 GB). 원본 자체를 작게 만들어야 한다.
#   8,000,000 read 쌍 = 32,000,000 줄. 전장 약 1.6x 이지만 array 평균에는 충분하다.
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
