#!/bin/bash
# 제브라피시 WGBS — 참조 준비 + read 내려받기
set -e
cd ~/zfwgbs
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p ref sra logs
if [ ! -f ref/zebrafish_GRCz12tu.fna ]; then
  cp ~/xsp/asm/zebrafish_GRCz12tu.fna ref/ 2>/dev/null || \
    scp -o StrictHostKeyChecking=no -i ~/.ssh/id_rsa ubuntu@$NODE_HOST_1:~/xsp/asm/zebrafish_GRCz12tu.fna ref/
fi
if [ ! -d ref/Bisulfite_Genome ]; then
  bismark_genome_preparation --parallel 8 --bowtie2 ref > logs/prep.log 2>&1
fi
echo PREP_DONE
