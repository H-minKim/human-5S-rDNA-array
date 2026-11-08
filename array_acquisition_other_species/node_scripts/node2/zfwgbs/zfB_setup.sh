#!/bin/bash
# Zebrafish WGBS preparation on node2: the assembly and the reads are fetched from source
#   노드끼리 ssh 가 안 되므로 각자 받는 편이 빠르다.
#   read 는 fastq-dump -X 로 앞 N spot 만 받아 큰 중간파일을 만들지 않는다.
set -u
cd ~/zfwgbs
mkdir -p ref sra logs bam tmp
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH

# 1) 조립본 (GCA_049306965.1 = GRCz12tu)
if [ ! -f ref/zebrafish_GRCz12tu.fa ]; then
  curl -sSL -o ncbi.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_049306965.1/download?include_annotation_type=GENOME_FASTA" \
    && unzip -o -q ncbi.zip -d ncbi \
    && find ncbi -name '*.fna' -exec cp {} ref/zebrafish_GRCz12tu.fa \; \
    && rm -rf ncbi ncbi.zip
fi
ls -la ref/ | tail -2
grep -c '^>' ref/zebrafish_GRCz12tu.fa

# 2) read — 앞 800 만 spot 만
if ! command -v fastq-dump >/dev/null 2>&1; then
  sudo apt-get install -y -qq sra-toolkit >/dev/null 2>&1
fi
for S in SRR14703401 SRR14703403; do
  [ -f sra/${S}_1.fastq.gz ] && { echo "SKIP $S"; continue; }
  fastq-dump -X 8000000 --split-files --gzip -O sra "$S" > logs/${S}.dl.log 2>&1
  ls -la sra/${S}_*.fastq.gz 2>/dev/null | awk '{printf "%s %.2f GB\n", $NF, $5/1e9}'
done
df -h / | tail -1
echo SETUP_DONE
