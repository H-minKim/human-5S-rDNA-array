#!/bin/bash
# GSE122695 / PRJNA506038 — vasa:eGFP 생식세포(GSC) · 체세포(Control) 저심도 PBAT WGBS 정렬
#   인자 : 목록파일(각 줄 "run<TAB>시료명<TAB>ftp")
#   PBAT 라 bismark --pbat · 단일말단
#   디스크 25 GB 가드
set -u
cd ~/zfwgbs || exit 1
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p gsc/fq gsc/bam gsc/log tmp
LIST=$1
guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP — 디스크 여유 ${F}G"; exit 1; fi
}
while IFS=$'\t' read -r R NAME URL; do
  [ -z "$R" ] && continue
  [ -f gsc/bam/${R}.done ] && { echo "SKIP $R $NAME"; continue; }
  guard
  [ -f gsc/fq/${R}.fastq.gz ] || curl -s -o gsc/fq/${R}.fastq.gz "https://$URL"
  python3 gsc_trim.py gsc/fq/${R}.fastq.gz gsc/fq/${R}.trim.fq.gz 2>> gsc/log/trim.log
  bismark --pbat --genome ref -p 4 --bam --temp_dir tmp -o gsc/bam gsc/fq/${R}.trim.fq.gz > gsc/log/${R}.log 2>&1
  B=$(ls gsc/bam/${R}.trim*bismark*.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 4 -o gsc/bam/${R}.sorted.bam "$B" && samtools index gsc/bam/${R}.sorted.bam && rm -f "$B"
    M=$(grep -m1 "Mapping efficiency" gsc/bam/${R}.trim*_report.txt 2>/dev/null | tr -d '\t' | sed 's/.*://')
    echo -e "$R\t$NAME\t$M" | tee -a gsc/mapping.tsv
    touch gsc/bam/${R}.done
  else
    echo "ALIGN_FAIL $R $NAME"; tail -3 gsc/log/${R}.log
  fi
  rm -f gsc/fq/${R}.trim.fq.gz tmp/* 2>/dev/null
done < "$LIST"
echo GSC_DONE
