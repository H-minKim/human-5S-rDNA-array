#!/bin/bash
set -u
cd ~/zfwgbs || exit 1
rm -f ZFPARTS.tsv
for X in "egg_deep 알_Jiang" "SRR800059_deep 알_Potok" "SRR837878 16세포_1.5hpf" \
         "SRR837894 1k세포_3hpf_MBT" "SRR800077 sphere_4hpf" "SRR837895 Germring_6hpf" \
         "SRR11614917 배아24hpf" "SRR800081_deep 근육성체_Potok" "SRR12797061 장성체" \
         "SRR12797059 심장성체" "SRR837872 정자_Jiang" "SRR800056 정자_Potok" "SRR837899 정소"; do
  set -- $X
  [ -f bam/$1.bam ] && B=bam/$1.bam || B=bam/$1.sorted.bam
  [ -f "$B" ] && python3 bismark_parts_meth.py "$B" "$2" zf_regions_parts.tsv ZFPARTS.tsv 2000
done
wc -l ZFPARTS.tsv
echo PARTS_DONE
