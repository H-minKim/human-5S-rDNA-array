#!/bin/bash
# 깊이 올린 BAM 으로 read 단위 메틸을 다시 잰다
set -u
cd ~/zfwgbs || exit 1
until grep -q DEEPEN_DONE deepen.out 2>/dev/null; do sleep 90; done
echo "깊이 올리기 끝 — 재측정 $(date -u +%H:%M)"
rm -f ZFREAD_DEEP.tsv
for X in "SRR800059_deep P_알_Potok_deep" "SRR800081_deep P_근육성체_Potok_deep" \
         "egg_deep J_알_Jiang_deep" "SRR12797061 O_장" "SRR837895 J_Germring" "SRR800077 P_sphere_Potok"; do
  set -- $X
  [ -f bam/$1.bam ] && B=bam/$1.bam || B=bam/$1.sorted.bam
  [ -f "$B" ] && python3 bismark_read_meth.py "$B" "$2" zf_regions.tsv ZFREAD_DEEP.tsv 3
done
wc -l ZFREAD_DEEP.tsv
echo DEEP_METH_DONE
