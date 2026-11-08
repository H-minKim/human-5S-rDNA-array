#!/bin/bash
# Potok 2013 계열 정렬이 끝나면 곧바로 구역별 메틸을 잰다
#   같은 표(ZFMETH.tsv)에 이어 붙인다 — 라벨 앞에 P_ 를 붙여 Jiang 계열과 구분한다.
set -u
cd ~/zfwgbs || exit 1
until grep -q POTOK_DONE potok.out 2>/dev/null; do sleep 60; done
echo "정렬 끝 — 측정 시작 $(date -u +%H:%M)"

m () {
  R=$1; L=$2
  [ -f bam/${R}.sorted.bam ] || { echo "BAM 없음 $R ($L)"; return; }
  python3 bismark_region_meth.py bam/${R}.sorted.bam "$L" zf_regions.tsv ZFMETH.tsv
}

m SRR800059 P_알_Potok
m SRR800081 P_근육성체_Potok
m SRR800056 P_정자_Potok
m SRR800077 P_sphere_Potok
m SRR800061 P_2-16세포_Potok
m SRR800066 P_64세포_Potok
m SRR800070 P_256세포_Potok
wc -l ZFMETH.tsv
echo POTOK_METH_DONE
