#!/bin/bash
# As soon as the Potok 2013 alignments finish, measure the methylation per region.
#   The results are appended to the same table, ZFMETH.tsv, with the labels prefixed P_ to distinguish them from the Jiang series.
set -u
cd ~/zfwgbs || exit 1
until grep -q POTOK_DONE potok.out 2>/dev/null; do sleep 60; done
echo "alignment finished: measurement start $(date -u +%H:%M)"

m () {
  R=$1; L=$2
  [ -f bam/${R}.sorted.bam ] || { echo "BAM none $R ($L)"; return; }
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
