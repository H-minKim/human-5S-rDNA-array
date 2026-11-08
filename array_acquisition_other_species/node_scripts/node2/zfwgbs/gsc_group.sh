#!/bin/bash
# 시료를 군별로 합쳐 5S 구역 메틸을 잰다(시료 하나로는 read 가 너무 적다)
#   군 : F_GSC(35–70d) · M_GSC(35–70d) · early_GSC(1–28d, 성 미분화) · F_Control · M_Control · early_Control
#   + 시기별(35·40·45·50·70d) F_GSC / F_Control
set -u
cd ~/zfwgbs || exit 1
until grep -q GSC_DONE gsc/run_all.out 2>/dev/null; do sleep 60; done
mkdir -p gsc/grp
rm -f gsc/GSC_PARTS.tsv
# run → 시료명
awk -F'\t' '{print $1"\t"$2}' all.list > gsc/names.tsv
mk () {  # 군이름 · 이름정규식
  G=$1; P=$2
  L=$(awk -F'\t' -v p="$P" '$2 ~ p {print "gsc/bam/"$1".sorted.bam"}' gsc/names.tsv | while read f; do [ -f "$f" ] && echo "$f"; done)
  N=$(echo "$L" | grep -c bam)
  [ "$N" -eq 0 ] && { echo "빈 군 $G"; return; }
  samtools merge -f -@ 4 gsc/grp/$G.bam $L && samtools index gsc/grp/$G.bam
  echo "$G 시료 $N · read $(samtools view -c gsc/grp/$G.bam)"
  python3 bismark_parts_meth.py gsc/grp/$G.bam "$G" zf_regions_parts.tsv gsc/GSC_PARTS.tsv 2000 >/dev/null
}
mk F_GSC_35to70      '_(35|40|45|50|70)d_F_GSC_'
mk M_GSC_35to70      '_(35|40|45|50|70)d_M_GSC_'
mk early_GSC_1to28   '_(01|02|04|07|11|14|21|25|28)d_GSC_'
mk F_Control_35to70  '_(35|40|45|50|70)d_F_Control_'
mk M_Control_35to70  '_(35|40|45|50|70)d_M_Control_'
mk early_Control_1to28 '_(01|02|04|07|11|14|21|25|28)d_Control_'
for D in 35 40 45 50 70; do
  mk F_GSC_${D}d "_${D}d_F_GSC_"
  mk M_GSC_${D}d "_${D}d_M_GSC_"
done
wc -l gsc/GSC_PARTS.tsv
echo GROUP_DONE
