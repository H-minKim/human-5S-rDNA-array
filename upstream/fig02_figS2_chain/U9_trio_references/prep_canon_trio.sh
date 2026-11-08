#!/bin/bash
# Build the trio reference folder and list (2026-09-18)
#   ~/trio_v11/cmp4d_trio = the 616 files of cmp4d_v11_fix1 as real copies + new reference files
#   ~/trio_v11/hap_seq_node_trio.tsv = hap_seq_node_v8 with the 4 changed slots replaced + new rows for the parents
#   ~/trio_v11/new_rows.tsv = only the changed or new slots (the targets for unit mapping)
# File naming rule : <source>__<assembly node>__<original name> (arr__58__hiarr_X = ~/hiasm/arr/hiarr_X)
set -eu
T=~/trio_v11; C=$T/cmp4d_trio
[ -e $C ] && { echo "$C 이미 있음"; exit 1; }
mkdir -p $C
cp -L ~/cmp4d_v11_fix1/*.fa $C/
n0=$(ls $C | wc -l)
add() {  # source new_name
  [ -e "$C/$2" ] && { echo "이름 충돌 $2"; exit 1; }
  cp "$1" "$C/$2"; cmp -s "$1" "$C/$2" || { echo "복사 불일치 $2"; exit 1; }
  echo "  추가 $2  $(md5sum < "$C/$2" | cut -c1-8)  $(grep -v '>' "$C/$2" | tr -d '\n' | wc -c) bp"
}
add $T/asmarr/arr_HG002_HPRCy1_mat.fa          asmy1__58__arr_hg002_mat.fa
add ~/hiasm/arr/hiarr_hg002_3.fa               arr__58__hiarr_hg002_3.fa
add $T/asmarr/arr_HG003_GIABont_2.fa           giabont__58__arr_hg003_2.fa
add $T/hifi_HG01785/arr/hiarr_HG01785_hifi_2.fa hifiarr__58__hiarr_HG01785_hifi_2.fa
NEWP=""
for S in HG01250 HG01251 HG03815 HG00556 HG01190 HG00704 HG01973 HG02614 HG03668; do
  for f in ~/hiasm/arr/hiarr_${S}_*.fa; do add "$f" "arr__58__$(basename $f)"; done
done
echo "  정본 폴더 : 기존 $n0 + 새 $(( $(ls $C | wc -l) - n0 )) = $(ls $C | wc -l) · 링크 $(find $C -type l | wc -l)"
# list
L=$T/hap_seq_node_trio.tsv; R=$T/new_rows.tsv
[ -e $L ] && { echo "$L 이미 있음"; exit 1; }
awk -F'\t' 'BEGIN{OFS="\t"}
  $1=="hg002" && $2=="1" {$3="asmy1__58__arr_hg002_mat.fa"}
  $1=="hg002" && $2=="2" {$3="arr__58__hiarr_hg002_3.fa"}
  $1=="hg003" && $2=="1" {$3="giabont__58__arr_hg003_2.fa"}
  $1=="HG01785" && $2=="1" {$3="hifiarr__58__hiarr_HG01785_hifi_2.fa"}
  {print}' ~/hap_seq_node_v8.tsv > $L
: > $R
printf "hg002\t1\tasmy1__58__arr_hg002_mat.fa\nhg002\t2\tarr__58__hiarr_hg002_3.fa\nhg003\t1\tgiabont__58__arr_hg003_2.fa\nHG01785\t1\thifiarr__58__hiarr_HG01785_hifi_2.fa\n" >> $R
for S in HG01250 HG01251 HG03815 HG00556 HG01190 HG00704 HG01973 HG02614 HG03668; do
  grep -q -P "^$S\t" ~/hap_seq_node_v8.tsv && { echo "$S 이미 목록에 있음"; exit 1; }
  for f in ~/hiasm/arr/hiarr_${S}_*.fa; do k=$(basename $f .fa); k=${k##*_}; printf "%s\t%s\tarr__58__%s\n" $S $k "$(basename $f)"; done
done | tee -a $R >> $L
echo "  목록 : v8 $(wc -l < ~/hap_seq_node_v8.tsv) 행 → trio $(wc -l < $L) 행 · 지도 대상 $(wc -l < $R) 행"
diff <(sort ~/hap_seq_node_v8.tsv) <(sort $L) | grep -c '^[<>]' | sed 's/^/  v8 대비 달라진 줄(교체 4 x 2 + 추가) : /'
cat $R
echo PREP_OK
