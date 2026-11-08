#!/bin/bash
# v11-fix2 ARRAY_ONLY re-cut + ARRAY_ONLY methylation. Only the G6 map is replaced by the sequence version (MAP_G6C2_f2u: unit sequence check)
#   other inputs are the same as v11-fix1: canonical cmp4d_v11_fix1, G8 map g8canonv11f1, stored placement AOJV, original v11 code
#   Arguments: tag individuals...   Output: ~/fix2/<tag>/{aog8,aoblast,bfinal,final,placed_bam_g8,meth_ao}  (stops if present)
set -eu
TAG="${1:?태그}"; shift
[ $# -ge 1 ] || { echo "개체 없음"; exit 1; }
cd ~
H=$(hostname); F=~/fix2; W=$F/$TAG
MAP6=$F/MAP_G6C2_f2u.tsv
EXP6="${EXP6:?기대 md5}"
[ -e $W ] && { echo "$W 이미 있음 - 중단"; exit 1; }
[ "$(md5sum < $MAP6 | cut -c1-32)" = "$EXP6" ] || { echo "G6 지도 md5 다름 - 중단"; exit 1; }
for x in ao_g8cut ao_blastcut build_ao_final_g recall_edge recall_mid recall_pad2 mkbam_ao_v11 meth_units_ao_v11; do
  echo "   코드 $x $(md5sum < ~/$x.py | cut -c1-8)"
done
A=$W/aog8; B=$W/aoblast; C=$W/bfinal; R=$W/final; PB=$W/placed_bam_g8; M=$W/meth_ao
mkdir -p $A/cmp4d_v7
cp ~/aog8v10/cmp4d_v7/*.fa $A/cmp4d_v7/
cp ~/cmp4d_v11_fix1/asmunits3__58__arr_HG00133_2.fa ~/cmp4d_v11_fix1/asmunits3__57__arr_HG03270_2.fa $A/cmp4d_v7/
# canonical text check
#   (a) every file with the same name is character-identical to the new fix1 canonical
#   (b) two slots of the run individual: the file the code reads (g8in/hap_seq_node_v7) = the current canonical file (hap_seq_node_v8 -> cmp4d_v11_fix1)
#   extra files (names not in the new canonical folder) are only listed. In v11 these are old hifiasm files of hg002 h2 (ARRAY_ONLY reads 0)
nd=0; nf=0; nx=0
for f in $A/cmp4d_v7/*.fa; do
  b=$(basename "$f")
  if [ -e ~/cmp4d_v11_fix1/$b ]; then nf=$((nf+1)); cmp -s "$f" ~/cmp4d_v11_fix1/$b || nd=$((nd+1)); else nx=$((nx+1)); echo "   여분 파일 $b"; fi
done
echo "   정본 대조 (가) : 같은 이름 $nf 파일 · 다른 파일 $nd · 여분 $nx"
[ "$nd" -eq 0 ] || { echo "정본 다름 - 중단"; exit 1; }
for S in "$@"; do for sl in 1 2; do
  fa=$(awk -F'\t' -v s=$S -v l=$sl '$1==s && $2==l {print $3}' ~/g8in/hap_seq_node_v7.tsv)
  fb=$(awk -F'\t' -v s=$S -v l=$sl '$1==s && $2==l {print $3}' ~/hap_seq_node_v8.tsv)
  [ -z "$fa" ] && [ -z "$fb" ] && continue
  if cmp -s "$A/cmp4d_v7/$fa" ~/cmp4d_v11_fix1/"$fb"; then echo "   정본 대조 (나) $S h$sl : 동일 ($fa)"; else echo "   정본 대조 (나) $S h$sl : 다름 ($fa vs $fb) - 중단"; exit 1; fi
done; done
cp ~/aog8v10/AOJV3.tsv ~/aog8v10/AOJV5.tsv ~/aog8v10/AOJV8.tsv ~/aog8v10/AOJV_GIAB.tsv $A/
cp $MAP6 $A/MAP_G6C2.tsv
cp ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv $A/MAP_G8.tsv
cp ~/g8canonv11f1/CONSENSUS_SUM_G8_full.tsv $A/SUM_G8.tsv
ln -s ~/s5s_tuc $A/s5s_tuc
ln -s ~/g8in $A/g8in
ln -s ~/aoreads_v10 $A/aoreads
printf "%s\n" "$@" > $A/samples.txt
python3 ~/ao_g8cut.py $A $A/samples.txt 4 > $A/run.log 2>&1
echo "   ao_g8cut : 유닛 $(wc -l < $A/AO_G8_UNITS.tsv) · read $(wc -l < $A/AO_G8_READS.tsv) · $(tail -1 $A/run.log)"
mkdir -p $B/cmp/g8v3 $B/ref $B/aoreads
cp ~/aoblast_v10/ref/* $B/ref/
for S in "$@"; do for f in ~/aoreads_v10/$S.*.fastq.gz; do ln -s "$f" $B/aoreads/; done; done
cp ~/aog8v10/AOJV3.tsv ~/aog8v10/AOJV5.tsv ~/aog8v10/AOJV8.tsv ~/aog8v10/AOJV_GIAB.tsv $B/cmp/
cp $MAP6 $B/cmp/CONSENSUS_MAP_G6C2.tsv
cp ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv $B/cmp/g8v3/CONSENSUS_MAP_G8_full.tsv
python3 ~/ao_blastcut.py $B $B/cmp 4 > $B/run.log 2>&1
echo "   ao_blastcut : 유닛 $(wc -l < $B/AO_BLAST_UNITS.tsv) · $(tail -1 $B/run.log)"
mkdir -p $C/cmp/g8cut_mid/run3 $C/cmp/g8cut_mid/final $C/lb/full
cp $A/AO_G8_UNITS.tsv $C/cmp/g8cut_mid/run3/AO_G8_UNITS_"$H".raw.tsv
cp $A/AO_G8_READS.tsv $C/cmp/g8cut_mid/run3/AO_G8_READS_"$H".raw.tsv
cp $B/AO_BLAST_UNITS.tsv $C/lb/full/AO_BLAST_UNITS.tsv
PYTHONIOENCODING=utf-8 python3 ~/build_ao_final_g.py $C/cmp $C/lb > $C/run.log 2>&1
cp $C/cmp/g8cut_mid/final/AO_G8_UNITS_"$H".final0.tsv $A/AO_G8_UNITS.final0f2.tsv
echo "   build_ao_final_g : $(wc -l < $A/AO_G8_UNITS.final0f2.tsv) 행"
mkdir -p $R
for f in recall_edge recall_mid recall_pad2; do
  sed "s|~/aog8/sam/|~/fix2/$TAG/aog8/sam/|" ~/$f.py > $W/${f}_f2.py
  [ "$(diff ~/$f.py $W/${f}_f2.py | grep -c '^>')" -ge 1 ] || { echo "$f sed 실패"; exit 1; }
done
python3 $W/recall_edge_f2.py AO $A/AO_G8_UNITS.final0f2.tsv $R/AO_G8_UNITS_"$H".edge.tsv > $R/ao_edge.log 2>&1
python3 $W/recall_mid_f2.py  AO $R/AO_G8_UNITS_"$H".edge.tsv $R/AO_G8_UNITS_"$H".mid.tsv > $R/ao_mid.log 2>&1
PAD=100 CTX=200 python3 $W/recall_pad2_f2.py AO $R/AO_G8_UNITS_"$H".mid.tsv $R/AO_G8_UNITS_"$H".pad.tsv > $R/ao_pad.log 2>&1
echo "   AO 최종 : $(wc -l < $R/AO_G8_UNITS_$H.pad.tsv) 행"
# ARRAY_ONLY methylation
sed -e "s|D = H(\"~/aog8v10\")|D = H(\"~/fix2/$TAG/aog8\")|" -e "s|OUT = H(\"~/placed_bam_g8_v11\")|OUT = H(\"~/fix2/$TAG/placed_bam_g8\")|" ~/mkbam_ao_v11.py > $W/mkbam_ao_f2.py
sed -e 's|MAPF = os.path.join(H, "g8canonv10", "CONSENSUS_MAP_G8_full.tsv")|MAPF = os.path.join(H, "g8canonv11f1", "CONSENSUS_MAP_G8_full.tsv")|' \
    -e "s|os.path.join(H, \"placed_bam_g8_v11\")|os.path.join(H, \"fix2\", \"$TAG\", \"placed_bam_g8\")|" ~/meth_units_ao_v11.py > $W/meth_units_ao_f2.py
n1=$(diff ~/mkbam_ao_v11.py $W/mkbam_ao_f2.py | grep -c '^>'); n2=$(diff ~/meth_units_ao_v11.py $W/meth_units_ao_f2.py | grep -c '^>')
echo "   sed 바뀐 줄 : mkbam $n1 · meth_ao $n2"
[ "$n1" -eq 2 ] && [ "$n2" -eq 2 ] || { echo "sed 줄 수 다름 - 중단"; exit 1; }
printf "%s\n" "$@" > $W/meth_list.txt
python3 $W/mkbam_ao_f2.py > $W/mkbam.log 2>&1
echo "   AO BAM : $(tr '\n' ' ' < $W/mkbam.log)"
python3 $W/meth_units_ao_f2.py $M $W/meth_list.txt > $W/meth_ao.log 2>&1
echo "   메틸 AO : $(tail -2 $W/meth_ao.log | tr '\n' ' ')"
echo "FIX2DONE $H $TAG $*"
