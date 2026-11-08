#!/bin/bash
# v11-fix1 ARRAY_ONLY resume run, the first run missed the s5s_tuc link in the work folder and ao_g8cut failed.
#   Add the s5s_tuc link as in setup_v10cut.sh and rerun only ao_g8cut -> build_ao_final_g -> recall (blastcut uses the first-run output as is)
set -eu
S="${1:?개체}"; NP="${2:-4}"
cd ~
H=$(hostname); F=~/fix1; A=$F/aog8f1; B=$F/aoblastf1; C=$F/bfinalf1; R=$F/finalf1
[ -s $B/AO_BLAST_UNITS.tsv ] || { echo "blastcut 산출 없음"; exit 1; }
[ -e $A/s5s_tuc ] || ln -s ~/s5s_tuc $A/s5s_tuc
[ -e $A/AO_G8_UNITS.tsv ] && { echo "AO_G8_UNITS 이미 있음 - 중단"; exit 1; }
rm -rf $A/sam $A/tmp
python3 ~/ao_g8cut.py $A $A/samples.txt "$NP" > $A/run.log 2>&1
echo "   ao_g8cut : 유닛 $(wc -l < $A/AO_G8_UNITS.tsv) · read $(wc -l < $A/AO_G8_READS.tsv) · $(tail -1 $A/run.log)"
rm -rf $C; mkdir -p $C/cmp/g8cut_mid/run3 $C/cmp/g8cut_mid/final $C/lb/full
cp $A/AO_G8_UNITS.tsv $C/cmp/g8cut_mid/run3/AO_G8_UNITS_"$H".raw.tsv
cp $A/AO_G8_READS.tsv $C/cmp/g8cut_mid/run3/AO_G8_READS_"$H".raw.tsv
cp $B/AO_BLAST_UNITS.tsv $C/lb/full/AO_BLAST_UNITS.tsv
PYTHONIOENCODING=utf-8 python3 ~/build_ao_final_g.py $C/cmp $C/lb > $C/run.log 2>&1
cp $C/cmp/g8cut_mid/final/AO_G8_UNITS_"$H".final0.tsv $A/AO_G8_UNITS.final0f1.tsv
echo "   build_ao_final_g : $(wc -l < $A/AO_G8_UNITS.final0f1.tsv) 행 · $(tail -2 $C/run.log | tr '\n' ' ')"
python3 $F/recall_edge_f1.py AO $A/AO_G8_UNITS.final0f1.tsv $R/AO_G8_UNITS_"$H".edge.tsv > $R/ao_edge.log 2>&1
python3 $F/recall_mid_f1.py  AO $R/AO_G8_UNITS_"$H".edge.tsv $R/AO_G8_UNITS_"$H".mid.tsv > $R/ao_mid.log 2>&1
PAD=100 CTX=200 python3 $F/recall_pad2_f1.py AO $R/AO_G8_UNITS_"$H".mid.tsv $R/AO_G8_UNITS_"$H".pad.tsv > $R/ao_pad.log 2>&1
echo "   AO 최종 : $(wc -l < $R/AO_G8_UNITS_$H.pad.tsv) 행 · $(cat $R/ao_pad.log)"
echo "FIX1AODONE $H $S"
