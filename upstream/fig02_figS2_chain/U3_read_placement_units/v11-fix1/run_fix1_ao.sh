#!/bin/bash
# v11-fix1 ARRAY_ONLY: ao_g8cut -> build_ao_final_g -> recall.
#   Adds the s5s_tuc link as in setup_v10cut.sh. The ao_blastcut output must already exist and is used as is.
set -eu
S="${1:?individual}"; NP="${2:-4}"
cd ~
H=$(hostname); F=~/fix1; A=$F/aog8f1; B=$F/aoblastf1; C=$F/bfinalf1; R=$F/finalf1
[ -s $B/AO_BLAST_UNITS.tsv ] || { echo "no blastcut output"; exit 1; }
[ -e $A/s5s_tuc ] || ln -s ~/s5s_tuc $A/s5s_tuc
[ -e $A/AO_G8_UNITS.tsv ] && { echo "AO_G8_UNITS already exists - aborting"; exit 1; }
rm -rf $A/sam $A/tmp
python3 ~/ao_g8cut.py $A $A/samples.txt "$NP" > $A/run.log 2>&1
echo "   ao_g8cut : units $(wc -l < $A/AO_G8_UNITS.tsv), read $(wc -l < $A/AO_G8_READS.tsv), $(tail -1 $A/run.log)"
rm -rf $C; mkdir -p $C/cmp/g8cut_mid/run3 $C/cmp/g8cut_mid/final $C/lb/full
cp $A/AO_G8_UNITS.tsv $C/cmp/g8cut_mid/run3/AO_G8_UNITS_"$H".raw.tsv
cp $A/AO_G8_READS.tsv $C/cmp/g8cut_mid/run3/AO_G8_READS_"$H".raw.tsv
cp $B/AO_BLAST_UNITS.tsv $C/lb/full/AO_BLAST_UNITS.tsv
PYTHONIOENCODING=utf-8 python3 ~/build_ao_final_g.py $C/cmp $C/lb > $C/run.log 2>&1
cp $C/cmp/g8cut_mid/final/AO_G8_UNITS_"$H".final0.tsv $A/AO_G8_UNITS.final0f1.tsv
echo "   build_ao_final_g : $(wc -l < $A/AO_G8_UNITS.final0f1.tsv) rows, $(tail -2 $C/run.log | tr '\n' ' ')"
python3 $F/recall_edge_f1.py AO $A/AO_G8_UNITS.final0f1.tsv $R/AO_G8_UNITS_"$H".edge.tsv > $R/ao_edge.log 2>&1
python3 $F/recall_mid_f1.py  AO $R/AO_G8_UNITS_"$H".edge.tsv $R/AO_G8_UNITS_"$H".mid.tsv > $R/ao_mid.log 2>&1
PAD=100 CTX=200 python3 $F/recall_pad2_f1.py AO $R/AO_G8_UNITS_"$H".mid.tsv $R/AO_G8_UNITS_"$H".pad.tsv > $R/ao_pad.log 2>&1
echo "   AO final : $(wc -l < $R/AO_G8_UNITS_$H.pad.tsv) rows, $(cat $R/ao_pad.log)"
echo "FIX1AODONE $H $S"
