#!/bin/bash
# v11-fix1 ARRAY_ONLY edition 2, redone with G6 map correction edition 2 (mkg6map_fix1b: removes only units that were actually lost).
#   Edition 1 outputs (the AO parts of ~/fix1/aog8f1, aoblastf1, bfinalf1, finalf1) are kept, not deleted. Edition 2 goes to new folders *_b.
set -eu
S="${1:?개체}"; NP="${2:-4}"
cd ~
H=$(hostname); F=~/fix1
A=$F/aog8f1b; B=$F/aoblastf1b; C=$F/bfinalf1b; R=$F/finalf1b
for d in $A $B $C $R; do [ -e $d ] && { echo "$d 이미 있음 - 중단"; exit 1; }; done
python3 $F/mkg6map_fix1b.py ~/aog8v10/MAP_G6C2.tsv ~/g8canonv10/CONSENSUS_MAP_G8_full.tsv ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv \
    $F/CANON_V11_FIX1_PATCH.tsv $F/MAP_G6C2_f1b.tsv "HG00133:2:122:125,HG03270:1:43:43"
mkdir -p $A/cmp4d_v7
cp ~/aog8v10/cmp4d_v7/*.fa $A/cmp4d_v7/
cp ~/cmp4d_v11_fix1/asmunits3__58__arr_HG00133_2.fa ~/cmp4d_v11_fix1/asmunits3__57__arr_HG03270_2.fa $A/cmp4d_v7/
cp ~/aog8v10/AOJV3.tsv ~/aog8v10/AOJV5.tsv ~/aog8v10/AOJV8.tsv ~/aog8v10/AOJV_GIAB.tsv $A/
cp $F/MAP_G6C2_f1b.tsv $A/MAP_G6C2.tsv
cp ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv $A/MAP_G8.tsv
cp ~/g8canonv11f1/CONSENSUS_SUM_G8_full.tsv $A/SUM_G8.tsv
ln -s ~/s5s_tuc $A/s5s_tuc
ln -s ~/g8in $A/g8in
ln -s ~/aoreads_v10 $A/aoreads
echo "$S" > $A/samples.txt
python3 ~/ao_g8cut.py $A $A/samples.txt "$NP" > $A/run.log 2>&1
echo "   ao_g8cut : 유닛 $(wc -l < $A/AO_G8_UNITS.tsv) · read $(wc -l < $A/AO_G8_READS.tsv) · $(tail -1 $A/run.log)"
mkdir -p $B/cmp/g8v3 $B/ref $B/aoreads
cp ~/aoblast_v10/ref/* $B/ref/
for f in ~/aoreads_v10/$S.*.fastq.gz; do ln -s "$f" $B/aoreads/; done
cp ~/aog8v10/AOJV3.tsv ~/aog8v10/AOJV5.tsv ~/aog8v10/AOJV8.tsv ~/aog8v10/AOJV_GIAB.tsv $B/cmp/
cp $F/MAP_G6C2_f1b.tsv $B/cmp/CONSENSUS_MAP_G6C2.tsv
cp ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv $B/cmp/g8v3/CONSENSUS_MAP_G8_full.tsv
python3 ~/ao_blastcut.py $B $B/cmp "$NP" > $B/run.log 2>&1
echo "   ao_blastcut : 유닛 $(wc -l < $B/AO_BLAST_UNITS.tsv) · $(tail -1 $B/run.log)"
mkdir -p $C/cmp/g8cut_mid/run3 $C/cmp/g8cut_mid/final $C/lb/full
cp $A/AO_G8_UNITS.tsv $C/cmp/g8cut_mid/run3/AO_G8_UNITS_"$H".raw.tsv
cp $A/AO_G8_READS.tsv $C/cmp/g8cut_mid/run3/AO_G8_READS_"$H".raw.tsv
cp $B/AO_BLAST_UNITS.tsv $C/lb/full/AO_BLAST_UNITS.tsv
PYTHONIOENCODING=utf-8 python3 ~/build_ao_final_g.py $C/cmp $C/lb > $C/run.log 2>&1
cp $C/cmp/g8cut_mid/final/AO_G8_UNITS_"$H".final0.tsv $A/AO_G8_UNITS.final0f1.tsv
echo "   build_ao_final_g : $(wc -l < $A/AO_G8_UNITS.final0f1.tsv) 행"
mkdir -p $R
for f in recall_edge recall_mid recall_pad2; do
  sed "s|~/aog8/sam/|~/fix1/aog8f1b/sam/|" ~/$f.py > $F/${f}_f1b.py
  grep -q "fix1/aog8f1b/sam" $F/${f}_f1b.py
done
python3 $F/recall_edge_f1b.py AO $A/AO_G8_UNITS.final0f1.tsv $R/AO_G8_UNITS_"$H".edge.tsv > $R/ao_edge.log 2>&1
python3 $F/recall_mid_f1b.py  AO $R/AO_G8_UNITS_"$H".edge.tsv $R/AO_G8_UNITS_"$H".mid.tsv > $R/ao_mid.log 2>&1
PAD=100 CTX=200 python3 $F/recall_pad2_f1b.py AO $R/AO_G8_UNITS_"$H".mid.tsv $R/AO_G8_UNITS_"$H".pad.tsv > $R/ao_pad.log 2>&1
echo "   AO 최종 : $(wc -l < $R/AO_G8_UNITS_$H.pad.tsv) 행 · $(cat $R/ao_pad.log)"
echo "FIX1AOBDONE $H $S"
