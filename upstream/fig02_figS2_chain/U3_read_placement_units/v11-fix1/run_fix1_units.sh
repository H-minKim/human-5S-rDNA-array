#!/bin/bash
# v11-fix1 : regenerate one individual's read unit tables (ARRAY_ONLY + junction), same v11 chain, only the inputs are the new canonical and map
#   ARRAY_ONLY : ao_g8cut → ao_blastcut → build_ao_final_g → recall_edge/mid/pad2 (only the AO SAM path points to the new folder)
#   junction   : mkjxmap_v10 (repair record v10 7 rows + fix1 2 rows) → jx_g8cut (only map path and input SAM folder) → recall_edge/mid/pad2 (original)
#   changed inputs : two canonical files, map/summary (g8canonv11f1), G6 map (coordinate correction for two slots), junction map (g8jxmap_f1)
#   Existing folders are not touched. Everything goes to new folders under ~/fix1/.
#   arguments : individual, number of parallel jobs
set -u
S="${1:?개체}"; NP="${2:-2}"
cd ~ || exit 1
H=$(hostname)
F=~/fix1
for x in ~/aog8v10 ~/aoblast_v10 ~/cmp4d_v11_fix1 ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv ~/g8canonv11f1/CONSENSUS_SUM_G8_full.tsv \
         $F/CANON_V11_FIX1_PATCH.tsv ~/CANON_V10_PATCH.tsv ~/ao_g8cut.py ~/ao_blastcut.py ~/build_ao_final_g.py \
         ~/jx_g8cut_v10b.py ~/recall_edge.py ~/recall_mid.py ~/recall_pad2.py ~/rjrepro/$S.placed.sam $F/mkg6map_fix1.py ~/mkjxmap_v10.py; do
  [ -e "$x" ] || { echo "$H 없음 $x"; exit 1; }
done
for d in $F/aog8f1 $F/aoblastf1 $F/bfinalf1 $F/finalf1 $F/jxf1 $F/rjsub; do [ -e "$d" ] && { echo "$H $d 이미 있음 - 중단"; exit 1; }; done
echo "== $H $S 코드 md5 : ao_g8cut $(md5sum < ~/ao_g8cut.py | cut -c1-8) · ao_blastcut $(md5sum < ~/ao_blastcut.py | cut -c1-8) · build_ao_final_g $(md5sum < ~/build_ao_final_g.py | cut -c1-8) · jx_g8cut_v10b $(md5sum < ~/jx_g8cut_v10b.py | cut -c1-8) · recall $(md5sum < ~/recall_edge.py | cut -c1-8)/$(md5sum < ~/recall_mid.py | cut -c1-8)/$(md5sum < ~/recall_pad2.py | cut -c1-8)"

# ── common : G6 map correction and junction map
python3 $F/mkg6map_fix1.py ~/aog8v10/MAP_G6C2.tsv $F/CANON_V11_FIX1_PATCH.tsv $F/MAP_G6C2_f1.tsv
mkdir -p ~/g8jxmap_f1_tmp
( cat ~/CANON_V10_PATCH.tsv; tail -n +2 $F/CANON_V11_FIX1_PATCH.tsv ) > $F/CANON_V11_FIX1_PATCH_ALL.tsv
mkdir -p $F/g8jxmap_f1
python3 ~/mkjxmap_v10.py ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv $F/CANON_V11_FIX1_PATCH_ALL.tsv $F/g8jxmap_f1/CONSENSUS_MAP_G8_full.tsv
rmdir ~/g8jxmap_f1_tmp

# ── ARRAY_ONLY ① ao_g8cut : work folder has the same layout as aog8v10 (canonical folder is actually copied)
A=$F/aog8f1
mkdir -p $A/cmp4d_v7
cp ~/aog8v10/cmp4d_v7/*.fa $A/cmp4d_v7/
cp ~/cmp4d_v11_fix1/asmunits3__58__arr_HG00133_2.fa ~/cmp4d_v11_fix1/asmunits3__57__arr_HG03270_2.fa $A/cmp4d_v7/
cp ~/aog8v10/AOJV3.tsv ~/aog8v10/AOJV5.tsv ~/aog8v10/AOJV8.tsv ~/aog8v10/AOJV_GIAB.tsv $A/
cp $F/MAP_G6C2_f1.tsv $A/MAP_G6C2.tsv
cp ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv $A/MAP_G8.tsv
cp ~/g8canonv11f1/CONSENSUS_SUM_G8_full.tsv $A/SUM_G8.tsv
ln -s ~/g8in $A/g8in
ln -s ~/aoreads_v10 $A/aoreads
echo "$S" > $A/samples.txt
echo "   aog8f1 정본 $(ls $A/cmp4d_v7 | wc -l) (링크 $(find $A/cmp4d_v7 -type l | wc -l))"
python3 ~/ao_g8cut.py $A $A/samples.txt "$NP" > $A/run.log 2>&1
echo "   ao_g8cut : 유닛 $(wc -l < $A/AO_G8_UNITS.tsv) · read $(wc -l < $A/AO_G8_READS.tsv) · $(tail -1 $A/run.log)"

# ── ARRAY_ONLY ② ao_blastcut
B=$F/aoblastf1
mkdir -p $B/cmp/g8v3 $B/ref $B/aoreads
cp ~/aoblast_v10/ref/* $B/ref/
for f in ~/aoreads_v10/$S.*.fastq.gz; do ln -s "$f" $B/aoreads/; done
cp ~/aog8v10/AOJV3.tsv ~/aog8v10/AOJV5.tsv ~/aog8v10/AOJV8.tsv ~/aog8v10/AOJV_GIAB.tsv $B/cmp/
cp $F/MAP_G6C2_f1.tsv $B/cmp/CONSENSUS_MAP_G6C2.tsv
cp ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv $B/cmp/g8v3/CONSENSUS_MAP_G8_full.tsv
python3 ~/ao_blastcut.py $B $B/cmp "$NP" > $B/run.log 2>&1
echo "   ao_blastcut : 유닛 $(wc -l < $B/AO_BLAST_UNITS.tsv) · $(tail -1 $B/run.log)"

# ── ARRAY_ONLY ③ build_ao_final_g (same folder layout as run_buildfinal.sh)
C=$F/bfinalf1
mkdir -p $C/cmp/g8cut_mid/run3 $C/cmp/g8cut_mid/final $C/lb/full
cp $A/AO_G8_UNITS.tsv $C/cmp/g8cut_mid/run3/AO_G8_UNITS_"$H".raw.tsv
cp $A/AO_G8_READS.tsv $C/cmp/g8cut_mid/run3/AO_G8_READS_"$H".raw.tsv
cp $B/AO_BLAST_UNITS.tsv $C/lb/full/AO_BLAST_UNITS.tsv
PYTHONIOENCODING=utf-8 python3 ~/build_ao_final_g.py $C/cmp $C/lb > $C/run.log 2>&1
cp $C/cmp/g8cut_mid/final/AO_G8_UNITS_"$H".final0.tsv $A/AO_G8_UNITS.final0f1.tsv
echo "   build_ao_final_g : $(wc -l < $A/AO_G8_UNITS.final0f1.tsv) 행"

# ── ARRAY_ONLY ④ recall (only the AO SAM path points to the new folder)
R=$F/finalf1
mkdir -p $R
for f in recall_edge recall_mid recall_pad2; do
  sed "s|~/aog8/sam/|~/fix1/aog8f1/sam/|" ~/$f.py > $F/${f}_f1.py
  grep -q "fix1/aog8f1/sam" $F/${f}_f1.py || { echo "치환 실패 $f"; exit 1; }
done
python3 $F/recall_edge_f1.py AO $A/AO_G8_UNITS.final0f1.tsv $R/AO_G8_UNITS_"$H".edge.tsv > $R/ao_edge.log 2>&1
python3 $F/recall_mid_f1.py  AO $R/AO_G8_UNITS_"$H".edge.tsv $R/AO_G8_UNITS_"$H".mid.tsv > $R/ao_mid.log 2>&1
PAD=100 CTX=200 python3 $F/recall_pad2_f1.py AO $R/AO_G8_UNITS_"$H".mid.tsv $R/AO_G8_UNITS_"$H".pad.tsv > $R/ao_pad.log 2>&1
echo "   AO 최종 : $(wc -l < $R/AO_G8_UNITS_$H.pad.tsv) 행 · $(cat $R/ao_pad.log)"

# ── junction ① jx_g8cut (only map path and input SAM folder)
mkdir -p $F/rjsub
cp ~/rjrepro/$S.placed.sam $F/rjsub/
sed -e 's|H("~/g8jxmap_v10/CONSENSUS_MAP_G8_full.tsv")|H("~/fix1/g8jxmap_f1/CONSENSUS_MAP_G8_full.tsv")|' \
    -e 's|src = H("~/rjrepro") if PASS == "main" else H("~/rjrepro_c2")|src = H("~/fix1/rjsub") if PASS == "main" else H("~/rjrepro_c2")|' \
    ~/jx_g8cut_v10b.py > $F/jx_g8cut_f1.py
echo "   jx_g8cut 바뀐 줄 $(diff ~/jx_g8cut_v10b.py $F/jx_g8cut_f1.py | grep -c '^[<>]')"
mkdir -p $F/jxf1
python3 $F/jx_g8cut_f1.py main $F/jxf1 "$NP" > $F/jxf1/main.log 2>&1
echo "   jx_g8cut : $(wc -l < $F/jxf1/JX_G8_UNITS_main_$H.tsv) 행 · $(tail -1 $F/jxf1/main.log)"
# ── junction ② recall (original code, SAM is the original ~/rjrepro)
J=$F/jxf1/JX_G8_UNITS_main_"$H".tsv
python3 ~/recall_edge.py JX "$J" $R/JX_G8_UNITS_main_"$H".edge.tsv > $R/jx_main_edge.log 2>&1
python3 ~/recall_mid.py  JX $R/JX_G8_UNITS_main_"$H".edge.tsv $R/JX_G8_UNITS_main_"$H".mid.tsv > $R/jx_main_mid.log 2>&1
PAD=100 CTX=200 python3 ~/recall_pad2.py JX $R/JX_G8_UNITS_main_"$H".mid.tsv $R/JX_G8_UNITS_main_"$H".pad.tsv > $R/jx_main_pad.log 2>&1
echo "   JX 최종 : $(wc -l < $R/JX_G8_UNITS_main_$H.pad.tsv) 행 · $(cat $R/jx_main_pad.log)"
echo "FIX1UNITSDONE $H $S"
