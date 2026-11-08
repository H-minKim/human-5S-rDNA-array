#!/bin/bash
# v11-fix1 methylation for one individual: ARRAY_ONLY BAM (mkbam_ao_v11) -> unit methylation AO (meth_units_ao_v11) and junction (meth_units_jx_v11)
#   Changes: work folder (aog8f1b), BAM output folder, maps (AO=g8canonv11f1, JX=g8jxmap_f1). The junction BAM is the original ~/placed_bam unchanged.
set -eu
S="${1:?개체}"
cd ~
H=$(hostname); F=~/fix1
for x in ~/mkbam_ao_v11.py ~/meth_units_ao_v11.py ~/meth_units_jx_v11.py; do [ -s $x ] || { echo "없음 $x"; exit 1; }; done
for d in $F/placed_bam_g8_f1 $F/methf1_ao $F/methf1_jx; do [ -e $d ] && { echo "$d 이미 있음"; exit 1; }; done
echo "== $H $S md5 : mkbam $(md5sum < ~/mkbam_ao_v11.py | cut -c1-8) · meth_ao $(md5sum < ~/meth_units_ao_v11.py | cut -c1-8) · meth_jx $(md5sum < ~/meth_units_jx_v11.py | cut -c1-8)"
sed -e 's|D = H("~/aog8v10")|D = H("~/fix1/aog8f1b")|' -e 's|OUT = H("~/placed_bam_g8_v11")|OUT = H("~/fix1/placed_bam_g8_f1")|' ~/mkbam_ao_v11.py > $F/mkbam_ao_f1.py
sed -e 's|MAPF = os.path.join(H, "g8canonv10", "CONSENSUS_MAP_G8_full.tsv")|MAPF = os.path.join(H, "g8canonv11f1", "CONSENSUS_MAP_G8_full.tsv")|' \
    -e 's|os.path.join(H, "placed_bam_g8_v11")|os.path.join(H, "fix1", "placed_bam_g8_f1")|' ~/meth_units_ao_v11.py > $F/meth_units_ao_f1.py
sed -e 's|MAPF = os.path.join(H, "g8jxmap_v10", "CONSENSUS_MAP_G8_full.tsv")|MAPF = os.path.join(H, "fix1", "g8jxmap_f1", "CONSENSUS_MAP_G8_full.tsv")|' ~/meth_units_jx_v11.py > $F/meth_units_jx_f1.py
for p in "mkbam_ao_v11 mkbam_ao_f1" "meth_units_ao_v11 meth_units_ao_f1" "meth_units_jx_v11 meth_units_jx_f1"; do set -- $p; echo "   $1 -> $2 바뀐 줄 $(diff ~/$1.py $F/$2.py | grep -c '^[<>]')"; done
echo "$S" > $F/meth_list.txt
python3 $F/mkbam_ao_f1.py > $F/mkbam_f1.log 2>&1
echo "   AO BAM : $(cat $F/mkbam_f1.log | tr '\n' ' ')"
python3 $F/meth_units_ao_f1.py $F/methf1_ao $F/meth_list.txt > $F/methf1_ao.log 2>&1
python3 $F/meth_units_jx_f1.py $F/methf1_jx $F/meth_list.txt > $F/methf1_jx.log 2>&1
echo "   메틸 AO : $(tail -2 $F/methf1_ao.log | tr '\n' ' ') · 행 $(zcat $F/methf1_ao/METH_READUNIT_$H.tsv.gz | wc -l)"
echo "   메틸 JX : $(tail -2 $F/methf1_jx.log | tr '\n' ' ') · 행 $(zcat $F/methf1_jx/METH_READUNIT_$H.tsv.gz | wc -l)"
echo "FIX1METHDONE $H $S"
