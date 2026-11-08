#!/bin/bash
# trio junction placement: v11 rjf_repro.py (unchanged, md5 checked against the local copy) in an isolated HOME
#   HAP_READ_FINAL = trio group table, G6 map/summary = fghome(v11 + trio + trio2), G8 map = fix-1 map with only the changed slots replaced by trio
#   Arguments: individuals...
set -eu
T=~/trio_v11; H=$T/rjhome
if [ ! -e $H ]; then
  mkdir -p $H/s5s_tuc/ref
  cp $T/HAP_READ_FINAL_trio.tsv $H/HAP_READ_FINAL.tsv
  cp ~/s5s_tuc/ref/unit_canon.fa ~/s5s_tuc/ref/FLANK_proximal_50000.fa ~/s5s_tuc/ref/FLANK_distal_50000.fa $H/s5s_tuc/ref/
  for d in kgp_ont giab5s_seq s5s_ont; do [ -e ~/$d ] && ln -s ~/$d $H/$d; done
  DROP='($1=="hg002") || ($1=="hg003") || ($1=="HG01785" && $2=="1")'
  { head -1 ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv
    awk -F'\t' "NR>1 && !($DROP)" ~/g8canonv11f1/CONSENSUS_MAP_G8_full.tsv
    tail -n +2 $T/g8trio/CONSENSUS_MAP_G8_trio.tsv
    tail -n +2 $T/g8trio2/CONSENSUS_MAP_G8_trio2.tsv; } > $H/G8MAP_trio.tsv
  echo "G8 합본 $(wc -l < $H/G8MAP_trio.tsv) 행"
fi
cd $H
HOME=$H HSN=$T/hap_seq_node_trio2.tsv RCFF=$T/fghome/CONSENSUS_SUM_G6.tsv RCFCOL=5 MAPF=$T/fghome/CONSENSUS_MAP_G6.tsv \
  G8MAP=$H/G8MAP_trio.tsv CMPD=$T/cmp4d_trio C2SKIP= RJW=$H/rj TH=4 TAG=trio \
  python3 $T/rjf_repro.py <(printf "%s\n" "$@") > $H/rjf_trio.log 2>&1
echo RJFDONE
