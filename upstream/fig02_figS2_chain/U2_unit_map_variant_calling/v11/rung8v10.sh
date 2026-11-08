#!/bin/bash
# G8v5 map, canonical v10 (surgical repair at 7 sites), full run, 2026-09-13
#   code : consmapG8v4_fixdedup.py (version with the gene HSP dedup F fixed)
#   canonical : ~/cmp4d_v10  (v8 + 7 break sites repaired as [assembly anchor][read consensus][assembly anchor])
#   output : ~/g8canonv10  (existing folder is not overwritten)
set -u
cd ~ || exit 1
L=~/hap_seq_node_v8.tsv
[ -s "$L" ] || { echo "없음 $L"; exit 1; }
[ -d ~/cmp4d_v10 ] || { echo "cmp4d_v10 없음"; exit 1; }
[ -e ~/g8canonv10 ] && { echo "~/g8canonv10 이미 있음 - 중단"; exit 1; }
mkdir -p ~/g8canonv10
rm -f ~/g8v10_part_*.txt
awk -v n=4 -v p="$HOME/g8v10_part_" 'NF>=3 && $1!="sample"{print > (p (NR%n) ".txt")}' "$L"
echo "  정본 $(ls ~/cmp4d_v10 | wc -l) 파일 · 분할 $(cat ~/g8v10_part_*.txt | wc -l) 행"
for k in 0 1 2 3; do
  f=~/g8v10_part_"$k".txt
  HSN="$L" CMPD="$HOME/cmp4d_v10" WDIR="$HOME/g8canonv10" TAG="full_p$k" \
    PRESET=asm20 TH=8 LO=2100 HI=2400 FRAG=5 PAD=100 CTX=200 PYTHONIOENCODING=utf-8 \
    setsid nohup python3 ~/consmapG8v4_fixdedup.py "$f" \
    > ~/g8canonv10/full_p$k.out 2> ~/g8canonv10/full_p$k.err < /dev/null &
done
sleep 5
echo "  기동 $(pgrep -f 'python3.*consmapG8v4_fixdedup[.]py' | wc -l) 개"
