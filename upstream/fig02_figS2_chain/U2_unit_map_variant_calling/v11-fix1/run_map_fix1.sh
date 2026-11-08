#!/bin/bash
# v11-fix1 map, same code and settings as rung8v10.sh. Only CMPD (=cmp4d_v11_fix1) and WDIR (=g8canonv11f1) differ
set -u
cd ~ || exit 1
L=~/hap_seq_node_v8.tsv
[ -s "$L" ] || { echo "없음 $L"; exit 1; }
[ -d ~/cmp4d_v11_fix1 ] || { echo "cmp4d_v11_fix1 없음"; exit 1; }
[ -e ~/g8canonv11f1 ] && { echo "~/g8canonv11f1 이미 있음 - 중단"; exit 1; }
mkdir -p ~/g8canonv11f1
rm -f ~/g8f1_part_*.txt
awk -v n=4 -v p="$HOME/g8f1_part_" 'NF>=3 && $1!="sample"{print > (p (NR%n) ".txt")}' "$L"
for k in 0 1 2 3; do
  f=~/g8f1_part_"$k".txt
  HSN="$L" CMPD="$HOME/cmp4d_v11_fix1" WDIR="$HOME/g8canonv11f1" TAG="full_p$k" \
    PRESET=asm20 TH=8 LO=2100 HI=2400 FRAG=5 PAD=100 CTX=200 PYTHONIOENCODING=utf-8 \
    setsid nohup python3 ~/consmapG8v4_fixdedup.py "$f" \
    > ~/g8canonv11f1/full_p$k.out 2> ~/g8canonv11f1/full_p$k.err < /dev/null &
done
sleep 5
echo "  기동 $(pgrep -f 'python3.*consmapG8v4_fixdedup[.]py' | wc -l) 개"
