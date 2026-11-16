#!/bin/bash
# K562 UL2 assembly (hifiasm --ul: 24,103 ONT 5S reads of 25 runs + HiFi array-region reads) -> K562ul2 hap1/hap2 contigs
# used for Supplementary Fig 2O.
# hifiasm 0.25.0-r726.
set -u
W=~/p11_encode/k562ul2; mkdir -p $W; cd $W
rm -f *.bin *.gfa *.fa UL2.done
~/hifiasm -o K562ul2 -t 14 --ul ~/p11_encode/k562ont/K562_ONT_ALL.hits_5S.fastq.gz ~/p11_encode/k562hifi_dl/k562_region.fastq.gz > hifiasm.log 2>&1
rc=$?
for g in *p_ctg.gfa; do [ -e "$g" ] && awk '/^S/{print ">"$2"\n"$3}' "$g" > "${g%.gfa}.fa"; done
echo "rc=$rc fa=$(ls *p_ctg.fa 2>/dev/null | wc -l)" > UL2.done
