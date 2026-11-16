#!/bin/bash
# A screening check on bisulfite data: conversion alone is tested on a small subset before the main run.
#   Measured: the PRJNA734348 oocyte samples are registered as Bisulfite-Seq, but
#   78-86% of their reads are in fact unconverted genomic DNA (mapping rate 1.2-2.1%).
#   So for each candidate only 200,000 spots are fetched and conversion is judged first from the C fraction.
#     A properly converted read carries 2-5% C, most having become T.
#     An unconverted read carries about 20% C, the normal genomic composition.
cd ~/zfwgbs
mkdir -p qc logs
echo -e "run\tdesc\tR1_C_pct\tR1_T_pct\tconverted_read_pct\t판정"
for X in "SRR837875 egg_Jiang2013" "ERR10788129 egg_PRJEB59047" "SRR19364705 control_egg_PRJNA841653" \
         "SRR11614917 embryo_24hpf_EMseq" "SRR12797061 intestine_WGBS" "SRR12797059 heart_WGBS"; do
  set -- $X
  R=$1; D=$2
  if [ ! -f qc/${R}_1.fastq.gz ]; then
    fastq-dump -X 200000 --split-files --gzip -O qc "$R" > logs/qc_${R}.log 2>&1
  fi
  F=qc/${R}_1.fastq.gz
  [ -f "$F" ] || { echo -e "${R}\t${D}\t-\t-\t-\t내려받기 실패"; continue; }
  zcat "$F" | head -400000 | awk -v r="$R" -v d="$D" '
    NR%4==2{ n=length($0); tot+=n;
      c=gsub(/C/,"C"); t=gsub(/T/,"T"); C+=c; T+=t;
      if(100*c/n<=5) k++; m++ }
    END{ cp=100*C/tot; tp=100*T/tot; kp=100*k/m;
      v = (kp>=70) ? "쓸 만함" : ((kp>=40) ? "의심" : "불가(변환 안 됨)");
      printf "%s\t%s\t%.1f\t%.1f\t%.1f\t%s\n", r, d, cp, tp, kp, v }'
done
df -h / | tail -1
echo QC_DONE
