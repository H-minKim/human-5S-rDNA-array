#!/bin/bash
# Conversion check: both strand signatures are examined.
#   Looking only at the C fraction of R1 is not enough, because a bisulfite read looks different depending on the strand:
#     the C-to-T strand is low in C and high in T, which is what a C-only check detects
#     the G-to-A strand is low in G and high in A, and its C fraction looks like a normal 20%
#   With a non-directional or PBAT library, or with R1 and R2 submitted the other way round,
#   properly converted data would be rejected as unconverted.
#   The PRJNA734348 oocytes are the only oocyte bisulfite data available here, so they are checked explicitly.
set -u
cd ~/zfwgbs || exit 1
mkdir -p qc logs
printf "run\t설명\t파일\tC%%\tT%%\tG%%\tA%%\tCT형read%%\tGA형read%%\t변환read합%%\t판정\n"
for X in "SRR14703401 난모세포_ctrl_rep1" "SRR14703403 난모세포_ctrl_rep2" \
         "SRR837875 알_Jiang_대조" "SRR800059 알_Potok_대조"; do
  set -- $X
  R=$1; D=$2
  [ -f qc/${R}_1.fastq.gz ] || fastq-dump -X 200000 --split-files --gzip -O qc "$R" > logs/qcs_${R}.log 2>&1
  for M in 1 2; do
    F=qc/${R}_${M}.fastq.gz
    [ -f "$F" ] || { printf "%s\t%s\tR%s\t-\t-\t-\t-\t-\t-\t-\t없음\n" "$R" "$D" "$M"; continue; }
    zcat "$F" | head -400000 | awk -v r="$R" -v d="$D" -v m="$M" '
      NR%4==2{ n=length($0); tot+=n
        c=gsub(/C/,"C"); t=gsub(/T/,"T"); g=gsub(/G/,"G"); a=gsub(/A/,"A")
        C+=c; T+=t; G+=g; A+=a
        if(100*c/n<=5) ct++          # C→T strand signature
        if(100*g/n<=5) ga++          # G→A strand signature
        k++ }
      END{ if(tot==0){exit}
        cp=100*C/tot; tp=100*T/tot; gp=100*G/tot; ap=100*A/tot
        ctp=100*ct/k; gap=100*ga/k; sum=ctp+gap
        v = (sum>=70) ? "쓸 만함" : ((sum>=40) ? "의심" : "불가")
        printf "%s\t%s\tR%s\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%s\n", r,d,m,cp,tp,gp,ap,ctp,gap,sum,v }'
  done
done
rm -f qc/SRR14703401_*.fastq.gz qc/SRR14703403_*.fastq.gz 2>/dev/null
echo QC_STRAND_DONE
