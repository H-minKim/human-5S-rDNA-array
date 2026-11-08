#!/bin/bash
cd ~/zfwgbs
mkdir -p qc logs
echo -e 'run\tdesc\tR1_C%\tR1_T%\t변환read%\t판정'
for X in 'SRR837878 16세포_ZGA전' 'SRR837895 Germring_ZGA후6hpf' 'SRR837899 정소_조직'; do
  set -- $X; R=$1; D=$2
  [ -f qc/${R}_1.fastq.gz ] || fastq-dump -X 200000 --split-files --gzip -O qc $R > logs/qc3_${R}.log 2>&1
  F=qc/${R}_1.fastq.gz
  [ -f "$F" ] || { echo -e "${R}\t${D}\t-\t-\t-\t실패"; continue; }
  zcat "$F" | head -400000 | awk -v r="$R" -v d="$D" 'NR%4==2{n=length($0); tot+=n; c=gsub(/C/,"C"); t=gsub(/T/,"T"); C+=c; T+=t; if(100*c/n<=5) k++; m++} END{kp=100*k/m; v=(kp>=70)?"쓸만함":((kp>=40)?"의심":"불가"); printf "%s\t%s\t%.1f\t%.1f\t%.1f\t%s\n", r,d,100*C/tot,100*T/tot,kp,v}'
  rm -f qc/${R}_*.fastq.gz
done
echo QC3_DONE
