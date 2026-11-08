#!/bin/bash
cd ~/zfwgbs
mkdir -p qc logs
echo -e 'run\tdesc\tR1_C%\tR1_T%\t변환read%\t판정'
for X in 'SRR837876 알rep2_같은연구' 'ERR10788130 알rep2_PRJEB59047' 'ERR10788127 알_bouncerKO_rep1' 'SRR837872 정자_같은연구' 'SRR837894 1k세포_같은연구' 'ERR10788159 배아3hpf_WGBS'; do
  set -- $X; R=$1; D=$2
  [ -f qc/${R}_1.fastq.gz ] || fastq-dump -X 200000 --split-files --gzip -O qc $R > logs/qc2_${R}.log 2>&1
  F=qc/${R}_1.fastq.gz
  [ -f "$F" ] || { echo -e "${R}\t${D}\t-\t-\t-\t내려받기실패"; continue; }
  zcat "$F" | head -400000 | awk -v r="$R" -v d="$D" 'NR%4==2{n=length($0); tot+=n; c=gsub(/C/,"C"); t=gsub(/T/,"T"); C+=c; T+=t; if(100*c/n<=5) k++; m++} END{cp=100*C/tot; tp=100*T/tot; kp=100*k/m; v=(kp>=70)?"쓸만함":((kp>=40)?"의심":"불가"); printf "%s\t%s\t%.1f\t%.1f\t%.1f\t%s\n", r,d,cp,tp,kp,v}'
done
df -h / | tail -1
echo QC2_DONE
