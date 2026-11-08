#!/bin/bash
# 변환 검사 v2 — 양쪽 가닥 서명을 다 본다
#   지난 검사는 R1 의 C 비율만 봤다. 그런데 바이설파이트 read 는 가닥에 따라 모양이 다르다:
#     C→T 가닥 : C 가 적고 T 가 많다   (지난 검사가 잡던 것)
#     G→A 가닥 : G 가 적고 A 가 많다 (C 는 정상 20% 로 보인다)
#   비방향성·PBAT 라이브러리이거나 R1/R2 가 바뀌어 올라와 있으면
#   변환이 잘 된 자료를 "미변환" 으로 거짓 탈락시킨다.
#   PRJNA734348 난모세포는 이 프로젝트에서 유일한 난모세포 바이설파이트 자료라 반드시 확인한다.
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
        if(100*c/n<=5) ct++          # C→T 가닥 서명
        if(100*g/n<=5) ga++          # G→A 가닥 서명
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
