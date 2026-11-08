#!/bin/bash
# 알의 생물학적 반복 후보 — 20만 spot 변환 검사
#   PRJNA193544 = Potok 2013(Cell). 우리가 쓴 PRJNA188516(Jiang 2013)과
#   같은 호에 나란히 실린 자매 논문이고, 알·정자·2~16세포·64·256·sphere·근육 계열이 통째로 있다.
#   알 BioSample 이 SAMN01991170 로 우리 것(SAMN02117044)과 다르다 → 진짜 생물학적 반복 후보.
#   함께 ERR10788130(PRJEB59047 Egg_rep2·다른 BioSample)과
#   SRR19364710(PRJNA841653 control eggs 1)도 마지막으로 확인한다.
#
#   판정 : 제대로 변환된 read 는 C 2~5%(대부분 T 로 바뀜) / 안 된 read 는 C 20% 안팎.
set -u
cd ~/zfwgbs || exit 1
mkdir -p qc logs
printf "run\tdesc\tR1_C_pct\tR1_T_pct\tconverted_read_pct\t판정\n"
for X in "SRR800059 알_Potok2013_PRJNA193544" \
         "SRR800056 정자_Potok2013" \
         "SRR800061 2-16세포_Potok2013" \
         "SRR800066 64세포_Potok2013" \
         "SRR800070 256세포_Potok2013" \
         "SRR800077 sphere_Potok2013" \
         "SRR800081 근육성체_Potok2013" \
         "ERR10788130 알rep2_PRJEB59047" \
         "SRR19364710 control_eggs1_PRJNA841653"; do
  set -- $X
  R=$1; D=$2
  if [ ! -f qc/${R}_1.fastq.gz ]; then
    fastq-dump -X 200000 --split-files --gzip -O qc "$R" > logs/qc_${R}.log 2>&1
  fi
  F=qc/${R}_1.fastq.gz
  [ -f "$F" ] || { printf "%s\t%s\t-\t-\t-\t내려받기 실패\n" "$R" "$D"; continue; }
  zcat "$F" | head -400000 | awk -v r="$R" -v d="$D" '
    NR%4==2{ n=length($0); tot+=n;
      c=gsub(/C/,"C"); t=gsub(/T/,"T"); C+=c; T+=t;
      if(100*c/n<=5) k++; m++ }
    END{ if(tot==0){printf "%s\t%s\t-\t-\t-\t빈 파일\n", r, d; exit}
      cp=100*C/tot; tp=100*T/tot; kp=100*k/m;
      v = (kp>=70) ? "쓸 만함" : ((kp>=40) ? "의심" : "불가(변환 안 됨)");
      printf "%s\t%s\t%.1f\t%.1f\t%.1f\t%s\n", r, d, cp, tp, kp, v }'
done
rm -f qc/*.fastq.gz
df -h / | tail -1
echo QC2_DONE
