#!/bin/bash
# Candidate biological replicates for the eggs: a conversion check on 200,000 spots.
#   PRJNA193544 is Potok 2013 (Cell), the companion paper to PRJNA188516 (Jiang 2013)
#   published alongside it in the same issue, and it holds the whole series: egg, sperm, 2-16 cell, 64, 256, sphere and muscle.
#   Its egg BioSample is SAMN01991170, different from SAMN02117044, so it is a candidate genuine biological replicate.
#   ERR10788130 (PRJEB59047 Egg_rep2, a different BioSample) and
#   SRR19364710 (PRJNA841653 control eggs 1) are checked as well.
#
#   The criterion: a properly converted read carries 2-5% C, most having become T, while an unconverted read carries about 20% C.
set -u
cd ~/zfwgbs || exit 1
mkdir -p qc logs
printf "run\tdesc\tR1_C_pct\tR1_T_pct\tconverted_read_pct\tverdict\n"
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
  [ -f "$F" ] || { printf "%s\t%s\t-\t-\t-\tdownload failed\n" "$R" "$D"; continue; }
  zcat "$F" | head -400000 | awk -v r="$R" -v d="$D" '
    NR%4==2{ n=length($0); tot+=n;
      c=gsub(/C/,"C"); t=gsub(/T/,"T"); C+=c; T+=t;
      if(100*c/n<=5) k++; m++ }
    END{ if(tot==0){printf "%s\t%s\t-\t-\t-\tempty file\n", r, d; exit}
      cp=100*C/tot; tp=100*T/tot; kp=100*k/m;
      v = (kp>=70) ? "usable" : ((kp>=40) ? "doubtful" : "unusable(not converted)");
      printf "%s\t%s\t%.1f\t%.1f\t%.1f\t%s\n", r, d, cp, tp, kp, v }'
done
rm -f qc/*.fastq.gz
df -h / | tail -1
echo QC2_DONE
