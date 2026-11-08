#!/bin/bash
# 바이설파이트 자료 선별 검사 — 본 작업 전에 소량으로 변환 여부만 본다
#   실패 경험 : PRJNA734348 난자 시료는 'Bisulfite-Seq' 로 등록돼 있으나
#   실제로는 read 의 78~86% 가 변환되지 않은 유전체 DNA 였다(매핑률 1.2~2.1%).
#   그래서 후보마다 20 만 spot 만 받아 C 비율로 변환 여부를 먼저 판정한다.
#     제대로 변환된 read : C 가 2~5% (대부분 T 로 바뀜)
#     변환 안 된 read    : C 가 20% 안팎(정상 유전체 조성)
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
