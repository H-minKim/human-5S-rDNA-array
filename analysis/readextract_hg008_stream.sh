#!/bin/bash
# s3://giab/data_somatic/HG008/Liss_lab/Northeastern_ONT-std_20240422/
#   단계 : ① read 추출 (HG008)
# HG008 3샘플 — 전량 스트리밍으로 통일 (길이필터본 금지: T 와 N 의 조건이 달라져 편향)
#   T   : CHM13v2.0 정렬 BAM 149GB
#   N-D : ubam 4 개 (90GB급)
#   N-P : ubam 3 개 (46-50GB급)
set -u
LOG=~/hg008_queue.log
echo "[$(date -u +%F\ %H:%M)] v2 대기 시작 — 전량 스트리밍 통일" >> "$LOG"
while :; do
  n=$(ls ~/s5s_ont/done 2>/dev/null | wc -l)
  r=$(pgrep -f "campaign_run_ont_v[5]" 2>/dev/null | wc -l)
  g=$(pgrep -f "giab_queu[e]" 2>/dev/null | wc -l)
  if { [ "$n" -ge 264 ] || { [ "$r" -eq 0 ] && [ "$n" -ge 100 ]; }; } && [ "$g" -le 1 ]; then
    echo "[$(date -u +%F\ %H:%M)] 시작 (HPRC $n · 스트림 $r)" >> "$LOG"; break
  fi
  sleep 600
done
W=~/hg008; mkdir -p $W/{out,log,done}; cd $W || exit 9
BB=$HOME/opt/bbmap/bbduk.sh; [ -x "$BB" ] || BB=$HOME/s5s/opt/bbmap/bbduk.sh
cp -f ~/s5s_ont/{5S_unit_discseg.fa,5S_gene_query.fa,campaign_extract_v5.py,tandem40.fa,summarize4.py} . 2>/dev/null
cp -f ~/s5s_ont/gene_db.* . 2>/dev/null
BASE="https://giab.s3.amazonaws.com/data_somatic/HG008/Liss_lab"
DIR="Northeastern_ONT-std_20240422"
one(){ # $1=라벨 $2=URL
  L=$1; U=$2
  [ -f "done/${L}.done" ] && return
  while [ "$(df -BG / | tail -1 | awk '{print $4}' | tr -d G)" -lt 40 ]; do sleep 300; done
  t0=$(date +%s)
  set -o pipefail
  samtools view -h -@2 -F 0x900 "$U" 2> "log/${L}.src.err" \
    | samtools fastq -T MM,ML -@2 - 2>/dev/null \
    | python3 campaign_extract_v5.py --nts 5S_unit_discseg.fa --gene 5S_gene_query.fa \
        --gene-db gene_db --bbduk "$BB" --out "out/${L}" --chunk 16000 --threads 6 \
        --k 25 --hdist 2 --tmp "$W" > "log/${L}.log" 2>&1
  st=("${PIPESTATUS[@]}"); set +o pipefail
  bad=""; for j in "${!st[@]}"; do [ "${st[$j]}" -ne 0 ] && bad="$bad s$j=${st[$j]}"; done
  t1=$(date +%s)
  nr=$(( $(wc -l < "out/${L}.hits.fastq" 2>/dev/null || echo 0) / 4 ))
  if [ -n "$bad" ]; then echo "[$(date -u +%H:%M)] ${L} 실패:$bad" >> "$LOG"
  else
    pigz -p4 -f "out/${L}.hits.fastq"
    echo "${L} ok reads=$nr sec=$((t1-t0))" > "done/${L}.done"
    echo "[$(date -u +%H:%M)] ${L} read=$nr $((t1-t0))s" >> "$LOG"
  fi
}
one "HG008-T" "$BASE/Northeastern-ONT-UL-20241216/HG008-T_CHM13v2.0_ONT-UL-R10.4.1-dorado_0.8.1_sup.5mC_5hmC_54x_20241216.bam"
for k in HG008-N-D HG008-N-P; do
  i=0
  for f in $(timeout 120 s5cmd --no-sign-request ls "s3://giab/data_somatic/HG008/Liss_lab/${DIR}/${k}_ubams/" 2>/dev/null | awk '{print $NF}' | grep "\.bam$"); do
    i=$((i+1)); one "${k}_${i}" "$BASE/${DIR}/${k}_ubams/${f}"
  done
done
echo "[$(date -u +%F\ %H:%M)] HG008_DONE 산출 $(ls out/*.gz 2>/dev/null | wc -l)" >> "$LOG"
