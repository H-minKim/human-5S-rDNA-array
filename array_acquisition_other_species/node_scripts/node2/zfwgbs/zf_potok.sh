#!/bin/bash
# 알의 생물학적 반복 — Potok 2013(PRJNA193544) 계열을 정렬한다
#   우리가 쓴 Jiang 2013(PRJNA188516)과 같은 호 Cell 에 실린 자매 논문이고,
#   알 BioSample 이 SAMN01991170 로 다르다(우리 것 SAMN02117044).
#   → 다른 연구실·다른 알 풀·다른 라이브러리 = 진짜 생물학적 반복.
#   같은 연구 안에 정자·2~16세포·64·256·sphere·근육 이 다 있어 계열째 대조가 된다.
#
#   변환 검사(zf_qc_eggs2.sh) 통과한 것만 돌린다.
#   8M read 쌍씩 하나하나 · 매번 임시 삭제 · 디스크 25GB 미만이면 중단.
set -u
cd ~/zfwgbs || exit 1
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p sra bam logs tmp
N=8000000

guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP — 디스크 여유 ${F}G"; exit 1; fi
}

one () {
  R=$1; D=$2
  [ -f bam/${R}_done ] && { echo "SKIP $R ($D)"; return; }
  guard
  [ -f sra/${R}_1.fastq.gz ] || fastq-dump -X $N --split-files --gzip -O sra "$R" > logs/${R}.dl.log 2>&1
  [ -f sra/${R}_1.fastq.gz ] || { echo "DL_FAIL $R"; return; }
  # 2013 자료라 단일말단으로 올라온 run 이 섞여 있다 — 짝 파일이 없으면 건너뛴다
  [ -f sra/${R}_2.fastq.gz ] || { echo "SINGLE_END $R ($D) — 건너뜀"; rm -f sra/${R}_*.fastq.gz; return; }
  echo "정렬 시작 $R ($D) $(date -u +%H:%M)"
  bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp -o bam \
    -1 sra/${R}_1.fastq.gz -2 sra/${R}_2.fastq.gz > logs/${R}.align.log 2>&1
  B=$(ls bam/${R}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o bam/${R}.sorted.bam "$B" && samtools index bam/${R}.sorted.bam && rm -f "$B" && touch bam/${R}_done
    echo "OK $R  $(grep 'Mapping efficiency' logs/${R}.align.log | head -1)"
  else
    echo "ALIGN_FAIL $R"; tail -3 logs/${R}.align.log
  fi
  rm -f sra/${R}_*.fastq.gz tmp/* 2>/dev/null
  df -h / | tail -1
}

# 중요한 순서대로 — 알이 본 목적, 근육·정자가 같은 연구 안 대조
one SRR800059 알_Potok2013
one SRR800081 근육성체_Potok2013
one SRR800056 정자_Potok2013
one SRR800077 sphere_Potok2013
one SRR800061 2-16세포_Potok2013
one SRR800066 64세포_Potok2013
one SRR800070 256세포_Potok2013
echo POTOK_DONE
