#!/bin/bash
# 금화조 : 5S array read 만 CPU dorado 재콜링 (5 kHz flowcell 만 — 4 kHz 는 dorado 2.1.2 에 모델 없음)
#   dorado basecaller hac@v5.0.0,5mCG_5hmCG (@v3) · --read-ids · pod5 는 mount-s3 로 필요한 부분만 읽는다
#     CPU 실측 (같은 20 read · 8 코어): hac@v5.0.0 2.64e4 samples/s vs 자동선택 hac@v6.0.0 9.5e3 -> v5 가 2.8 배 빠름
#   -> samtools fastq -T MM,ML -> minimap2 -y -ax map-ont --secondary=no (v0.3 이배체) -> 정렬 BAM (+MM/ML)
#      (제브라피시와 같은 정렬 방식 · MAPQ 로 거르지 않는다)
#   인자 : flowcell 목록 파일 · taskset 코어 범위 (예 0-15) · [read id 폴더 = ~/xsp/zfont/sel] · [모델 = hac,5mCG_5hmCG]
#   산출 : bc/<fc>.ubam (재콜링 원본) · bc/<fc>.bam(+bai) · bc/<fc>.done (종료코드 · read 수 · MM 있는 read · 초)
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
Z=~/xsp/zfont
A=~/xsp/asm
D=~/xsp/tools/dorado-2.1.2-linux-x64/bin/dorado
MNT=~/xsp/ms3mnt
L="${1:?flowcell 목록}"
CORES="${2:?코어 범위}"
IDD="${3:-$Z/sel}"
MODEL="${4:-hac@v5.0.0,5mCG_5hmCG}"
BC="${BC:-$Z/bc}"
mkdir -p "$BC" ~/xsp/dmodels
for b in $(cat "$L"); do
  [ -s $BC/$b.done ] && continue
  p=$(ls $MNT/*/raw/uncorrected/$b.pod5 2> /dev/null | head -1)
  [ -n "$p" ] || { printf '%s\tNO_POD5\n' "$b" > $Z/bc/$b.done; continue; }
  [ -s "$IDD/$b.ids" ] || { printf '%s\tNO_IDS\n' "$b" > $Z/bc/$b.done; continue; }
  s=$(date +%s)
  taskset -c "$CORES" $D basecaller "$MODEL" "$p" --read-ids "$IDD/$b.ids" --device cpu --models-directory ~/xsp/dmodels > $BC/$b.ubam 2> $BC/$b.dorado.log
  rc=$?
  # primary 만 (-F 0x900 · 제브라피시와 같게 · 보조정렬 섞임 방지)
  samtools fastq -T MM,ML $BC/$b.ubam 2> /dev/null | minimap2 -y -ax map-ont -t 6 --secondary=no $A/zfinch_v03_dip.map-ont.mmi - 2> /dev/null | samtools view -b -F 0x900 - | samtools sort -@2 -o $BC/$b.bam - && samtools index $BC/$b.bam
  n=$(samtools view -c $BC/$b.ubam)
  m=$(samtools view $BC/$b.ubam | grep -c "MM:Z")
  printf '%s\t%s\t%s\t%s\t%s\n' "$b" "$rc" "$n" "$m" "$(( $(date +%s) - s ))" > $BC/$b.done
done
echo "재콜링 끝 $(date -u +%H:%M) UTC"
cat $Z/bc/*.done
