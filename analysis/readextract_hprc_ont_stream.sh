#!/bin/bash
# ~/ao_urls.tsv 의 URL 을 흘려 읽어 개체별 hits 를 만든다 (디스크에 원본을 받지 않는다)
#   단계 : ① read 추출 (HPRC ONT)
# 한 개체의 BAM 을 병렬로 처리해 마지막 개체를 빨리 끝낸다 (2026-09-04)
#
#  왜 필요한가 (실측)
#    extract72.sh 는 한 개체의 BAM 을 `for u in $urls` 로 순차 처리한다.
#    그래서 개체당 스트림이 1 개뿐이고 14~25 MB/s 에 묶여, 개체당 중앙 251 분이 걸린다.
#    남은 2 개체(NA21102 10 BAM · NA20503 6 BAM)가 전체 완료를 4 시간 늘리고 있다.
#    BAM 을 병렬로 쪼개면 그만큼 줄어든다.
#
#  출력은 extract72.sh 와 같은 자리·같은 형식이어야 한다
#    (`~/s5s_new/hits/<sample>.hits.fastq.gz` + `~/s5s_new/done/<sample>.done`).
#    추출기(campaign_extract_v5.py)와 옵션은 그대로 쓴다 — 바꾸면 264 개체와 비교가 깨진다.
#  MM/ML 보존(-T MM,ML)과 -F 0x900 은 필수.
set -u
S=${1:?개체}
PAR=${2:-5}
W=~/s5s_new
G=~/giab5s_seq
BB=~/s5s/opt/bbmap/bbduk.sh
MAN=~/ao_urls.tsv
mkdir -p $W/hits $W/done $W/log $W/tmp/fx_$S

[ -f "$W/done/$S.done" ] && { echo "  $S 이미 완료"; exit 0; }
urls=$(awk -F'\t' -v s="$S" '$1==s{print $2}' "$MAN")
N=$(echo "$urls" | grep -c .)
[ "$N" -gt 0 ] || { echo "  $S URL 없음"; exit 1; }
echo "[$(date +%T)] $S : BAM $N 개를 최대 $PAR 병렬로"

i=0
for u in $urls; do
  i=$((i+1))
  (
    curl -s -L --retry 3 --retry-delay 5 "$u" 2>/dev/null \
      | samtools fastq -@ 2 -F 0x900 -T MM,ML - 2>/dev/null \
      | python3 $G/campaign_extract_v5.py \
          --nts $G/5S_unit_discseg.fa --gene $G/5S_gene_query.fa \
          --gene-db $G/gene_db --bbduk $BB \
          --out "$W/tmp/fx_$S/p$i" --threads 3 --tmp "$W/tmp/fx_$S" \
          >> "$W/log/$S.fx.log" 2>&1
    echo "$i" >> "$W/tmp/fx_$S/done.list"
  ) &
  # 동시 수 제한
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 5; done
done
wait
echo "[$(date +%T)] $S : BAM 처리 끝 ($(wc -l < "$W/tmp/fx_$S/done.list" 2>/dev/null || echo 0)/$N)"

# 조각을 번호 순서대로 이어붙인다 (순서가 섞이면 안 되는 건 아니지만 재현성을 위해)
: > "$W/tmp/$S.fx.fastq"
for k in $(seq 1 $N); do
  f="$W/tmp/fx_$S/p$k.hits.fastq"
  [ -s "$f" ] && cat "$f" >> "$W/tmp/$S.fx.fastq"
done
if [ -s "$W/tmp/$S.fx.fastq" ]; then
  gzip -c "$W/tmp/$S.fx.fastq" > "$W/hits/$S.hits.fastq.gz"
  R=$(( $(wc -l < "$W/tmp/$S.fx.fastq") / 4 ))
  echo "$R" > "$W/done/$S.done"
  echo "[$(date +%T)] $S 완료 · 5S read $R · $(ls -lh $W/hits/$S.hits.fastq.gz | awk '{print $5}')" >> $W/progress.log
  echo "  $S 완료 · read $R"
else
  echo "0" > "$W/done/$S.done"
  echo "[$(date +%T)] $S 결과 없음" >> $W/progress.log
  echo "  $S 결과 없음"
fi
rm -rf "$W/tmp/fx_$S" "$W/tmp/$S.fx.fastq"
