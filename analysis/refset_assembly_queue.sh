#!/bin/bash
# hifiasm --ont 로 5S array 재조립 (노드 큐)
#   단계 : 조립
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md
# hifiasm --ont 로 5S array 를 재조립한다 (1kGP · GIAB · HG008 · 암)
#
#  검증된 전제 (2026-09-02 실측):
#    HPRC 개체 3명 5 하플로타입에서 우리 hifiasm --ont 조립본이
#    HPRC 공식 조립본과 bp 단위로 일치했다 (200410/88 · 117677/51 · 126824/55 ·
#    167029/73 · 185018 vs 185024/81). 방법 자체는 유효하다.
#
#  디스크 방어: hifiasm 한 개체가 676 MB 를 쓴다. 노드에 43 GB 밖에 없으므로
#    개체마다 끝나는 즉시 중간산물을 지우고 array 서열과 요약만 남긴다.
#    (OOM/디스크 부족으로 조용히 손상된 전례가 있어 반드시 확인하며 돈다)
#
#  사용: hifi_queue.sh <목록파일> <출력태그> <스레드>
set -u
LIST=$1
TAG=$2
TH=${3:-6}
OUT=~/hiasm
mkdir -p $OUT/arr $OUT/log
SUM=$OUT/summary_$TAG.tsv

find_reads () {
  local s=$1
  for p in ~/kgp_ont/out/$s.array.fastq.gz \
           ~/s5s_ont/hits/$s.hits.fastq.gz \
           ~/giab_hits/$s.hits.fastq.gz \
           ~/cancer5s/hits/$s.hits.fastq.gz ; do
    [ -f "$p" ] && { echo "$p"; return; }
  done
  echo ""
}

for S in $(tr -d '\r' < "$LIST"); do
  # 이미 한 개체는 건너뛴다 (재시작 안전)
  grep -q "^$S	" $SUM 2>/dev/null && continue
  Q=$(find_reads "$S")
  if [ -z "$Q" ]; then
    printf "%s\tNO_READS\t-\t-\t-\t-\n" "$S" >> $SUM
    continue
  fi
  # 디스크 여유 확인 — 5 GB 미만이면 멈춘다 (조용한 손상 방지)
  FREE=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$FREE" -lt 5 ]; then
    echo "디스크 부족 ${FREE}G — 중단" >> $OUT/log/$TAG.log
    break
  fi
  D=$(mktemp -d ~/hitmp_XXXXXX)
  ( cd "$D" && timeout 7200 ~/hifiasm --ont -o "$S" -t "$TH" "$Q" > hifiasm.log 2>&1 )
  RC=$?
  NC=$(ls "$D"/*.p_ctg.gfa 2>/dev/null | wc -l)
  if [ "$RC" -ne 0 ] || [ "$NC" -eq 0 ]; then
    printf "%s\tFAIL_rc%s\t-\t-\t-\t-\n" "$S" "$RC" >> $SUM
    tail -3 "$D"/hifiasm.log >> $OUT/log/$TAG.log 2>/dev/null
    rm -rf "$D"
    continue
  fi
  # array 추출 (hicheck 와 같은 방법)
  ~/miniconda3/bin/python ~/hiarr.py "$S" "$D" "$OUT/arr" >> $SUM 2>> $OUT/log/$TAG.log
  rm -rf "$D"
done
echo "QUEUE_DONE $TAG"
