#!/bin/bash
# cancer17.tsv 의 URL 을 흘려 읽는다 (BAM·fastq.gz 두 종류)
#   단계 : ① read 추출 (암 17 시료)
# 암 17 샘플에서 5S rDNA read 추출 — HPRC/GIAB 와 동일한 v5 관문(BBDuk k25 hdist2)
# 원본을 디스크에 받지 않는다(compute-to-data 스트리밍).
# 입력이 BAM(PacBio/GIAB) 과 fastq.gz(ENA CASTLE) 두 종류다.
set -u
W=~/cancer5s; mkdir -p $W/{hits,done,log,tmp}; cd $W || exit 9
G=~/giab5s_seq
BB=~/s5s/opt/bbmap/bbduk.sh
MAN=~/spanfork/cancer17.tsv
for f in $G/campaign_extract_v5.py $G/5S_unit_discseg.fa $G/5S_gene_query.fa $G/gene_db.nin $BB; do
  [ -e "$f" ] || { echo "FATAL 없음 $f"; exit 9; }
done
: > $W/progress.log
for S in $(cut -f1 $MAN | sort -u); do
  [ -f $W/done/$S.done ] && continue
  mkdir $W/done/$S.lock 2>/dev/null || continue
  urls=$(awk -F'\t' -v s="$S" '$1==s{print $2}' $MAN)
  [ -z "$urls" ] && { rmdir $W/done/$S.lock; continue; }
  echo "[$(date +%T)] $S 시작 ($(echo "$urls"|wc -l) URL)" >> $W/progress.log
  rm -f $W/tmp/$S.hits.fastq
  ok=1
  for u in $urls; do
    case "$u" in
      *.bam)      SRC="curl -s -L --retry 3 '$u' | samtools fastq -@ 4 -F 0x900 -T MM,ML -" ;;
      *.fastq.gz) SRC="curl -s -L --retry 3 '$u' | zcat" ;;
      *)          SRC="curl -s -L --retry 3 '$u'" ;;
    esac
    eval "$SRC" 2>/dev/null | python3 $G/campaign_extract_v5.py \
      --nts $G/5S_unit_discseg.fa --gene $G/5S_gene_query.fa --gene-db $G/gene_db \
      --bbduk $BB --out $W/tmp/${S}_part --threads 6 --tmp $W/tmp >> $W/log/$S.log 2>&1
    rc=$?
    [ $rc -ne 0 ] && { echo "[$(date +%T)] $S URL 실패 rc=$rc" >> $W/progress.log; ok=0; }
    [ -s $W/tmp/${S}_part.hits.fastq ] && cat $W/tmp/${S}_part.hits.fastq >> $W/tmp/$S.hits.fastq
    rm -f $W/tmp/${S}_part.hits.fastq $W/tmp/${S}_part.5s.tsv
  done
  if [ -s $W/tmp/$S.hits.fastq ]; then
    gzip -c $W/tmp/$S.hits.fastq > $W/hits/$S.hits.fastq.gz
    n=$(( $(wc -l < $W/tmp/$S.hits.fastq) / 4 ))
    echo "[$(date +%T)] $S 완료 · 5S read $n · $(ls -lh $W/hits/$S.hits.fastq.gz|awk '{print $5}')" >> $W/progress.log
  else
    echo "[$(date +%T)] $S 5S read 0" >> $W/progress.log
  fi
  rm -f $W/tmp/$S.hits.fastq
  touch $W/done/$S.done; rmdir $W/done/$S.lock 2>/dev/null
done
touch $W/C5S_DONE
