#!/bin/bash
# Stream the URLs in ao_urls.tsv to build the hits per individual, without writing the source files to disk.
#   Step 1: read extraction, HPRC ONT.
# Process one individual's BAMs in parallel, so that the last individuals finish sooner.
#
#  Measured:
#    one stream is capped at 14-25 MB/s, so processing an individual's BAMs one after another
#    takes a median of 251 minutes per individual.
#    Individuals with many BAMs (NA21102 with 10 and NA20503 with 6) add 4 hours to the total.
#    Splitting the BAMs across parallel streams removes that.
#
#  The output is `~/s5s_new/hits/<sample>.hits.fastq.gz` + `~/s5s_new/done/<sample>.done`.
#    The extractor (campaign_extract_v5.py) and its options must be the same for every individual; altering them would break comparability between individuals.
#  Preserving MM/ML with -T MM,ML and passing -F 0x900 are both required.
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
  # cap on the number of concurrent streams
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 5; done
done
wait
echo "[$(date +%T)] $S : BAM 처리 끝 ($(wc -l < "$W/tmp/fx_$S/done.list" 2>/dev/null || echo 0)/$N)"

# the chunks are concatenated in numerical order; the order does not matter for correctness, but it makes the output reproducible
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
