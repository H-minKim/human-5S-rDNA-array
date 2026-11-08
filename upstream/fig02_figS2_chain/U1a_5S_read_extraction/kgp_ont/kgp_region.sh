#!/bin/bash
# 1kGP-ONT: query **only the regions** from the CHM13-aligned BAM (the whole BAM is not downloaded).
#   HPRC ONT is raw BAM, so full streaming (141GB per individual) was unavoidable,
#   but this one has a .bai index, so only the needed regions are fetched by range request. Measured 12 s per individual.
#
# Regions to extract
#   ① 5S array (chr1:227,744,873-228,024,958) + flank 100kb
#   (2) top ectopic loci (those with read>=100 in loci.tsv)
# Exclude secondary/supplementary (-F 0x900): aligned BAMs contain duplicate records
#   (precedent: in 45S, 81% of the records in the array region were secondary)
set -u
LIST="${1:?sample list}"
W=~/kgp_ont
cd "$W" || exit 9
mkdir -p out log done
ARRAY="chr1:227644873-228124958"      # array ±100kb
REG="regions.txt"
[ -s "$REG" ] || echo "$ARRAY" > "$REG"

while read -r S; do
  [ -n "$S" ] || continue
  [ -f "done/${S}.done" ] && { echo "  skip $S"; continue; }
  U=$(grep -P "^${S}\t" manifest.tsv | cut -f3)
  [ -n "$U" ] || { echo "  $S URL 없음"; echo "$S" >> log/no_url.txt; continue; }
  H=$(echo "$U" | sed 's|s3://1000g-ont/|https://1000g-ont.s3.amazonaws.com/|')
  t0=$(date +%s)
  # Index query -> primary only -> fastq (MM/ML/HP kept)
  set -o pipefail
  samtools view -@2 -F 0x900 -u "$H" $(cat "$REG" | tr '\n' ' ') 2> "log/${S}.err" \
    | samtools fastq -T MM,ML,HP -@2 - 2>/dev/null \
    | pigz -p4 -c > "out/${S}.array.fastq.gz"
  st=("${PIPESTATUS[@]}"); set +o pipefail
  bad=""; for j in "${!st[@]}"; do [ "${st[$j]}" -ne 0 ] && bad="$bad s$j=${st[$j]}"; done
  if [ -n "$bad" ]; then
    echo "  $S 파이프 실패:$bad"; echo "$S$bad" >> log/failed.txt
    rm -f "out/${S}.array.fastq.gz"; continue
  fi
  gz=$(stat -c%s "out/${S}.array.fastq.gz" 2>/dev/null || echo 0)
  [ "$gz" -gt 1000 ] || { echo "  $S 산출 없음"; echo "$S" >> log/empty.txt; continue; }
  nr=$(zcat "out/${S}.array.fastq.gz" | awk 'END{print NR/4}')
  nm=$(zcat "out/${S}.array.fastq.gz" | awk 'NR%4==1 && /MM:Z:/{n++} END{print n+0}')
  t1=$(date +%s)
  echo "$S ok reads=$nr mm=$nm gz=$gz sec=$((t1-t0))" > "done/${S}.done"
  echo "  [$(date -u +%H:%M:%S)] $S read=$nr MM=$nm $((gz/1048576))MB $((t1-t0))s"
done < "$LIST"
echo "[$(date -u +%H:%M:%S)] KGP_DONE 완료=$(ls done | wc -l)"
