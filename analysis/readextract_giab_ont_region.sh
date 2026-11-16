#!/bin/bash
# Extract from s3://ont-open-data/giab_2023.05/analysis/<sample>/sup/*.pass.cram by region query over GRCh38 chr1:228,460,268-228,796,040.
#   Step 1: read extraction, GIAB ONT.
# Extraction of the GIAB trio ONT data, by region query against the indexed GRCh38-aligned CRAM.
# The GRCh38 5S array is chr1:228,610,268-228,646,040, holding 17 5S genes, taken with +-150 kb of slack.
set -u
W=~/giab5s; mkdir -p $W/{out,log,done}; cd $W || exit 9
REG="chr1:228460268-228796040"
BASE="https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis"
declare -A FC
FC[hg001]="PAO83395 PAO89685"
FC[hg002]="PAO83395 PAO89685"
FC[hg003]="PAO83395 PAO89685"
FC[hg004]="PAO83395 PAO89685"
for S in hg001 hg002 hg003 hg004; do
  # the actual flowcell list is queried and used
  L=$(timeout 90 s5cmd --no-sign-request ls "s3://ont-open-data/giab_2023.05/analysis/${S}/sup/" 2>/dev/null | awk '{print $NF}' | grep "pass.cram$" | head -4)
  [ -n "$L" ] || { echo "  $S CRAM 없음"; continue; }
  n=0
  for f in $L; do
    n=$((n+1))
    O="out/${S}_${n}.fastq.gz"
    [ -s "$O" ] && continue
    U="${BASE}/${S}/sup/${f}"
    t0=$(date +%s)
    set -o pipefail
    timeout 900 samtools view -@2 -F 0x900 -u "$U" "$REG" 2> "log/${S}_${n}.err" \
      | samtools fastq -T MM,ML,HP,PS -@2 - 2>/dev/null | pigz -p4 -c > "$O"
    st=("${PIPESTATUS[@]}"); set +o pipefail
    bad=""; for j in "${!st[@]}"; do [ "${st[$j]}" -ne 0 ] && bad="$bad s$j=${st[$j]}"; done
    t1=$(date +%s)
    nr=$(zcat "$O" 2>/dev/null | awk 'END{print int(NR/4)}')
    mm=$(zcat "$O" 2>/dev/null | awk 'NR%4==1 && /MM:Z/{m++} END{print m+0}')
    printf "  %-6s %-24s read %6s · MM %6s · %ss%s\n" "$S" "$f" "${nr:-0}" "${mm:-0}" "$((t1-t0))" "$bad"
    [ -n "$bad" ] && rm -f "$O"
  done
done
echo "GIAB_DONE 산출 $(ls out/*.fastq.gz 2>/dev/null | wc -l)"
