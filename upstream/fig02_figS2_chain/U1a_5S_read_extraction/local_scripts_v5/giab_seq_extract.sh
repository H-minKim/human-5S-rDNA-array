#!/bin/bash
# Stream the full CRAMs of the 4 GIAB individuals -> BBDuk NTS gate (same method as HPRC v5)
set -u
W=~/giab5s_seq; mkdir -p $W/{out,log,done}; cd $W || exit 9
export PATH=$HOME/s5s_ont:$PATH
BB=~/opt/bbmap/bbduk.sh; [ -x "$BB" ] || BB=~/s5s/opt/bbmap/bbduk.sh
for f in ~/s5s_ont/5S_unit_discseg.fa ~/s5s_ont/5S_gene_query.fa ~/s5s_ont/gene_db.nin ~/s5s_ont/campaign_extract_v5.py; do
  [ -s "$f" ] || { echo "FATAL 없음 $f"; exit 9; }
done
for S in hg001 hg002 hg003 hg004; do
  L=$(timeout 120 s5cmd --no-sign-request ls "s3://ont-open-data/giab_2023.05/analysis/${S}/sup/" 2>/dev/null | awk '{print $NF}' | grep "pass.cram$")
  n=0
  for f in $L; do
    n=$((n+1))
    [ -f "done/${S}_${n}.done" ] && continue
    U="https://ont-open-data.s3.amazonaws.com/giab_2023.05/analysis/${S}/sup/${f}"
    t0=$(date +%s)
    set -o pipefail
    samtools view -@2 -u "$U" 2> "log/${S}_${n}.err" \
      | samtools fastq -T MM,ML -@2 - 2>/dev/null \
      | python3 ~/s5s_ont/campaign_extract_v5.py --nts ~/s5s_ont/5S_unit_discseg.fa \
          --gene ~/s5s_ont/5S_gene_query.fa --gene-db ~/s5s_ont/gene_db --bbduk "$BB" \
          --out "out/${S}_${n}" --chunk 16000 --threads 6 --k 25 --hdist 2 --tmp $W \
          > "log/${S}_${n}.log" 2>&1
    st=("${PIPESTATUS[@]}"); set +o pipefail
    bad=""; for j in "${!st[@]}"; do [ "${st[$j]}" -ne 0 ] && bad="$bad s$j=${st[$j]}"; done
    t1=$(date +%s)
    nr=$(awk 'END{print int(NR/4)}' "out/${S}_${n}.hits.fastq" 2>/dev/null)
    echo "  $S $f read=${nr:-0} $((t1-t0))s$bad"
    [ -z "$bad" ] && echo "ok" > "done/${S}_${n}.done"
  done
done
echo "GIABSEQ_DONE"
