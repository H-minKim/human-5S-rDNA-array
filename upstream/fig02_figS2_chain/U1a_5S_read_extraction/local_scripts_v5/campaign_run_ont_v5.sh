#!/bin/bash
# ONT campaign runner v5, BBDuk single gate plus TUC tandem map (CHM13 alignment and tRNA removed)
#
# Basis (all measured)
#  (1) Gate replacement: taking the 2,262 reads that the current pipeline called ARRAY and kept as the truth,
#     BBDuk k31 hdist1 (NTS) 188MB/s recall 1.0000, blastn ws16 8MB/s recall 1.0000 (first measurement used to choose the gate, one 383MB input. The value actually run is k=25 hdist=2 = line 53, ONT tandem span>500 recall k31h1 0.9158, k25h2 1.0000)
#     and megablast 14MB/s recall 0.9881. -> 23 times faster with no loss of recall.
#     End-to-end check on the full HG00117 BAM: of 283 current ARRAY reads, v5 recovers 283, misses 0.
#  (2) CHM13 alignment (44-54% of per-sample time) replaced by the tandem map. Hits are only 0.036%, so
#     the downstream shrinks from 37GB to 360MB and becomes practically free.
#  (3) tRNA excluded. BBDuk runs on the 5S discriminating segment, so no separate step is needed.
#  (4) All hit reads are kept, about 360MB per sample (gz about 150MB), 34GB for 228 samples.
#     The old pipeline kept only ARRAY, but v5 is already selective enough that keeping everything is feasible.
#
# Usage: bash campaign_run_ont_v5.sh <sample_list>
set -u
LIST="${1:?sample list}"
W=~/s5s_ont
mkdir -p "$W"/{done,sum,array,hits,log,tmp} && cd "$W" || exit 9
for t in blastn samtools minimap2 s5cmd java pigz; do
  command -v $t >/dev/null || { echo "FATAL: $t not on PATH"; exit 9; }
done
for f in 5S_gene_query.fa gene_db.nin 5S_unit_discseg.fa tandem40.fa \
         campaign_extract_v5.py summarize4.py manifest.tsv; do
  [ -s "$f" ] || { echo "FATAL: $f 없음/빈파일"; exit 9; }
done
BB=~/opt/bbmap/bbduk.sh; [ -x "$BB" ] || BB=~/s5s/opt/bbmap/bbduk.sh
[ -x "$BB" ] || { echo "FATAL: bbduk.sh 없음"; exit 9; }
THR=${THR:-6}; CHUNK=${CHUNK:-16000}; MINFREE=${MINFREE:-30}; PIGZ_T=${PIGZ_T:-4}

echo "[$(date -u +%H:%M:%S)] v5 시작 개체=$(wc -l < "$LIST") THR=$THR"

while read -r S; do
  [ -n "$S" ] || continue
  [ -f "done/${S}.done" ] && continue
  while [ "$(df -BG / | tail -1 | awk '{print $4}' | tr -d G)" -lt "$MINFREE" ]; do
    echo "  [$(date -u +%H:%M:%S)] 디스크 부족 — 대기"; sleep 300
  done
  t0=$(date +%s)
  rm -f "tmp/${S}"*
  grep -P "^${S}\t" manifest.tsv | cut -f3 > "tmp/${S}.bams"
  nb=$(wc -l < "tmp/${S}.bams")
  [ "$nb" -gt 0 ] || { echo "  $S BAM 없음"; echo "$S" >> log/no_bam.txt; continue; }

  ok=1; : > "tmp/${S}.hits.fastq"; : > "tmp/${S}.5s.tsv"
  for i in $(seq 1 "$nb"); do
    U=$(sed -n "${i}p" "tmp/${S}.bams")
    set -o pipefail
    s5cmd --no-sign-request cat -c 2 -p 8 "$U" 2> "log/${S}_${i}.s5err" \
      | samtools fastq -T MM,ML -@ 2 - 2>/dev/null \
      | python3 campaign_extract_v5.py --nts 5S_unit_discseg.fa --gene 5S_gene_query.fa \
          --gene-db gene_db --bbduk "$BB" --out "tmp/${S}_${i}" \
          --chunk "$CHUNK" --threads "$THR" --k "${BBK:-25}" --hdist "${BBH:-2}" --tmp "$W/tmp" > "log/${S}_${i}.log" 2>&1
    st=("${PIPESTATUS[@]}")
    set +o pipefail
    bad=""
    for j in "${!st[@]}"; do [ "${st[$j]}" -ne 0 ] && bad="$bad stage$j=${st[$j]}"; done
    if [ -n "$bad" ]; then
      echo "  $S BAM$i 파이프 실패:$bad"; echo "$S BAM$i$bad" >> log/pipe_failed.txt
      ok=0; break
    fi
    cat "tmp/${S}_${i}.hits.fastq" >> "tmp/${S}.hits.fastq"
    tail -n +2 "tmp/${S}_${i}.5s.tsv" >> "tmp/${S}.5s.body"
    rm -f "tmp/${S}_${i}.hits.fastq" "tmp/${S}_${i}.5s.tsv"
    fr=$(df -BG / | tail -1 | awk '{print $4}' | tr -d G)
    if [ "$fr" -lt 10 ]; then
      echo "  $S 디스크 $fr GB — 중단"; echo "$S" >> log/disk_abort.txt
      rm -f "tmp/${S}"*; ok=0; break
    fi
  done
  [ $ok -eq 1 ] || { echo "$S" >> log/failed.txt; continue; }

  nr=$(( $(wc -l < "tmp/${S}.hits.fastq") / 4 ))
  [ "$nr" -gt 0 ] || { echo "  $S 히트 0"; echo "$S" >> log/zero_hits.txt; continue; }
  { echo -e "read_id\tread_len\tn_hsp\tn_gene\tbest_pid\tspacings\thas_MM"; cat "tmp/${S}.5s.body"; } > "tmp/${S}.5s.tsv"

  # tandem map alignment, the reference is 90kb so no lock or slot is needed (CHM13 needed serialization because of its 7GB index)
  awk 'NR%4==1{sub(/^@/,">"); print $1} NR%4==2{print}' "tmp/${S}.hits.fastq" > "tmp/${S}.fa"
  minimap2 -x map-ont -t "$THR" -c --secondary=no tandem40.fa "tmp/${S}.fa" > "tmp/${S}.paf" 2>/dev/null
  [ -s "tmp/${S}.paf" ] || { echo "  $S PAF 빈파일"; echo "$S" >> log/empty_paf.txt; continue; }

  python3 summarize4.py --sample "$S" --tech ont --hits "tmp/${S}.5s.tsv" \
      --paf "tmp/${S}.paf" --fastq "tmp/${S}.hits.fastq" --nts 5S_unit_discseg.fa \
      --out "sum/${S}" --array "array/${S}.array.fastq.gz" --threads "$THR" \
      > "log/${S}.sum.log" 2>&1
  [ $? -eq 0 ] && [ -s "sum/${S}.summary.tsv" ] || { echo "  $S 요약 실패"; echo "$S" >> log/sum_failed.txt; continue; }

  # keep all hit reads (feasible because v5 keeps only 0.036%)
  pigz -p "$PIGZ_T" -c "tmp/${S}.hits.fastq" > "hits/${S}.hits.fastq.gz"
  hz=$(stat -c%s "hits/${S}.hits.fastq.gz" 2>/dev/null || echo 0)
  arr=$(stat -c%s "array/${S}.array.fastq.gz" 2>/dev/null || echo 0)
  [ "$hz" -gt 0 ] || { echo "  $S hits 보관 실패"; echo "$S" >> log/no_hits.txt; continue; }

  rm -f "tmp/${S}"*
  t1=$(date +%s)
  echo "$S ok v5 reads=$nr hits_gz=$hz array=$arr sec=$((t1-t0))" > "done/${S}.done"
  echo "  [$(date -u +%H:%M:%S)] $S 완료 read=$nr hits=$((hz/1048576))MB array=$((arr/1048576))MB $((t1-t0))s"
done < "$LIST"
echo "[$(date -u +%H:%M:%S)] V5_DONE 완료=$(ls done | wc -l)"
