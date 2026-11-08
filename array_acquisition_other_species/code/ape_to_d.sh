#!/bin/bash
# 유인원 산출물을 통일 규격으로 D 에 가져온다 (03_reads · 04_bam · 05_bedmethyl · 06_aux)
#   표에는 species_key · sample 열을 앞에 붙인다 (sample = 종키와 같음 : 종마다 개체 1 명)
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
K=$NODE_KEY
O="-q -i $K -o StrictHostKeyChecking=no"
V="${RDNA_OTHERSP_DIR:?set RDNA_OTHERSP_DIR}/data_v2"
T="$SCRATCH_DIR/relay/ape"
mkdir -p "$T"
get() {
  ip=$1; nm=$2; sp=$3
  d="$T/$sp"
  mkdir -p "$d"
  scp $O "ubuntu@$ip:xsp/ape/out_$nm/units/*" "$d/" 2> /dev/null
  scp $O "ubuntu@$ip:xsp/ape/out_$nm/$sp.5S.bam" "ubuntu@$ip:xsp/ape/out_$nm/$sp.5S.bam.bai" "$d/" 2> /dev/null
  scp $O "ubuntu@$ip:xsp/ape/out_$nm/$sp.5S.cpg.bed" "$d/" 2> /dev/null
  scp $O "ubuntu@$ip:xsp/ape/out_$nm/use.txt" "$d/flowcells_used.txt" 2> /dev/null
  mkdir -p "$V/03_reads/$sp" "$V/04_bam/$sp" "$V/05_bedmethyl/$sp" "$V/06_aux/$sp"
  for f in "$d"/*.read_unit.tsv "$d"/5S_unit_summary.tsv "$d"/5S_array_summary.tsv; do
    [ -s "$f" ] || continue
    b=$(basename "$f")
    case "$b" in
      *.read_unit.tsv) out="$V/03_reads/$sp/${sp}__${sp}.read_unit.tsv" ;;
      *) out="$V/03_reads/$sp/$b" ;;
    esac
    awk -F'\t' -v OFS='\t' -v sp="$sp" 'NR==1{print "species_key", "sample", $0; next} {print sp, sp, $0}' "$f" > "$out"
    echo "  $(basename "$out")  $(( $(wc -l < "$out") - 1 )) 행"
  done
  [ -s "$d/$sp.5S.bam" ] && cp "$d/$sp.5S.bam" "$d/$sp.5S.bam.bai" "$V/04_bam/$sp/" && echo "  BAM $(du -h "$d/$sp.5S.bam" | cut -f1)"
  if [ -s "$d/$sp.5S.cpg.bed" ]; then gzip -f -1 -k "$d/$sp.5S.cpg.bed"; cp "$d/$sp.5S.cpg.bed.gz" "$V/05_bedmethyl/$sp/"; fi
  [ -s "$d/flowcells_used.txt" ] && cp "$d/flowcells_used.txt" "$V/06_aux/$sp/"
}
echo "== 유인원 가져오기 $(date -u +%H:%M) UTC"
get $NODE_HOST_1 chimp chimp_mPanTro3
get $NODE_HOST_1 siamang siamang_mSymSyn1
get $NODE_HOST_3 bonobo bonobo_mPanPan1
get $NODE_HOST_3 orang_bornean orang_bornean_mPonPyg2
get $NODE_HOST_2 gorilla gorilla_mGorGor1
get $NODE_HOST_2 orang_sumatran orang_sumatran_mPonAbe1
echo "== 끝 $(date -u +%H:%M) UTC"
