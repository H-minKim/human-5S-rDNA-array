#!/bin/bash
# v12 reassembly, re-assemble the 112 samples of our hifiasm slots (188) while keeping the flanks (2026-09-15)
#
#   Existing outputs are never touched. Writes only to the new folder ~/reasm_v12.
#   Same command as the original : hifiasm --ont -o <sample> -t <TH>   (a different mode would make the comparison invalid)
# Disk : when a sample finishes, intermediate files such as .bin are deleted immediately
#             kept = p_ctg-type gfa converted to fasta + log + lowQ.bed
#   Arguments : list file (sample<TAB>read path), number of parallel jobs, threads
set -u
LIST=${1:?목록}
NP=${2:-2}
TH=${3:-12}
W=~/reasm_v12
mkdir -p "$W"/out "$W"/log "$W"/tmp
[ -x ~/hifiasm ] || { echo "$(hostname) hifiasm 없음"; exit 1; }

one() {
  S=$1; Q=$2
  T="$W/tmp/$S"
  [ -s "$W/out/$S.p_ctg.fa" ] && { echo "  $S 이미 있음"; return; }
  rm -rf "$T"; mkdir -p "$T"
  t0=$(date +%s)
  ( cd "$T" && timeout 21600 ~/hifiasm --ont -o "$S" -t "$TH" "$Q" ) > "$W/log/$S.log" 2>&1
  rc=$?
  t1=$(date +%s)
  # gfa -> fasta (only records with sequence)
  n=0
  for g in "$T"/*.p_ctg.gfa; do
    [ -s "$g" ] || continue
    b=$(basename "$g" .gfa)
    awk '/^S/{print ">"$2"\n"$3}' "$g" > "$W/out/$b.fa"
    [ -s "$W/out/$b.fa" ] && n=$((n+1)) || rm -f "$W/out/$b.fa"
  done
  # one more copy under the representative name (bp.p_ctg -> p_ctg)
  [ -s "$W/out/$S.bp.p_ctg.fa" ] && cp "$W/out/$S.bp.p_ctg.fa" "$W/out/$S.p_ctg.fa"
  cp "$T"/*.lowQ.bed "$W/out/" 2>/dev/null
  rm -rf "$T"                       # delete intermediate files immediately
  echo "  $S rc=$rc · $((t1-t0))초 · fasta $n · $(df -h / | tail -1 | awk '{print $4}') 남음"
}

i=0
while IFS=$'\t' read -r S Q; do
  [ -n "$S" ] || continue
  [ -s "$Q" ] || { echo "  $S read 없음 $Q"; continue; }
  one "$S" "$Q" &
  i=$((i+1))
  if [ $((i % NP)) -eq 0 ]; then wait; fi
done < "$LIST"
wait
echo "REASMDONE $(hostname) · fasta $(ls $W/out/*.p_ctg.fa 2>/dev/null | wc -l)"
