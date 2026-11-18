#!/bin/bash
# Fiber-seq placement -> the AOJV table, under exactly the same conditions as the v11 placement.
#   The map is CONSENSUS_MAP_G6.tsv with the summary CONSENSUS_SUM_G6.tsv.
#   The canonical sequences are hap_seq_node.tsv plus cmp4d.
#   The placement is decided in the G6 frame alone, and the coordinate-frame reconciliation is done by ao_g8cut_v12.
#
#   $1 = FS or PB, $2 = the number of parallel jobs, $3 = the tag prefix.
#   FS, ONT Fiber-seq: K = 12, margin >= 10, MINARR = 26880, MAXFL unlimited, read names truncated to 16 characters.
#   PB, PacBio: K = 4, margin >= 5, MINARR = 8960, MAXFL unlimited, full read names.
set -e
KIND="$1"; NP="${2:-4}"; PFX="${3:-n}"

if [ "$KIND" = "FS" ]; then
  W=~/fsv12_run; CODE=~/fsplace_v12.py; SRC=~/fs41/out
  GATE="K=12 DMIN=10 MINARR=26880 MAXFL=999999999"
  PAT='s/\.hit\.fq$//'; FILT='.'
else
  W=~/pbv12_run; CODE=~/pbfsplace_v11.py; SRC=~/pbfs/out
  GATE="K=4 DMIN=5 MINARR=8960 MAXFL=999999999"
  PAT='s/\.hit\.fq$//'; FILT='(__PBFS|__MO|__FIRE)$'
fi
mkdir -p "$W"

cd "$SRC"
ls -1 *.hit.fq 2>/dev/null | sed "$PAT" | grep -E "$FILT" | sort -u > "$W/all.txt"
if [ -s "$W/done.txt" ]; then
  sort -u "$W/done.txt" > "$W/.d"; comm -23 "$W/all.txt" "$W/.d" > "$W/todo.txt"
else cp "$W/all.txt" "$W/todo.txt"; fi

n=$(wc -l < "$W/todo.txt")
echo "$(hostname) [$KIND] available $(wc -l < "$W/all.txt"), to do $n, parallel $NP"
[ "$n" -eq 0 ] && exit 0

split -n "l/$NP" -d "$W/todo.txt" "$W/part_"
i=0
for p in "$W"/part_*; do
  [ -s "$p" ] || continue
  env W="$W" TH=6 TAG="${PFX}${i}" $GATE \
    MAPF=~/CONSENSUS_MAP_G6.tsv SUMF=~/CONSENSUS_SUM_G6.tsv \
    HSNF=~/hap_seq_node.tsv CMPD=~/cmp4d \
    nohup python3 "$CODE" "$p" > "$W/log_${PFX}${i}.txt" 2>&1 &
  echo "  launched ${PFX}${i} : $(wc -l < "$p") individuals"
  i=$((i + 1))
done
sleep 3
echo "  running $(ps -eo args | grep -cE '^python3 .*(fsplace_v12|pbfsplace_v11)')"
