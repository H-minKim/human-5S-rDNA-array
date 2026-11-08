#!/bin/bash
# Fiber-seq 배치 -> AOJV 표 (v11 aoplace 와 완전히 같은 조건)
#   지도  ~/CONSENSUS_MAP_G6.tsv · 요약 ~/CONSENSUS_SUM_G6.tsv
#   정본  ~/hap_seq_node.tsv + ~/cmp4d
#   → 배치는 G6 틀로만 정하고, 좌표틀 정합은 ao_g8cut_v12 가 한다.
#
#   $1 = FS | PB      $2 = 병렬 수      $3 = 태그 접두
#   FS : ONT Fiber-seq  K=12 Δ>=10 MINARR=26880 MAXFL=무제한 · read 이름 16 자
#   PB : PacBio        K=4  Δ>=5  MINARR=8960  MAXFL=무제한 · read 이름 전체
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
echo "$(hostname) [$KIND] 보유 $(wc -l < "$W/all.txt") · 할 것 $n · 병렬 $NP"
[ "$n" -eq 0 ] && exit 0

split -n "l/$NP" -d "$W/todo.txt" "$W/part_"
i=0
for p in "$W"/part_*; do
  [ -s "$p" ] || continue
  env W="$W" TH=6 TAG="${PFX}${i}" $GATE \
    MAPF=~/CONSENSUS_MAP_G6.tsv SUMF=~/CONSENSUS_SUM_G6.tsv \
    HSNF=~/hap_seq_node.tsv CMPD=~/cmp4d \
    nohup python3 "$CODE" "$p" > "$W/log_${PFX}${i}.txt" 2>&1 &
  echo "  띄움 ${PFX}${i} : $(wc -l < "$p") 개체"
  i=$((i + 1))
done
sleep 3
echo "  실행중 $(ps -eo args | grep -cE '^python3 .*(fsplace_v12|pbfsplace_v11)')"
