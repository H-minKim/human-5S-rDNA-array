#!/bin/bash
# 애기장대 재콜링 배분기 (로컬에서 백그라운드) — 노는 코어를 모두 쓴다
#   흘려받기가 끝나기를 기다리지 않고, 필터가 끝난 pod5 조각을 ~300 read 묶음으로 노드에 맡긴다 (ath_claim.sh · flock)
# Cores : node1 2-15 (2 cores left for the streaming share)
# node3 24-31 while the zebra finch re-calling runs, 0-31 afterwards
# node2 0-15 after the zebra finch re-calling finishes
# the pod5 share is streamed straight from node1 to the other nodes (tar | ssh), and the finished BAM is sent back to node1
#   흘려받기가 끝나고 남은 조각이 없고 도는 묶음도 없으면 ath_finalize.sh
set -u
K=$NODE_KEY
S="-n -i $K -o StrictHostKeyChecking=no -o ConnectTimeout=20"
R="-i $K -o StrictHostKeyChecking=no -o ConnectTimeout=20"
O="-q -i $K -o StrictHostKeyChecking=no"
declare -A HOST=([A]=$NODE_HOST_2 [W]=$NODE_HOST_1 [B]=$NODE_HOST_3)
declare -A RUN=([A]="" [W]="" [B]="")
A=${HOST[A]}
P="$SCRATCH_DIR/xsp"
scp $O "$P/ath_claim.sh" "$P/ath_finalize.sh" "$P/ath_dorado.sh" ubuntu@$A:~/xsp/ath/ && echo "scripts uploaded $(date -u +%H:%M) UTC"
for n in W B; do scp $O "$P/ath_dorado.sh" ubuntu@${HOST[$n]}:~/xsp/ath/; done
fdone() {
  case $1 in
    W) ssh $S ubuntu@${HOST[W]} 'grep -c "재콜링 끝" ~/xsp/zfont/dorado_w3.log 2>/dev/null' ;;
    B) ssh $S ubuntu@${HOST[B]} 'grep -c "재콜링 끝" ~/xsp/zfont/dorado_b.log 2>/dev/null' ;;
  esac
}
cores() {
  case $1 in
    # other streaming jobs use about 6-8 cores per node, so that many are left free
    A) echo 2-13 ;;
    W) [ "$(fdone W)" -ge 1 ] 2> /dev/null && echo 0-23 || echo 24-31 ;;
    B) [ "$(fdone B)" -ge 1 ] 2> /dev/null && echo 0-11 || echo NONE ;;
  esac
}
while true; do
  sd=$(ssh $S ubuntu@$A 'grep -c "전체 끝" ~/xsp/ath/stream_full.log 2>/dev/null')
  sd=${sd:-0}
  none=0
  for n in A W B; do
    h=${HOST[$n]}
    id=${RUN[$n]}
    if [ -n "$id" ]; then
      if [ "$(ssh $S ubuntu@$h "ls ~/xsp/ath/bc/$id.done 2> /dev/null | wc -l")" -ge 1 ]; then
        if [ "$n" != A ]; then
          ssh $S ubuntu@$h "tar cf - -C ~/xsp/ath/bc $id.bam $id.bam.bai $id.done" | ssh $R ubuntu@$A 'tar xf - -C ~/xsp/ath/bc'
        fi
        echo "$(date -u +%H:%M) $n $id 끝 : $(ssh $S ubuntu@$A "cat ~/xsp/ath/bc/$id.done")"
        RUN[$n]=""
      else
        continue
      fi
    fi
    c=$(cores $n)
    [ "$c" = NONE ] && continue
    id=$(ssh $S ubuntu@$A "bash ~/xsp/ath/ath_claim.sh 300 $([ "$sd" -ge 1 ] && echo 1 || echo 0)")
    if [ "$id" = NONE ] || [ -z "$id" ]; then
      none=$((none + 1))
      continue
    fi
    if [ "$n" != A ]; then
      ssh $S ubuntu@$A "tar cf - -T ~/xsp/ath/batch/$id.txt ~/xsp/ath/batch/$id.txt 2> /dev/null" | ssh $R ubuntu@$h 'mkdir -p ~/xsp/ath/sub ~/xsp/ath/batch ~/xsp/ath/bc && tar xf - -C /'
    fi
    ssh $S ubuntu@$h "setsid nohup bash ~/xsp/ath/ath_dorado.sh ~/xsp/ath/batch/$id.txt $c $id > ~/xsp/ath/dorado_$id.log 2>&1 < /dev/null &"
    echo "$(date -u +%H:%M) $n $id 기동 · 코어 $c · 조각 $(ssh $S ubuntu@$A "wc -l < ~/xsp/ath/batch/$id.txt") · 흘려받기끝=$sd"
    RUN[$n]=$id
  done
  if [ "$sd" -ge 1 ] && [ -z "${RUN[A]}${RUN[W]}${RUN[B]}" ]; then
    left=$(ssh $S ubuntu@$A 'bash ~/xsp/ath/ath_claim.sh 1 1')
    if [ "$left" = NONE ]; then
      break
    fi
    ssh $S ubuntu@$A "setsid nohup bash ~/xsp/ath/ath_dorado.sh ~/xsp/ath/batch/$left.txt 2-15 $left > ~/xsp/ath/dorado_$left.log 2>&1 < /dev/null &"
    RUN[A]=$left
    echo "$(date -u +%H:%M) A $left 기동 (나머지)"
  fi
  sleep 60
done
echo "== 모든 묶음 끝 $(date -u +%H:%M) UTC"
ssh $S ubuntu@$A 'echo "맡긴 조각 $(wc -l < ~/xsp/ath/batch/claimed.txt) · 남김 read 있는 멤버 $(awk -F"\t" "NF>=5 && \$5>0" ~/xsp/ath/members.tsv | wc -l)"; bash ~/xsp/ath/ath_finalize.sh > ~/xsp/ath/final.log 2>&1; cat ~/xsp/ath/final.log'
