#!/bin/bash
# 유인원 마무리 지휘 (로컬 백그라운드)
#   ① 세 노드의 선별(각 2 종)이 끝나기를 기다린다
#   ② 그 노드가 맡은 종의 조립본·cmsearch 표·구조표·modkit 을 노드에 갖춰 준다 (선별 결과는 그 노드에 있으니 옮기지 않는다)
#   ③ 노드에서 ape_align.sh 두 종을 차례로 : fastq(MM/ML) -> 정렬 -> array read BAM -> modkit -> read x 유닛 행렬
set -u
K=$NODE_KEY
S="-n -i $K -o StrictHostKeyChecking=no -o ConnectTimeout=25"
R="-i $K -o StrictHostKeyChecking=no -o ConnectTimeout=25"
O="-q -i $K -o StrictHostKeyChecking=no"
W=$NODE_HOST_1
B=$NODE_HOST_3
A=$NODE_HOST_2
P="$SCRATCH_DIR/xsp"
for ip in $W $B $A; do scp $O "$P/ape_align.sh" "$P/ape_prep.sh" ubuntu@$ip:~/xsp/; done
echo "스크립트 배포 $(date -u +%H:%M) UTC"
for i in $(seq 1 400); do
  cw=$(ssh $S ubuntu@$W 'grep -c "선별 끝" ~/xsp/ape/run.log 2>/dev/null')
  cb=$(ssh $S ubuntu@$B 'grep -c "선별 끝" ~/xsp/ape/run.log 2>/dev/null')
  ca=$(ssh $S ubuntu@$A 'grep -c "선별 끝" ~/xsp/ape/run.log 2>/dev/null')
  echo "$(date -u +%H:%M) 선별 : node3 ${cw:-0}/2 · node2 ${cb:-0}/2 · node1 ${ca:-0}/2"
  [ "${cw:-0}" -ge 2 ] && [ "${cb:-0}" -ge 2 ] && [ "${ca:-0}" -ge 2 ] && break
  sleep 300
done
echo "== 선별 끝 $(date -u +%H:%M) UTC · flowcell 별 (이름 · 걸린 read · MM 있는 read · 종료코드 · 초)"
for ip in $W $B $A; do ssh $S ubuntu@$ip 'cat ~/xsp/ape/*/*.done 2>/dev/null | sort'; done
# 작은 파일 옮기기 (구조표 · cmsearch 표 · modkit)
for ip in $A $B; do
  ssh $S ubuntu@$W "tar cf - -C ~/xsp/struct ALL.struct.tsv" | ssh $R ubuntu@$ip 'mkdir -p ~/xsp/struct && tar xf - -C ~/xsp/struct'
  if [ "$(ssh $S ubuntu@$ip 'ls ~/xsp/tools/modkit 2>/dev/null | wc -l')" -lt 1 ]; then
    ssh $S ubuntu@$W "tar cf - -C ~/xsp/tools modkit" | ssh $R ubuntu@$ip 'mkdir -p ~/xsp/tools && tar xf - -C ~/xsp/tools'
  fi
done
ssh $S ubuntu@$W "tar cf - -C ~/xsp/cms gorilla_mGorGor1.tbl orang_sumatran_mPonAbe1.tbl" | ssh $R ubuntu@$A 'mkdir -p ~/xsp/cms && tar xf - -C ~/xsp/cms'
ssh $S ubuntu@$W "tar cf - -C ~/xsp/cms bonobo_mPanPan1.tbl orang_bornean_mPonPyg2.tbl" | ssh $R ubuntu@$B 'mkdir -p ~/xsp/cms && tar xf - -C ~/xsp/cms'
echo "표 · 도구 옮김 $(date -u +%H:%M) UTC"
ssh $S ubuntu@$A "bash ~/xsp/ape_prep.sh gorilla_mGorGor1 GCA_029281585.3; bash ~/xsp/ape_prep.sh orang_sumatran_mPonAbe1 GCA_028885655.3"
ssh $S ubuntu@$B "bash ~/xsp/ape_prep.sh bonobo_mPanPan1 GCA_029289425.3; bash ~/xsp/ape_prep.sh orang_bornean_mPonPyg2 GCA_028885625.3"
ssh $S ubuntu@$A "setsid nohup bash -c 'bash ~/xsp/ape_align.sh gorilla_mGorGor1 auto ~/xsp/ape/gorilla ~/xsp/ape/out_gorilla; bash ~/xsp/ape_align.sh orang_sumatran_mPonAbe1 auto ~/xsp/ape/orang_sumatran ~/xsp/ape/out_orang_sumatran' > ~/xsp/ape/align.log 2>&1 < /dev/null & sleep 2; echo node1 정렬 기동"
ssh $S ubuntu@$B "setsid nohup bash -c 'bash ~/xsp/ape_align.sh bonobo_mPanPan1 auto ~/xsp/ape/bonobo ~/xsp/ape/out_bonobo; bash ~/xsp/ape_align.sh orang_bornean_mPonPyg2 auto ~/xsp/ape/orang_bornean ~/xsp/ape/out_orang_bornean' > ~/xsp/ape/align.log 2>&1 < /dev/null & sleep 2; echo node2 정렬 기동"
ssh $S ubuntu@$W "setsid nohup bash -c 'bash ~/xsp/ape_align.sh chimp_mPanTro3 auto ~/xsp/ape/chimp ~/xsp/ape/out_chimp; bash ~/xsp/ape_align.sh siamang_mSymSyn1 auto ~/xsp/ape/siamang ~/xsp/ape/out_siamang' > ~/xsp/ape/align.log 2>&1 < /dev/null & sleep 2; echo node3 정렬 기동"
for i in $(seq 1 200); do
  n=0
  for ip in $W $B $A; do n=$((n + $(ssh $S ubuntu@$ip 'grep -c "^끝 " ~/xsp/ape/align.log 2>/dev/null'))); done
  echo "$(date -u +%H:%M) 정렬 끝난 종 $n/6"
  [ "$n" -ge 6 ] && break
  sleep 300
done
echo "== 유인원 전 과정 끝 $(date -u +%H:%M) UTC"
for ip in $W $B $A; do ssh $S ubuntu@$ip 'cat ~/xsp/ape/align.log'; done
