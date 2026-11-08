#!/bin/bash
# 애기장대 : array read 만 남긴 pod5 조각을 노드별로 나눈다 (크기 기준 · 코어 가중 A 16 · W 28 · B 16)
#   남긴 read 가 1 개 이상인 조각만 (members.tsv 5 열 > 0)
#   산출 : grp_A.txt · grp_W.txt · grp_B.txt (절대경로 목록 — 세 노드 모두 /home/ubuntu 라 경로가 같다)
set -u
H=~/xsp/ath
rm -f $H/grp_A.txt $H/grp_W.txt $H/grp_B.txt
awk -F'\t' '$5>0 {print $1}' $H/members.tsv | while read b; do
  f=$H/sub/$b
  [ -s "$f" ] && printf '%s\t%s\n' "$(stat -c %s "$f")" "$f"
done | sort -k1,1nr | awk -F'\t' -v H="$H" 'BEGIN{w["A"]=16; w["W"]=28; w["B"]=16}
  {best=""; for(k in w){r=(ld[k]+$1)/w[k]; if(best=="" || r<br){best=k; br=r}}; ld[best]+=$1; n[best]++; print $2 > (H "/grp_" best ".txt")}
  END{for(k in w) printf "%s : 조각 %d · %.2f GB\n", k, n[k], ld[k]/1e9}'
