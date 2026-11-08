#!/bin/bash
# 금화조 : 재콜링 read 가 원 fastq 정렬과 같은 array 에 붙었나 (유닛 번호보다 먼저 array 배정부터)
#   array 계열 = chr2_A · chr2_B · chr4 · chr9 (mat·pat 는 같은 계열로 본다)
#   입력 : final/orig_aln.tsv (원 PAF 전 줄, pat chr2 는 mat 좌표로 옮김) · final/new_aln.tsv (재콜링 primary)
#   분류 : 같은 array · 다른 array · 새 정렬이 array 밖 · 원 정렬이 array 밖
set -u
Z=~/xsp/zfont
F=$Z/final
awk -F'\t' 'FILENAME==ARGV[1]{n++; c[n]=$1; s[n]=$2; e[n]=$3; nm=$4; sub(/^(mat|pat)_/, "", nm); if(nm !~ /^chr2/) sub(/_[0-9]+$/, "", nm); f[n]=nm; next}
  FILENAME==ARGV[2]{for(i=1;i<=n;i++) if($2==c[i] && $3<e[i] && $4>s[i]) o[$1]=o[$1] "|" f[i] "|"; next}
  {nf=""; for(i=1;i<=n;i++) if($2==c[i] && $3<e[i] && $4>s[i]) nf=f[i];
   if(nf=="") cls="새정렬_array밖"; else if(!($1 in o)) cls="원정렬_array밖"; else if(index(o[$1], "|" nf "|")>0) cls="같은_array"; else cls="다른_array";
   k[cls]++; kk[nf "\t" cls]++}
  END{for(x in k) print k[x] "\t" x; print "--- 새 정렬 array 계열별"; for(x in kk) print kk[x] "\t" x}' \
  $Z/zfinch_v03_5S_arrays.bed $F/orig_aln.tsv $F/new_aln.tsv | sort -k2,2 -k1,1nr
# 같은 array 안에서 시작점 차이 (유닛 368 bp 배수 분포)
echo "--- 같은 array 안 시작점 차이 (poscmp best, 유닛=368bp)"
awk -F'\t' '$4>=0 {u=int($4/368); if(u>5)u=">5"; h[u]++} END{for(x in h) print h[x] "\t" x " 유닛"}' $F/poscmp.tsv | sort -k2,2n
