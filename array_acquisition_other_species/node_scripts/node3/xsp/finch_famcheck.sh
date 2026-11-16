#!/bin/bash
# Zebra finch: did the re-basecalled read land on the same array as the original fastq alignment? The array assignment is checked before the unit index.
#   Array families: chr2_A, chr2_B, chr4, chr9; the maternal and paternal copies count as the same family.
#   Inputs: final/orig_aln.tsv (all original PAF rows, with paternal chr2 shifted into maternal coordinates) and final/new_aln.tsv (the re-basecalled primary alignments).
#   Categories: same array, different array, the new alignment lies outside any array, the original alignment lies outside any array.
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
# The difference in start point within the same array, as a distribution over multiples of the 368 bp unit.
echo "--- 같은 array 안 시작점 차이 (poscmp best, 유닛=368bp)"
awk -F'\t' '$4>=0 {u=int($4/368); if(u>5)u=">5"; h[u]++} END{for(x in h) print h[x] "\t" x " 유닛"}' $F/poscmp.tsv | sort -k2,2n
