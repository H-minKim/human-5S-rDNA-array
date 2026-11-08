#!/bin/bash
# Regenerate the AWS own references. Input is CHM13v2.0 only. No SNU/local outputs are used at all.
#  canonical unit = revcomp(chr1:227,887,698-227,889,937)  <- confirmed with pident 100.00/2240bp
#  unit frame: AluY 1-318 / NTS1 319-1178 / GA-rich+(CA)n 1179-1551 / NTS2 1552-2121 / 5S gene 2122-2240
#  discriminating segment = NTS1 + NTS2 = 1430 bp
#  array(CHM13) chr1:227,744,873-228,024,958, TUC junction proximal 227,744,885 / distal 228,025,899
set -eu
W=~/s5s_tuc/ref; mkdir -p $W
CH=~/ectopic_sp/chm13v2.0.fa
U1=227887698; U2=227889937
PJ=227744885; DJ=228025899

samtools faidx -i $CH chr1:${U1}-${U2} | tail -n +2 | tr -d '\n' > $W/_unit.txt
python3 - <<PY
import io
W="$W"
u=io.open(W+"/_unit.txt").read().strip().upper()
assert len(u)==2240, "정본 유닛 길이 %d"%len(u)
io.open(W+"/unit_canon.fa","w",newline="\n").write(u">UNIT_CANON_2240 revcomp_CHM13_chr1_${U1}_${U2}\n%s\n"%u)
# 40 tandem repeats (TUC mapping reference)
io.open(W+"/tandem40.fa","w",newline="\n").write(u">tandem\n%s\n"%(u*40))
# discriminating segment = NTS1(319-1178) + NTS2(1552-2121), 1-based inclusive
d1,d2=u[318:1178],u[1551:2121]
io.open(W+"/5S_unit_discseg.fa","w",newline="\n").write(
    u">disc_unitpos_319_1178_len%d\n%s\n>disc_unitpos_1552_2121_len%d\n%s\n"%(len(d1),d1,len(d2),d2))
# 5S gene (2122-2240)
io.open(W+"/5S_gene_query.fa","w",newline="\n").write(u">5S_gene_unitpos_2122_2240\n%s\n"%u[2121:2240])
print("  unit=%d tandem=%d disc=%d+%d gene=%d"%(len(u),len(u)*40,len(d1),len(d2),len(u)-2121))
PY

# flank reference (50 kb outside the junction)
samtools faidx $CH chr1:$((PJ-50000))-$((PJ-1)) | sed "1s/.*/>FLANK_proximal_50000 CHM13_chr1_$((PJ-50000))_$((PJ-1))_junction_at_END/" > $W/FLANK_proximal_50000.fa
samtools faidx $CH chr1:$((DJ+1))-$((DJ+50000))  | sed "1s/.*/>FLANK_distal_50000 CHM13_chr1_$((DJ+1))_$((DJ+50000))_junction_at_START/" > $W/FLANK_distal_50000.fa

makeblastdb -in $W/5S_gene_query.fa -dbtype nucl -out $W/gene_db  >/dev/null
makeblastdb -in $W/tandem40.fa      -dbtype nucl -out $W/tandem_db >/dev/null
makeblastdb -in $W/unit_canon.fa    -dbtype nucl -out $W/unit_db   >/dev/null
rm -f $W/_unit.txt

echo "=== 재생성 결과 ==="
for f in unit_canon.fa tandem40.fa 5S_unit_discseg.fa 5S_gene_query.fa FLANK_proximal_50000.fa FLANK_distal_50000.fa; do
  printf "  %-28s %8d B  md5=%s\n" "$f" "$(stat -c%s $W/$f)" "$(md5sum $W/$f | cut -c1-12)"
done
echo "=== 기존 AWS ONT 캠페인본과 대조 (서열만) ==="
for p in tandem40.fa 5S_unit_discseg.fa 5S_gene_query.fa; do
  a=$(grep -v '^>' ~/s5s_ont/$p 2>/dev/null | tr -d '\n' | tr a-z A-Z | md5sum | cut -c1-12)
  b=$(grep -v '^>' $W/$p          | tr -d '\n' | tr a-z A-Z | md5sum | cut -c1-12)
  [ "$a" = "$b" ] && echo "  $p  일치 ($b)" || echo "  $p 불일치 ONT=$a 재생성=$b"
done
