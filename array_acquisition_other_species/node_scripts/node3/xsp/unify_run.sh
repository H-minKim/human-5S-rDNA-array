set -u
X=~/xsp
cd $X
for sp in human_CHM13v2 chimp_mPanTro3 bonobo_mPanPan1 gorilla_mGorGor1 orang_sumatran_mPonAbe1 orang_bornean_mPonPyg2 siamang_mSymSyn1 rhesus_T2T_MMU8v2 mouse_C57BL6J_T2T zebrafish_GRCz12tu arabidopsis_ColCC zfinch_v03_dip; do
  python3 $X/species_map.py $sp $X/asm/$sp.fna $X/cms/$sp.tbl $X/struct/ALL.struct.tsv $X/smap 20000
done
head -1 $X/smap/MAP_SUMMARY_human_CHM13v2.tsv > $X/smap/MAP_SUMMARY.tsv
for f in $X/smap/MAP_SUMMARY_*.tsv; do case $f in *MAP_SUMMARY.tsv) continue;; esac; tail -n +2 $f >> $X/smap/MAP_SUMMARY.tsv; done
echo "정본·지도 끝 $(date -u +%H:%M) UTC · array $(($(wc -l < $X/smap/MAP_SUMMARY.tsv)-1))"
for s in Control1 Control2 Control3 Control4 Control5 Control6; do
  x=$X/zf/meth/$s.5S.reads.tsv
  c=$(head -1 $x | tr "\t" "\n" | grep -n -x ref_position | cut -d: -f1)
  { head -1 $x; tail -n +2 $x | LC_ALL=C sort -t "	" -k1,1 -k${c},${c}n -S 2G --parallel 4 -T $X/zf; } > $x.s && mv $x.s $x
done
python3 $X/units_matrix.py --fa $X/asm/zebrafish_GRCz12tu.fna --bed $X/zf/grcz12_5S_arrays.bed --tbl $X/cms/zebrafish_GRCz12tu.tbl --samples Control1,Control2,Control3,Control4,Control5,Control6 --bam "$X/zf/bam/{s}.g12.5S.bam" --meth "$X/zf/meth/{s}.5S.reads.tsv" --out $X/zf/units2
echo "제브라피시 통일 끝 $(date -u +%H:%M) UTC"
