set -u
cd ~/xsp
rm -rf g8s; mkdir -p g8s
for sp in human_CHM13v2 chimp_mPanTro3 bonobo_mPanPan1 gorilla_mGorGor1 orang_sumatran_mPonAbe1 orang_bornean_mPonPyg2 siamang_mSymSyn1 rhesus_T2T_MMU8v2 mouse_C57BL6J_T2T zebrafish_GRCz12tu arabidopsis_ColCC zfinch_v03_dip; do
  python3 ~/xsp/species_g8map.py $sp ~/xsp/asm/$sp.fna ~/xsp/cms/$sp.tbl ~/xsp/ref2/REFUNIT_PARAMS.tsv ~/xsp/ref2/refunit ~/xsp/g8s
done
head -1 g8s/G8MAP_SUM_human_CHM13v2.tsv > g8s/G8MAP_SUM.tsv
for f in g8s/G8MAP_SUM_*.tsv; do case $f in *"/G8MAP_SUM.tsv") continue;; esac; tail -n +2 $f >> g8s/G8MAP_SUM.tsv; done
echo "끝 $(date -u +%H:%M) UTC · array $(($(wc -l < g8s/G8MAP_SUM.tsv)-1))"
awk -F"\t" "NR>1{u+=\$7; b+=\$8; o+=\$9; nf+=\$10; na+=\$11} END{printf \"유닛 %d · 끊김 %d · OK %d (%.2f%%) · NO_FULL %d · NOALN %d\n\", u, b, o, 100*o/u, nf, na}" g8s/G8MAP_SUM.tsv
