set -u
cd ~/xsp
rm -rf smap2; mkdir -p smap2
for sp in human_CHM13v2 chimp_mPanTro3 bonobo_mPanPan1 gorilla_mGorGor1 orang_sumatran_mPonAbe1 orang_bornean_mPonPyg2 siamang_mSymSyn1 rhesus_T2T_MMU8v2 mouse_C57BL6J_T2T zebrafish_GRCz12tu arabidopsis_ColCC zfinch_v03_dip; do
  python3 ~/xsp/species_map2.py $sp ~/xsp/asm/$sp.fna ~/xsp/cms/$sp.tbl ~/xsp/struct/ALL.struct.tsv ~/xsp/smap2 20000
done
head -1 smap2/MAP_SUMMARY_human_CHM13v2.tsv > smap2/MAP_SUMMARY.tsv
for f in smap2/MAP_SUMMARY_*.tsv; do case $f in *"/MAP_SUMMARY.tsv") continue;; esac; tail -n +2 $f >> smap2/MAP_SUMMARY.tsv; done
echo "2판 끝 $(date -u +%H:%M) UTC · array $(($(wc -l < smap2/MAP_SUMMARY.tsv)-1))"
awk -F"\t" "NR>1{t+=\$8; m+=\$9; b+=\$10; n+=\$11} END{printf \"유닛 합 %d : minimap2 %d · blastn %d · 실패 %d (%.2f%%)\n\", t, m, b, n, 100*n/t}" smap2/MAP_SUMMARY.tsv
