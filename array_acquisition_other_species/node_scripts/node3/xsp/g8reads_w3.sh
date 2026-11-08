set -u
X=~/xsp
for sp in siamang_mSymSyn1 zebrafish_GRCz12tu zfinch_v03_dip; do
  mkdir -p $X/g8s/md_$sp $X/g8reads/$sp
  cp $X/g8s/g8map/${sp}__*.g8map.tsv $X/g8s/md_$sp/ 2>/dev/null
done
python3 $X/units_matrix_g8.py --fa $X/asm/siamang_mSymSyn1.fna --mapdir $X/g8s/md_siamang_mSymSyn1 --samples siamang_mSymSyn1 --bam "$X/ape/out_siamang/{s}.5S.bam" --meth "$X/ape/out_siamang/{s}.5S.reads.tsv" --out $X/g8reads/siamang_mSymSyn1
python3 $X/units_matrix_g8.py --fa $X/asm/zebrafish_GRCz12tu.fna --mapdir $X/g8s/md_zebrafish_GRCz12tu --samples Control1,Control2,Control3,Control4,Control5,Control6 --bam "$X/zf/bam/{s}.g12.5S.bam" --meth "$X/zf/meth/{s}.5S.reads.tsv" --out $X/g8reads/zebrafish_GRCz12tu
python3 $X/units_matrix_g8.py --fa $X/asm/zfinch_v03_dip.fna --mapdir $X/g8s/md_zfinch_v03_dip --samples bTaeGut7_LIG,bTaeGut7_UL --bam "$X/zfont/final/{s}.5S.bam" --meth "$X/zfont/final/{s}.5S.reads.tsv" --out $X/g8reads/zfinch_v03_dip
echo "node3 끝 $(date -u +%H:%M) UTC"
