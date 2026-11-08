set -u
X=~/xsp
for pair in "bonobo_mPanPan1 bonobo" "orang_bornean_mPonPyg2 orang_bornean"; do
  set -- $pair
  sp=$1; nm=$2
  mkdir -p $X/g8s/md_$sp $X/g8reads/$sp
  cp $X/g8s/g8map/${sp}__*.g8map.tsv $X/g8s/md_$sp/ 2>/dev/null
  python3 $X/units_matrix_g8.py --fa $X/asm/$sp.fna --mapdir $X/g8s/md_$sp --samples $sp --bam "$X/ape/out_$nm/{s}.5S.bam" --meth "$X/ape/out_$nm/{s}.5S.reads.tsv" --out $X/g8reads/$sp
done
echo "node2 끝 $(date -u +%H:%M) UTC"
