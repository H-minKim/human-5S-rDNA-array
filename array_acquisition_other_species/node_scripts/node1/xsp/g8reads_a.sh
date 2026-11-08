set -u
X=~/xsp
for pair in "gorilla_mGorGor1 gorilla" "orang_sumatran_mPonAbe1 orang_sumatran"; do
  set -- $pair
  sp=$1; nm=$2
  mkdir -p $X/g8s/md_$sp $X/g8reads/$sp
  cp $X/g8s/g8map/${sp}__*.g8map.tsv $X/g8s/md_$sp/ 2>/dev/null
  python3 $X/units_matrix_g8.py --fa $X/asm/$sp.fna --mapdir $X/g8s/md_$sp --samples $sp --bam "$X/ape/out_$nm/{s}.5S.bam" --meth "$X/ape/out_$nm/{s}.5S.reads.tsv" --out $X/g8reads/$sp
done
mkdir -p $X/g8s/md_arabidopsis_ColCC $X/g8reads/arabidopsis_ColCC
cp $X/g8s/g8map/arabidopsis_ColCC__*.g8map.tsv $X/g8s/md_arabidopsis_ColCC/ 2>/dev/null
python3 $X/units_matrix_g8.py --fa $X/asm/arabidopsis_ColCC.fna --mapdir $X/g8s/md_arabidopsis_ColCC --samples ColCC_ATRun3 --bam "$X/ath/final/{s}.5S.bam" --meth "$X/ath/final/{s}.5S.reads.tsv" --out $X/g8reads/arabidopsis_ColCC
echo "node1 끝 $(date -u +%H:%M) UTC"
