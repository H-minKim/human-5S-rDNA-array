#!/bin/bash
# 금화조 : 제출본(v0.3) PacBio 메틸 bigwig 와 좌표가 같은 GenBank 판을 받아 5S array 를 찾는다
#   GCA_048771995.1 (bTaeGut7.mat) · GCA_048772025.1 (bTaeGut7.pat) — 2025-03 판
#   ① 받기 · 염색체 길이를 bigwig 와 대조 (이름표 = 제출본 seq_id_mapping / alias)
#   ② cmsearch RF00001 --cut_ga · array 목록
set -u
X=$HOME/xsp
Z=$X/zfinch
cd $Z || exit 1
getncbi() {
  a=$1; nm=$2
  [ -s $X/asm/$nm.fna.fai ] && return
  curl -sL --retry 3 -o $X/asm/$nm.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$a/download?include_annotation_type=GENOME_FASTA"
  unzip -o -q $X/asm/$nm.zip -d $X/asm/$nm.d && cat $(find $X/asm/$nm.d -name "*.fna") > $X/asm/$nm.fna && rm -rf $X/asm/$nm.d
  samtools faidx $X/asm/$nm.fna
}
getncbi GCA_048771995.1 zfinch_v03_mat
getncbi GCA_048772025.1 zfinch_v03_pat
[ -s sub_pat.bw ] || curl -sL -o sub_pat.bw "https://genomeark.s3.amazonaws.com/species/Taeniopygia_guttata/bTaeGut7/submission/v0.3/annotations/pat/bTaeGut7.pat%2BWZ.cur.20250313.PBmethylation.v0.1.bw"
$X/tools/bigWigInfo -chroms sub_pat.bw 2>&1 | awk 'NF==3{print $1"\t"$3}' > sub_pat_chroms.tsv
for h in mat pat; do
  # 길이로 bigwig 염색체 <-> GenBank 서열 대응 (1:1 만 인정)
  awk -F'\t' 'NR==FNR{L[$2]=L[$2] (L[$2]?",":"") $1; next} {print $1"\t"$2"\t"(($2 in L)?L[$2]:"-")}' $X/asm/zfinch_v03_$h.fna.fai sub_${h}_chroms.tsv > len_match_$h.tsv
  echo "$h : bigwig 염색체 $(wc -l < sub_${h}_chroms.tsv) · 길이 1:1 대응 $(awk -F'\t' '$3!="-" && $3 !~ /,/' len_match_$h.tsv | wc -l) · 대응 없음 $(awk -F'\t' '$3=="-"' len_match_$h.tsv | wc -l)"
done > len_match_summary.txt
for h in mat pat; do
  [ -s $X/cms/zfinch_v03_$h.tbl ] || cmsearch --cut_ga --cpu 4 --noali --tblout $X/cms/zfinch_v03_$h.tbl $X/rfam/RF00001.cm $X/asm/zfinch_v03_$h.fna > $X/cms/zfinch_v03_$h.log 2>&1
done
rm -f catalog_v03.tsv
python3 $X/array_catalog.py Taeniopygia_guttata bTaeGut7_mat_GCA_048771995.1 $X/cms/zfinch_v03_mat.tbl $X/asm/zfinch_v03_mat.fna catalog_v03.tsv > catalog_v03.txt
python3 $X/array_catalog.py Taeniopygia_guttata bTaeGut7_pat_GCA_048772025.1 $X/cms/zfinch_v03_pat.tbl $X/asm/zfinch_v03_pat.fna catalog_v03.tsv >> catalog_v03.txt
touch V03_DONE
