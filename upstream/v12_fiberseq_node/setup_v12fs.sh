#!/bin/bash
# Lay out the ao_g8cut_v12 working folder for ONT Fiber-seq, with cfg = AOJV_FS.
# the v11 files are read only; everything is made as a copy (no symlinks)
set -eu
W=~/aog8v12fs
mkdir -p "$W" "$W/g8in" "$W/aoreads" "$W/tmp"

# 1) the canonical sequences: the 616 v11 canonical sequences, under the name cmp4d_v7, which is the name the code expects
if [ ! -d "$W/cmp4d_v7" ]; then
  mkdir -p "$W/cmp4d_v7"
  cp ~/cmp4d_v11/*.fa "$W/cmp4d_v7/"
fi
echo "  정본 $(ls -1 "$W/cmp4d_v7"/*.fa | wc -l)개"

# 2) the map and the summary
cp -f ~/MAP_V11.tsv "$W/MAP_G8.tsv"
cp -f ~/SUM_V11.tsv "$W/SUM_G8.tsv"
cp -f ~/aog8v10/MAP_G6C2.tsv "$W/MAP_G6C2.tsv"
cp -f ~/hap_seq_node_v9.tsv "$W/g8in/hap_seq_node_v7.tsv"
echo "  지도 $(md5sum "$W/MAP_G8.tsv" | cut -c1-8) · 요약 $(md5sum "$W/SUM_G8.tsv" | cut -c1-8)"

# 3) the references (s5s_tuc/ref)
[ -d "$W/s5s_tuc" ] || { mkdir -p "$W/s5s_tuc"; cp -r ~/s5s_tuc/ref "$W/s5s_tuc/"; }

# 4) AOJV: array-only plus junction merged; where they overlap the answer is the same, so only duplicates are removed
cat ~/fsv12_run/fsjv_*.tsv ~/jxplace_all.tsv 2>/dev/null | sort -u > "$W/AOJV_FS.tsv"
echo "  AOJV_FS $(wc -l < "$W/AOJV_FS.tsv")행 · read $(cut -f1,2 "$W/AOJV_FS.tsv" | sort -u | wc -l) · 개체 $(cut -f1 "$W/AOJV_FS.tsv" | sort -u | wc -l)"

# 5) the read fastq: <individual>.<cfg>.<data_source>.fastq.gz, with MM/ML unchanged
cut -f1 "$W/AOJV_FS.tsv" | sort -u > "$W/samples.txt"
while read -r s; do
  [ -n "$s" ] || continue
  out="$W/aoreads/$s.AOJV_FS.UCSC_ONT_FIBERSEQ.fastq.gz"
  [ -s "$out" ] && continue
  src=~/fs41/out/"$s".hit.fq
  [ -s "$src" ] || continue
  awk -F'\t' -v s="$s" 'NR==FNR{if($1==s) K[$2]=1; next}
    {h=$0; getline sq; getline p; getline q
     id=substr(h,2); split(id,a,/[ \t]/); r16=substr(a[1],1,16)
     if(r16 in K){print h; print sq; print p; print q}}' \
    "$W/AOJV_FS.tsv" "$src" | gzip -c > "$out"
done < "$W/samples.txt"
echo "  aoreads $(ls -1 "$W/aoreads"/*.fastq.gz 2>/dev/null | wc -l)개 · $(du -sh "$W/aoreads" | cut -f1)"
echo "  samples $(wc -l < "$W/samples.txt")"
