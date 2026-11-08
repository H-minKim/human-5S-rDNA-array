#!/bin/bash
# Re-run at threshold 5; the earlier outputs aog8v12fs2, fsfinal_v12b, aog8v12pb and pbfinal_v12 are not touched
#   명령은 옛 실행과 같고 바뀐 것 : MINSITE=5 · 출력 폴더 이름
cd ~; H=$(hostname)
for d in ~/aog8v12fs3 ~/fsfinal_v12c ~/aog8v12pb2 ~/pbfinal_v12b; do [ -e $d ] && { echo "이미 있음 $d"; exit 1; }; done
# ── 1. 합치기
cat ~/fsv12_run/fsjv_*.tsv > /tmp/ont_ao_s5.tsv
MINSITE=5 KEYLEN=16 python3 ~/merge_aojv_v2.py /tmp/ont_ao_s5.tsv --jx ~/jxplace_all.tsv --vote ~/slotvote_qc2.tsv --out ~/AOJV_FS_merged_s5.tsv --log ~/AOJV_FS_merge_s5.log
cat ~/pbv12_run/pfjv_*.tsv 2>/dev/null | cut -f1-6 > /tmp/pb_ao_s5.tsv
MINSITE=5 KEYLEN=0 python3 ~/merge_aojv_v2.py /tmp/pb_ao_s5.tsv --jx ~/jxplace_pb2.tsv --vote ~/slotvote_pb.tsv --out ~/AOJV_PB_merged_s5.tsv --log ~/AOJV_PB_merge_s5.log
echo "MERGE_DONE"
# ── 2. ONT 작업폴더 (aog8v12fs2 와 같은 구성)
W=~/aog8v12fs3; mkdir -p $W/g8in $W/aoreads $W/tmp
cp -r ~/aog8v12fs/cmp4d_v7 $W/; cp -r ~/aog8v12fs/s5s_tuc $W/
for f in MAP_G8.tsv SUM_G8.tsv MAP_G6C2.tsv; do cp -f ~/aog8v12fs/$f $W/; done
cp -f ~/aog8v12fs/g8in/hap_seq_node_v7.tsv $W/g8in/
cp -f ~/AOJV_FS_merged_s5.tsv $W/AOJV_FS.tsv
for f in AOJV3 AOJV5 AOJV8 AOJV_GIAB AOJV_PB; do : > $W/$f.tsv; done
cut -f1 $W/AOJV_FS.tsv | sort -u > $W/samples.txt
for f in ~/aog8v12fs/aoreads/*.fastq.gz; do ln -sf "$f" $W/aoreads/ 2>/dev/null; done
# ── 3. PacBio 작업폴더 (aog8v12pb 와 같은 구성 · aoreads 는 새 합본으로 다시 뽑음)
P=~/aog8v12pb2; mkdir -p $P/g8in $P/aoreads $P/tmp
cp -r ~/aog8v12fs/cmp4d_v7 $P/
cp -f ~/MAP_V11.tsv $P/MAP_G8.tsv; cp -f ~/SUM_V11.tsv $P/SUM_G8.tsv
cp -f ~/aog8v10/MAP_G6C2.tsv $P/MAP_G6C2.tsv; cp -f ~/hap_seq_node_v9.tsv $P/g8in/hap_seq_node_v7.tsv
mkdir -p $P/s5s_tuc; cp -r ~/s5s_tuc/ref $P/s5s_tuc/
cp -f ~/AOJV_PB_merged_s5.tsv $P/AOJV_PB.tsv
for f in AOJV3 AOJV5 AOJV8 AOJV_GIAB AOJV_FS; do : > $P/$f.tsv; done
cut -f1 $P/AOJV_PB.tsv | sort -u > $P/samples.txt
while read -r s; do
  out=$P/aoreads/$s.AOJV_PB.PACBIO_FIBERSEQ.fastq.gz
  src=~/pbfs_bysample/$s.hit.fq; [ -s "$src" ] || continue
  awk -F"\t" -v s="$s" 'NR==FNR{if($1==s) K[$2]=1; next}
    {h=$0; getline sq; getline p; getline q; id=substr(h,2); split(id,a,/[ \t]/)
     if(a[1] in K){print h; print sq; print p; print q}}' $P/AOJV_PB.tsv "$src" | gzip -c > "$out"
done < $P/samples.txt
echo "SETUP_DONE ONT개체 $(wc -l < $W/samples.txt) aoreads $(ls $W/aoreads | wc -l) · PB개체 $(wc -l < $P/samples.txt) aoreads $(ls $P/aoreads | wc -l)"
# ── 4. 절단·재호출 (옛 실행과 같은 명령)
sed "s|~/aog8v10/sam|~/aog8v12fs3/sam|" ~/recall_edge_v10.py > ~/recall_edge_fs3.py
sed "s|~/aog8v10/sam|~/aog8v12fs3/sam|" ~/recall_mid_v10.py  > ~/recall_mid_fs3.py
sed "s|~/aog8v10/sam|~/aog8v12fs3/sam|" ~/recall_pad2_v10.py > ~/recall_pad2_fs3.py
sed "s|~/aog8v10/sam|~/aog8v12pb2/sam|" ~/recall_edge_v10.py > ~/recall_edge_pb2.py
sed "s|~/aog8v10/sam|~/aog8v12pb2/sam|" ~/recall_mid_v10.py  > ~/recall_mid_pb2.py
sed "s|~/aog8v10/sam|~/aog8v12pb2/sam|" ~/recall_pad2_v10.py > ~/recall_pad2_pb2.py
( python3 ~/ao_g8cut_v12.py $W $W/samples.txt 8 > $W/cut.log 2>&1
  cd $W && python3 ~/recall_edge_fs3.py AO AO_G8_UNITS.tsv AO_G8_UNITS.edge.tsv &&
  python3 ~/recall_mid_fs3.py AO AO_G8_UNITS.edge.tsv AO_G8_UNITS.mid0.tsv &&
  PAD=100 CTX=200 python3 ~/recall_pad2_fs3.py AO AO_G8_UNITS.mid0.tsv AO_G8_UNITS.pad.tsv
  mkdir -p ~/fsfinal_v12c
  python3 ~/fix_relabel_only.py $W/AO_G8_UNITS.pad.tsv $W/AO_G8_READS.tsv ~/fsfinal_v12c/AO_G8_UNITS_$H.final.tsv
  python3 ~/scan_shift_all.py ~/fsfinal_v12c/AO_G8_UNITS_$H.final.tsv ~/MAP_V11.tsv ~/fsfinal_v12c/SHIFT_$H.tsv
  echo FS3_DONE ) > ~/fs3run.log 2>&1 &
( python3 ~/ao_g8cut_v12.py $P $P/samples.txt 8 > $P/cut.log 2>&1
  cd $P && python3 ~/recall_edge_pb2.py AO AO_G8_UNITS.tsv AO_G8_UNITS.edge.tsv &&
  python3 ~/recall_mid_pb2.py AO AO_G8_UNITS.edge.tsv AO_G8_UNITS.mid0.tsv &&
  PAD=100 CTX=200 python3 ~/recall_pad2_pb2.py AO AO_G8_UNITS.mid0.tsv AO_G8_UNITS.pad.tsv
  mkdir -p ~/pbfinal_v12b
  python3 ~/fix_relabel_only.py $P/AO_G8_UNITS.pad.tsv $P/AO_G8_READS.tsv ~/pbfinal_v12b/AO_G8_UNITS_$H.final.tsv
  python3 ~/scan_shift_all.py ~/pbfinal_v12b/AO_G8_UNITS_$H.final.tsv ~/MAP_V11.tsv ~/pbfinal_v12b/SHIFT_$H.tsv
  echo PB2_DONE ) > ~/pb2run.log 2>&1 &
wait
echo S5_ALL_DONE >> ~/fs3run.log
