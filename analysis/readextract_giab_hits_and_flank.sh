#!/bin/bash
# 개체별 hits 합본 + CHM13 chr1 PROX/DIST 40 kb 창으로 junction 확정
#   단계 : ① read 추출 (GIAB 합치기·flank)
# GIAB + HG008 — HPRC ONT 정식 경로(ont_flank2.sh) 그대로. 입력만 다르다.
#   ① ont_flank_extract.py  ② probe5k → 40kb 창+tandem 정렬  ③ ont_flank_keep.py(거리<=50·mapq>=20·strand·탠덤제외)
set -u
W=$HOME/giab_flank; rm -rf $W; mkdir -p $W/{junc,jd,verify,log}; cd $W || exit 9
R=$HOME/ectopic_sp/chm13v2.0.fa
TAN=$HOME/s5s_tuc/ref/tandem40.fa
BD=$HOME/giab_tuc/out
HD=$HOME/giab_hits; mkdir -p $HD
# hits 를 개체별 단일 파일로 (HPRC 의 ~/s5s_ont/hits 형식과 동일하게)
for s in hg001 hg002 hg003 hg004; do
  [ -s $HD/$s.hits.fastq.gz ] || cat $HOME/giab5s_seq/out/${s}_*.hits.fastq.gz > $HD/$s.hits.fastq.gz 2>/dev/null
done
for s in HG008-N-D HG008-N-P; do
  [ -s $HD/$s.hits.fastq.gz ] || cat $HOME/hg008/out/${s}_*.hits.fastq.gz > $HD/$s.hits.fastq.gz 2>/dev/null
done
[ -s $HD/HG008-T.hits.fastq.gz ] || cp $HOME/hg008/out/HG008-T.hits.fastq.gz $HD/ 2>/dev/null
echo "  hits: $(ls $HD/*.hits.fastq.gz | wc -l) 개체"
samtools faidx $R chr1:227704289-227744288 | sed "1s/.*/>PROX_JUNC_AT_END_len40000/" > junc/ref_prox.fa
samtools faidx $R chr1:228026263-228066262 | sed "1s/.*/>DIST_JUNC_AT_START_len40000/" > junc/ref_dist.fa
cat junc/ref_prox.fa junc/ref_dist.fa $TAN > junc/ref3.fa
python3 $HOME/ont_flank_extract.py "$BD" "$HD" "$W" 5000 || exit 9
echo "  probe $(grep -c '>' junc/probe5k.fa) · full $(grep -c '>' junc/flank_full.fa)"
minimap2 -x map-ont -t 8 -c --secondary=no junc/ref3.fa junc/probe5k.fa > junc/probe5k.paf 2>/dev/null
python3 $HOME/ont_flank_keep.py "$W" || exit 9
echo "  확정 $(grep -c '>' jd/flank_confirmed_full.fa)"
awk -F"\t" '{c[$2]++} END{for(x in c) printf "    %s: %d\n", x, c[x]}' jd/keep.tsv
# 2단계
minimap2 -x map-ont -a -t 8 --secondary=no $HOME/allspan/ref3.fa jd/flank_confirmed_full.fa 2>mm.log | samtools sort -@2 -o all_flank.bam -
samtools index all_flank.bam
mkdir -p kids
samtools view all_flank.bam | cut -f1 | cut -d"|" -f1 | sort -u > samples.txt
samtools view -H all_flank.bam > hdr.sam; samtools view all_flank.bam > body.sam
while read s; do [ -n "$s" ] || continue
  { cat hdr.sam; grep "^${s}|" body.sam; } | samtools sort -o kids/$s.bam - 2>/dev/null; done < samples.txt
rm -f body.sam
python3 $HOME/hap_all_pri.py kids hap_reads_giab.tsv hap_summary_giab.tsv > hap.log 2>&1
echo "  $(tail -1 hap.log) · 배정 read $(( $(wc -l < hap_reads_giab.tsv) - 1 )) · 개체 $(tail -n+2 hap_reads_giab.tsv|cut -f1|sort -u|wc -l)"
awk -F"\t" 'NR>1{c[$1"/"$2]++} END{for(x in c) printf "    %-16s %d\n", x, c[x]}' hap_reads_giab.tsv | sort
touch GIABOK_DONE
