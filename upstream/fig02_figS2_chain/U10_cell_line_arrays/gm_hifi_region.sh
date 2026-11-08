#!/bin/bash
# Build the GM12878 HiFi assembly input properly, using the same gate as K562
#  Why change (measured): coordinate region extraction 715 reads vs 5S unit gate 2,373 reads, overlap only 346
#    - 369 only in the region = **pure flank reads** without 5S k-mers -> anchors at both array ends, needed for assembly, but the unit gate misses them
#    - 2,027 only in the full set = missed by region extraction (though many may be dispersed 5S copies)
#  => the correct gate is not the unit but BBDuk k=31 mkf=0.2 against the **whole region sequence** (array +-115 kb)
#     (= the same gate K562 used, `p11_k562_region.sh`. With mkf=0.2, 20 % of the read k-mers must match, so
#      a dispersed-copy read that has only the 120 bp 5S gene inside a 15 kb read is filtered out at about 1 %)
#  Coordinates are not used, so unaligned reads and reads primarily aligned outside the region also come in
set -u
W=~/p11_encode/gmhifi_region; mkdir -p $W; cd $W
BB=$HOME/opt/bbmap/bbduk.sh
U="https://giab.s3.amazonaws.com/data/NA12878/PacBio_CCS_15kb_20kb_chemistry2/CHM13v2.0/GIAB_5mC_CpG/HG001.chm13.deepvariant.haplotagged.bam"
log(){ echo "[$(date +%H:%M)] $*" >> $W/RUN.log; }
# region reference = the same coordinates as K562 (CHM13 chr1:227,626,780-228,142,150, array +-115 kb)
if [ ! -s $W/region_ref.fa ]; then
  cp ~/p11_encode/ref/array_pm50k.fa $W/region_ref.fa 2>/dev/null
  log "구간 참조 = array_pm50k.fa ($(grep -v '>' $W/region_ref.fa | tr -d '\n' | wc -c) bp)"
fi
log "HG001 HiFi 전량 → 구간 서열 관문(k=31 mkf=0.2) 시작"
s=$(date +%s)
samtools fastq -T '*' -F 0x900 "$U" 2> $W/samtools.err \
  | $BB -Xmx10g in=stdin.fq interleaved=f outm=$W/gm_region_full.fastq.gz \
      ref=$W/region_ref.fa k=31 mkf=0.2 threads=12 overwrite=t > $W/bbduk.log 2>&1
rc=$?; e=$(date +%s)
log "관문 종료 rc=$rc · $(( (e-s)/60 )) 분"
grep -E "^Input:|^Contaminants:" $W/bbduk.log >> $W/RUN.log
if [ -s $W/gm_region_full.fastq.gz ]; then
  N=$(( $(zcat $W/gm_region_full.fastq.gz | wc -l) / 4 ))
  log "구간관문 전량 read $N (좌표 구간뽑기 1차는 715) · $(zcat $W/gm_region_full.fastq.gz | awk 'NR%4==2{s+=length($0); if(length($0)>m)m=length($0)}END{printf "%.1f Mb · 최장 %d", s/1e6, m}')"
  zcat $W/gm_region_full.fastq.gz | awk 'NR%4==1{print substr($1,2)}' | sort -u > $W/regionfull.ids
  log "  고유 read $(wc -l < $W/regionfull.ids)"
else log "산출 없음"; exit 9; fi
# HiFi-only assembly, same conditions as round 1 (119/91), only the input changed for comparison
A=$W/asm; mkdir -p $A; cd $A; rm -f *.bin *.gfa *.fa
~/hifiasm -o HG001full -t 14 $W/gm_region_full.fastq.gz > hifiasm.log 2>&1
rc=$?
for g in *p_ctg.gfa; do [ -e "$g" ] && awk '/^S/{print ">"$2"\n"$3}' "$g" > "${g%.gfa}.fa"; done
log "HiFi 전용 재조립 rc=$rc"
for fa in $A/HG001full.bp.hap1.p_ctg.fa $A/HG001full.bp.hap2.p_ctg.fa; do
  [ -e "$fa" ] && log "  $(basename $fa) : $(awk '/^>/{n++;next}{L[n]+=length($0)}END{for(i=1;i<=n;i++) if(L[i]>100000) printf "%d bp ", L[i]}' $fa)"
done
log GMHIFI_REGION_DONE
