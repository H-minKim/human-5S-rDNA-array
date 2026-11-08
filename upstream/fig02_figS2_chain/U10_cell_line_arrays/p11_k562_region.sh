#!/bin/bash
# K562 HiFi BAM (local) → extract reads in the CHM13 5S array ±120 kb region (for assembly, reaching the 50 kb flank probes)
#   BBDuk k=31, only when at least 20 % of read k-mers match the reference (mkf=0.2) - keeps out reads pulled in by a single Alu piece in the flank
#   tags (MM, ML, m6A) preserved
set -u
W=~/p11_encode/k562hifi_dl; cd $W
[ -e region_ref.fa ] || samtools faidx ~/ectopic_sp/chm13v2.0.fa chr1:227626780-228142150 > region_ref.fa
samtools quickcheck k562_hifi.bam && echo "BAM 무결성 OK" || { echo "BAM 손상"; exit 1; }
samtools fastq -T '*' k562_hifi.bam 2> fastq.err | ~/opt/bbmap/bbduk.sh -Xmx12g in=stdin.fq interleaved=f outm=k562_region.fastq.gz ref=region_ref.fa k=31 mkf=0.2 threads=12 > region_bbduk.log 2>&1
grep -E "^Input|^Contam" region_bbduk.log
echo REGIONDONE
