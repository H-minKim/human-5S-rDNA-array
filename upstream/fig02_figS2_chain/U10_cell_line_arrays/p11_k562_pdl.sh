#!/bin/bash
# Parallel ranged download of the original submitted K562 HiFi BAM (32 connections, 256MB chunks) → concatenate → check md5 (ENA submitted_md5)
set -u
W=~/p11_encode/k562hifi_dl; mkdir -p $W/parts; cd $W
U=https://ftp.sra.ebi.ac.uk/vol1/run/ERR176/ERR17608832/m84259_250328_214634_s2.hifi_reads.bc2009.bam
SIZE=$(curl -sI "$U" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r')
MD5=$(curl -fs "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERR17608832&result=read_run&fields=submitted_md5&format=tsv" | tail -1 | cut -f2 | cut -d';' -f1)
echo "크기 $SIZE · ENA md5 $MD5" > info.txt
CH=268435456; N=$(( (SIZE + CH - 1) / CH ))
seq 0 $((N-1)) | xargs -P 32 -I{} bash -c 'i={}; s=$((i*'$CH')); e=$((s+'$CH'-1)); f=parts/p$(printf %05d $i); [ -e $f.ok ] && exit 0; for t in 1 2 3 4 5; do curl -s -r $s-$e -o $f "'$U'" && [ $(stat -c %s $f) -eq $(( e < '$SIZE' ? '$CH' : '$SIZE'-s )) ] && { touch $f.ok; exit 0; }; sleep 5; done; echo "FAIL $i" >> fail.log'
ls parts/*.ok | wc -l >> info.txt
cat $(ls parts/p????? | sort) > k562_hifi.bam
echo "이어붙인 크기 $(stat -c %s k562_hifi.bam) · md5 $(md5sum k562_hifi.bam | cut -c1-32)" >> info.txt
echo PDLDONE >> info.txt
