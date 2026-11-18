#!/bin/bash
# A biological replicate for the eggs: the Potok 2013 series (PRJNA193544) is aligned.
#   It is the companion paper to Jiang 2013 (PRJNA188516), published in the same issue of Cell,
#   and its egg BioSample is SAMN01991170, different from SAMN02117044.
#   A different laboratory, a different pool of eggs and a different library, so a genuine biological replicate.
#   The same study holds sperm, 2-16 cell, 64, 256, sphere and muscle, so the whole series serves as a control.
#
#   Only runs that pass the conversion check (zf_qc_eggs2.sh) are processed.
#   8M read pairs at a time, temporary files deleted after each, stopping if less than 25 GB of disk is free.
set -u
cd ~/zfwgbs || exit 1
export PATH=$HOME/zfwgbs/Bismark-0.24.2:$PATH
mkdir -p sra bam logs tmp
N=8000000

guard () {
  F=$(df --output=avail -BG / | tail -1 | tr -dc '0-9')
  if [ "$F" -lt 25 ]; then echo "STOP: free disk ${F}G"; exit 1; fi
}

one () {
  R=$1; D=$2
  [ -f bam/${R}_done ] && { echo "SKIP $R ($D)"; return; }
  guard
  [ -f sra/${R}_1.fastq.gz ] || fastq-dump -X $N --split-files --gzip -O sra "$R" > logs/${R}.dl.log 2>&1
  [ -f sra/${R}_1.fastq.gz ] || { echo "DL_FAIL $R"; return; }
  # These are 2013 data and some runs were submitted single-end, so a run with no mate file is skipped.
  [ -f sra/${R}_2.fastq.gz ] || { echo "SINGLE_END $R ($D): skipped"; rm -f sra/${R}_*.fastq.gz; return; }
  echo "alignment start $R ($D) $(date -u +%H:%M)"
  bismark --genome ref --parallel 2 -p 3 --bam --temp_dir tmp -o bam \
    -1 sra/${R}_1.fastq.gz -2 sra/${R}_2.fastq.gz > logs/${R}.align.log 2>&1
  B=$(ls bam/${R}_1*bismark*pe.bam 2>/dev/null | head -1)
  if [ -n "$B" ]; then
    samtools sort -@ 6 -o bam/${R}.sorted.bam "$B" && samtools index bam/${R}.sorted.bam && rm -f "$B" && touch bam/${R}_done
    echo "OK $R  $(grep 'Mapping efficiency' logs/${R}.align.log | head -1)"
  else
    echo "ALIGN_FAIL $R"; tail -3 logs/${R}.align.log
  fi
  rm -f sra/${R}_*.fastq.gz tmp/* 2>/dev/null
  df -h / | tail -1
}

# In order of importance: the eggs are the target, with muscle and sperm as controls within the same study.
one SRR800059 알_Potok2013
one SRR800081 근육성체_Potok2013
one SRR800056 정자_Potok2013
one SRR800077 sphere_Potok2013
one SRR800061 2-16세포_Potok2013
one SRR800066 64세포_Potok2013
one SRR800070 256세포_Potok2013
echo POTOK_DONE
