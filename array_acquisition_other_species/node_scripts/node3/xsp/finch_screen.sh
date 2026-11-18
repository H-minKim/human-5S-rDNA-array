#!/bin/bash
# Select only the reads containing a 5S gene k-mer from the zebra finch ONT fastq, before any alignment.
#   Streamed (S3 -> pigz -> join to one line -> grep -F), so only the candidate reads reach disk.
#   The k-mers are four 22-mers of the v0.3 maternal chr2 5S gene (cmsearch hit CM109790.1:97272500-97272617, minus strand = gene 5'->3') plus their reverse complements.
#     An array read carries hundreds of 5S copies, so an exact 22-mer match is effectively certain at R10 quality (~Q20). Dispersed 5S copies are removed by the later alignment.
#   The four fastq lines are joined with the control character 0x01, so that a tab in the header cannot break the record.
#   Arguments: a file listing the S3 keys, and the number of concurrent streams.
#   Outputs: cand/<flowcell>.cand.fq.gz and <flowcell>.done (candidate count, exit code of each stage of the pipe).
# MIT License. Copyright (c) 2026 Hyungmin Kim
set -u
L="${1:?key list}"
P="${2:-4}"
export O=~/xsp/zfont/cand
mkdir -p "$O"
K=~/xsp/zfont/k.txt
printf '%s\n' GATCTCGGAAGCTAAGCAGGGT GGGCCCGGTTAGTACTTGGATG GGAGACCGCCTGGGAATACCGG GCCGCCCTGAGAGCGGCCGGTC > "$K.fw"
rev "$K.fw" | tr ACGT TGCA > "$K.rc"
cat "$K.fw" "$K.rc" > "$K"
export K
one() {
  k="$1"
  b=$(basename "$k" .raw.fastq.gz)
  [ -s "$O/$b.done" ] && return 0
  aws s3 cp --no-sign-request --region us-east-1 "s3://genomeark/species/Taeniopygia_guttata/bTaeGut7/genomic_data/ont/$k" - 2> "$O/$b.err" \
    | pigz -dc -p 3 | paste -d $'\001\001\001' - - - - | grep -F -f "$K" | tr '\001' '\n' | gzip -1 > "$O/$b.cand.fq.gz"
  st="${PIPESTATUS[*]}"
  n=$(( $(zcat "$O/$b.cand.fq.gz" | wc -l) / 4 ))
  printf '%s\t%s\t%s\n' "$b" "$n" "$st" > "$O/$b.done"
}
export -f one
xargs -a "$L" -P "$P" -I{} bash -c 'one "$1"' _ {}
echo "screening finished $(date -u +%H:%M) UTC"
cat "$O"/*.done
