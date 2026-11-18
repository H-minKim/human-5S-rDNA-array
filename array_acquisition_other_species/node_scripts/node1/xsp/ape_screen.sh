#!/bin/bash
# Great apes: from the GenomeArk ONT BAMs (re-called with the guppy 6.3.8 methylation model, MM/ML present), keep only reads containing a 5S k-mer.
#   No re-basecalling is needed: the methylation tags are already present (GenomeArk README: "re-called using v6.3.8 with a model that included methylation").
#   Streamed: samtools view (https) -> grep -F with 8 k-mers (four 22-mers from the species 5S gene plus their reverse complements).
#   The k-mers are built from that species' own array units in the assembly (the 119 bp gene of the first unit in struct/<species>.units.fa).
#   Arguments: GenomeArk path (e.g. Pan_troglodytes/mPanTro3), k-mer file, output folder, number of concurrent streams.
#   Outputs: <flowcell>.hit.sam (tags kept) and <flowcell>.done (read count, reads with MM, pipe exit code, seconds).
# MIT License. Copyright (c) 2026 Hyungmin Kim
set -u
SPP="${1:?GenomeArk path}"
KF="${2:?k-mer file}"
export O="${3:?output folder}"
P="${4:-2}"
mkdir -p "$O"
export KF
aws s3 ls --no-sign-request --region us-east-1 --recursive "s3://genomeark/species/$SPP/genomic_data/ont/" \
  | awk '$4 ~ /[.]bam$/ {print $4}' > "$O/bams.txt"
echo "BAM $(wc -l < "$O/bams.txt")"
one() {
  key="$1"
  b=$(basename "$key" .bam)
  [ -s "$O/$b.done" ] && return 0
  s=$(date +%s)
  samtools view -@3 "https://genomeark.s3.amazonaws.com/$key" 2> "$O/$b.err" | grep -F -f "$KF" > "$O/$b.hit.sam"
  st="${PIPESTATUS[*]}"
  n=$(wc -l < "$O/$b.hit.sam")
  m=$(grep -c "MM:Z" "$O/$b.hit.sam" || true)
  printf '%s\t%s\t%s\t%s\t%s\n' "$b" "$n" "$m" "$st" "$(( $(date +%s) - s ))" > "$O/$b.done"
}
export -f one
xargs -a "$O/bams.txt" -P "$P" -I{} bash -c 'one "$1"' _ {}
echo "screening done $(date -u +%H:%M) UTC"
cat "$O"/*.done
