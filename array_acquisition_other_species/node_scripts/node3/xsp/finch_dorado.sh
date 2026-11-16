#!/bin/bash
# Zebra finch: only the 5S array reads are re-basecalled with CPU dorado, and only for the 5 kHz flowcells, because dorado 2.1.2 has no model for 4 kHz.
#   dorado basecaller hac@v5.0.0,5mCG_5hmCG (@v3), --read-ids, with the pod5 files read partially over mount-s3.
#     Measured on CPU (same 20 reads, 8 cores): hac@v5.0.0 2.64e4 samples/s against the auto-selected hac@v6.0.0 at 9.5e3, so v5 is 2.8x faster.
#   -> samtools fastq -T MM,ML -> minimap2 -y -ax map-ont --secondary=no (the v0.3 diploid assembly) -> aligned BAM with MM/ML.
#      (The same alignment as zebrafish; no MAPQ filter.)
#   Arguments: flowcell list file, taskset core range (e.g. 0-15), [read id folder], [model = hac,5mCG_5hmCG].
#   Outputs: bc/<fc>.ubam (the re-basecalled reads), bc/<fc>.bam (+bai), bc/<fc>.done (exit code, read count, reads with MM, seconds).
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
Z=~/xsp/zfont
A=~/xsp/asm
D=~/xsp/tools/dorado-2.1.2-linux-x64/bin/dorado
MNT=~/xsp/ms3mnt
L="${1:?flowcell 목록}"
CORES="${2:?코어 범위}"
IDD="${3:-$Z/sel}"
MODEL="${4:-hac@v5.0.0,5mCG_5hmCG}"
BC="${BC:-$Z/bc}"
mkdir -p "$BC" ~/xsp/dmodels
for b in $(cat "$L"); do
  [ -s $BC/$b.done ] && continue
  p=$(ls $MNT/*/raw/uncorrected/$b.pod5 2> /dev/null | head -1)
  [ -n "$p" ] || { printf '%s\tNO_POD5\n' "$b" > $Z/bc/$b.done; continue; }
  [ -s "$IDD/$b.ids" ] || { printf '%s\tNO_IDS\n' "$b" > $Z/bc/$b.done; continue; }
  s=$(date +%s)
  taskset -c "$CORES" $D basecaller "$MODEL" "$p" --read-ids "$IDD/$b.ids" --device cpu --models-directory ~/xsp/dmodels > $BC/$b.ubam 2> $BC/$b.dorado.log
  rc=$?
  # Primary only (-F 0x900, as for zebrafish), so that secondary alignments cannot mix in.
  samtools fastq -T MM,ML $BC/$b.ubam 2> /dev/null | minimap2 -y -ax map-ont -t 6 --secondary=no $A/zfinch_v03_dip.map-ont.mmi - 2> /dev/null | samtools view -b -F 0x900 - | samtools sort -@2 -o $BC/$b.bam - && samtools index $BC/$b.bam
  n=$(samtools view -c $BC/$b.ubam)
  m=$(samtools view $BC/$b.ubam | grep -c "MM:Z")
  printf '%s\t%s\t%s\t%s\t%s\n' "$b" "$rc" "$n" "$m" "$(( $(date +%s) - s ))" > $BC/$b.done
done
echo "재콜링 끝 $(date -u +%H:%M) UTC"
cat $Z/bc/*.done
