#!/bin/bash
# Arabidopsis: pod5 slices holding only the array reads -> CPU dorado (hac@v5.0.0,5mC_5hmC = all contexts CG, CHG, CHH) -> realignment to Col-CC (primary only).
#   Realignment = samtools fastq -T MM,ML | minimap2 -y -ax map-ont --secondary=no (the same as zebrafish and zebra finch); no MAPQ filter.
#   Arguments: pod5 list file (one path per line), taskset core range, label, [model].
#   Outputs: bc/<label>.ubam, bc/<label>.bam (+bai), bc/<label>.done (exit code, reads, reads with MM, seconds).
# MIT License · Copyright (c) 2026 Hungmin Kim
set -u
H=~/xsp/ath
A=~/xsp/asm
D=~/xsp/tools/dorado-2.1.2-linux-x64/bin/dorado
L="${1:?pod5 목록}"
CORES="${2:?코어 범위}"
TAG="${3:?이름표}"
# The model is given as a path, not a name: automatic selection fails with "Failed to determine sequencing chemistry" (the pod5 files measure 5 kHz, FLO-PRO114M, SQK-LSK114).
M=~/xsp/dmodels
SIMPLEX="${4:-dna_r10.4.1_e8.2_400bps_hac@v5.0.0}"
MODS="${5:-dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mC_5hmC@v3}"
B=$H/bc
mkdir -p $B/$TAG.in $M
for m in "$SIMPLEX" "$MODS"; do [ -d $M/$m ] || $D download --model "$m" --models-directory $M > /dev/null 2>&1; done
for f in $(cat "$L"); do ln -sf "$f" $B/$TAG.in/; done
s=$(date +%s)
taskset -c "$CORES" $D basecaller $M/$SIMPLEX $B/$TAG.in --modified-bases-models $M/$MODS --device cpu > $B/$TAG.ubam 2> $B/$TAG.dorado.log
rc=$?
samtools fastq -T MM,ML $B/$TAG.ubam 2> /dev/null | minimap2 -y -ax map-ont -t 6 --secondary=no $A/arabidopsis_ColCC.map-ont.mmi - 2> /dev/null | samtools view -b -F 0x900 - | samtools sort -@2 -o $B/$TAG.bam - && samtools index $B/$TAG.bam
n=$(samtools view -c $B/$TAG.ubam)
m=$(samtools view $B/$TAG.ubam | grep -c "MM:Z")
printf '%s\t%s\t%s\t%s\t%s\n' "$TAG" "$rc" "$n" "$m" "$(( $(date +%s) - s ))" > $B/$TAG.done
echo "재콜링 끝 $TAG $(date -u +%H:%M) UTC"
cat $B/$TAG.done
