#!/bin/bash
# On the node: the placed-read SAM from the reproduction release -> a sorted, indexed BAM holding the placed reads only.
#   For ARRAY_ONLY the QUAL is filled in from the original fastq, by qualfill.py, and only where the sequence matches character for character.
#   The junction reads were aligned from FASTQ and already carry a QUAL.
#   The MM/ML tags come straight from the original fastq. The data source, srcmap, goes into the file name.
#   One output line gives: node, job, individual, source, record count, records with a QUAL, records with MM.
set -u
mkdir -p ~/placed_bam /tmp/qf
H=$(hostname)
SM=~/srcmap_$H.tsv
[ -s "$SM" ] || { echo -e "$H\t-\t-\tNO_SRCMAP"; exit 1; }
for d in rjrepro rjrepro_c2 aorepro3 aorepro5 aorepro8 aoreproG; do
  case $d in
    rjrepro) j=JX ;; rjrepro_c2) j=JX_C2 ;; aorepro3) j=AOJV3 ;; aorepro5) j=AOJV5 ;; aorepro8) j=AOJV8 ;; aoreproG) j=AOJVG ;;
  esac
  for s in ~/$d/*.placed.sam; do
    [ -s "$s" ] || continue
    b=$(basename "$s" .placed.sam)
    lab=$(awk -F'\t' -v j="$j" -v b="$b" '$2==j && $3==b{print $4}' "$SM" | head -1)
    fq=$(awk -F'\t' -v j="$j" -v b="$b" '$2==j && $3==b{print $5}' "$SM" | head -1 | sed "s#~#$HOME#g")
    [ -n "$lab" ] || lab=UNKNOWN
    out=~/placed_bam/${b}.${j}.${lab}.bam
    in="$s"
    case $d in
      aorepro*)
        if ! python3 ~/qualfill.py "$s" "$fq" /tmp/qf/$b.$j.sam > /tmp/qf/log_$b.$j.txt 2>&1; then
          echo -e "$H\t$j\t$b\t$lab\tQUALFILL_FAIL"
          continue
        fi
        in=/tmp/qf/$b.$j.sam ;;
    esac
    if samtools sort -@2 -o "$out" "$in" 2>/dev/null && samtools index "$out"; then
      n=$(samtools view -c "$out")
      q=$(samtools view "$out" | awk -F'\t' '$11!="*"' | wc -l)
      m=$(samtools view "$out" | grep -c "MM:Z")
      echo -e "$H\t$j\t$b\t$lab\t$n\t$q\t$m"
    else
      echo -e "$H\t$j\t$b\t$lab\tSORT_FAIL"
    fi
    rm -f /tmp/qf/$b.$j.sam
  done
done
