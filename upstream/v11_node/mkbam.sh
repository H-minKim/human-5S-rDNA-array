#!/bin/bash
# 노드에서 : 재현판이 낸 배치 read SAM -> 정렬·색인 BAM (배치 read 만)
#   ARRAY_ONLY 는 원 fastq 의 QUAL 을 채운다 (qualfill.py · 서열 글자 대조 통과 시만)
#   junction 은 FASTQ 로 정렬했으므로 QUAL 이 이미 있다
#   MM/ML 은 원 fastq 태그 그대로. 파일 이름에 데이터 출처(srcmap)를 넣는다
#   출력 한 줄 = 노드 · 작업 · 개체 · 출처 · 기록 수 · QUAL 있는 기록 · MM 있는 기록
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
