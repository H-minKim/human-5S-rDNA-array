#!/bin/bash
# Stops only the slow single-connection wget that downloads ATRun3_SUP; only processes named wget are targeted
for p in $(ps -eo pid=,comm=,args= | awk '$2=="wget" && index($0,"ATRun3_SUP_all_reads")>0 {print $1}'); do
  kill "$p" && echo "wget $p 멈춤"
done
rm -f ~/xsp/ath/ATRun3_SUP_all_reads.fastq.gz ~/xsp/ath/wget_sup.log
