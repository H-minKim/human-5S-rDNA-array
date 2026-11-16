#!/bin/bash
# Other species: cmsearch --cut_ga with Rfam RF00001 (5S rRNA) against each assembly.
#   Assemblies are processed as they arrive, smallest first, and anything already done is skipped.
set -u
X=$HOME/xsp
mkdir -p $X/cms
for nm in arabidopsis_ColCC zebrafish_GRCz12tu zfinch_bTaeGut7_mat zfinch_bTaeGut7_pat maize_Mo17_T2T zebrafish_GRCz11_ens; do
  f=$X/asm/$nm.fna
  [ "$nm" = zebrafish_GRCz11_ens ] && f=$X/asm/$nm.fa
  while [ ! -s "$f.fai" ]; do
    if [ -f $X/log/SETUP_DONE ] && [ ! -s "$f.fai" ]; then
      echo "[$(date -u +%H:%M)] $nm 어셈블리 없음 - 건너뜀" >> $X/cms/progress.log
      continue 2
    fi
    sleep 30
  done
  [ -s $X/cms/$nm.tbl ] && continue
  cmsearch --cut_ga --cpu 6 --noali --tblout $X/cms/$nm.tbl $X/rfam/RF00001.cm "$f" > $X/cms/$nm.log 2>&1
  echo "[$(date -u +%H:%M)] $nm rc=$? · hit $(grep -vc '^#' $X/cms/$nm.tbl)" >> $X/cms/progress.log
done
touch $X/cms/CMS_DONE
