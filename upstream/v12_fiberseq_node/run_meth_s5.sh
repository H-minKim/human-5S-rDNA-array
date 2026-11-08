#!/bin/bash
# 문턱 5 판 m6A·CpG 재추출 — 코드 fsmeth_v12c.py (f1b83605) · 옛 폴더는 건드리지 않음
export PATH=$PATH:~/miniconda3/bin:~/bin; cd ~; H=$(hostname)
for d in ~/fsmeth_v12c ~/pbmeth_v12d ~/pbft_v12b; do [ -e $d ] && { echo "이미 있음 $d"; exit 1; }; done
D=~/aog8v12fs3 CFG=AOJV_FS LABEL=UCSC_ONT_FIBERSEQ OUT=~/fsmeth_v12c FINAL="~/fsfinal_v12c/AO_G8_UNITS_$H.final.tsv" python3 ~/fsmeth_v12c.py > ~/fsmeth_v12c.run.log 2>&1
D=~/aog8v12pb2 CFG=AOJV_PB LABEL=PACBIO_FIBERSEQ OUT=~/pbmeth_v12d FINAL="~/pbfinal_v12b/AO_G8_UNITS_$H.final.tsv" python3 ~/fsmeth_v12c.py > ~/pbmeth_v12d.run.log 2>&1
# ft 재호출판 (원자료 m6A 없는 두 개체가 이 노드에 있을 때만)
for s in HG03672 HG02624; do
  [ -f ~/aog8v12pb2/sam/$s.AOJV_PB.sam ] || continue
  mkdir -p ~/pbft_v12b; cp -r ~/pbft_v12/rg ~/pbft_v12b/ 2>/dev/null
  SRCD=~/aog8v12pb2/sam FTW=~/pbft_v12b python3 ~/ftprep_v12b.py $s >> ~/pbft_v12b.run.log 2>&1
  D=~/aog8v12pb2 SAMD=~/pbft_v12b/sam CFG=AOJV_PB LABEL=PACBIO_FIBERSEQ_FT OUT=~/pbft_v12b/meth FINAL="~/pbfinal_v12b/AO_G8_UNITS_$H.final.tsv" SMPS=$s python3 ~/fsmeth_v12c.py >> ~/pbft_v12b.run.log 2>&1
done
echo METH_S5_DONE >> ~/pbmeth_v12d.run.log
