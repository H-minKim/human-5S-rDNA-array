#!/bin/bash
# Zebrafish 5S read and methylation extraction (run where CPU is free)
#   ① GRCz12tu(NCBI GCA_049306965.1) · modkit 준비  ② zf_extract.sh 실행
# The two BED files must be uploaded to ~/xsp/zf/ beforehand (made by the preparation script)
set -u
X=$HOME/xsp
mkdir -p $X/asm $X/tools $X/zf $X/log
if [ ! -s $X/asm/zebrafish_GRCz12tu.fna.fai ]; then
  curl -sL --retry 3 -o $X/asm/grcz12.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_049306965.1/download?include_annotation_type=GENOME_FASTA"
  unzip -o -q $X/asm/grcz12.zip -d $X/asm/grcz12.d && cat $(find $X/asm/grcz12.d -name "*.fna") > $X/asm/zebrafish_GRCz12tu.fna && rm -rf $X/asm/grcz12.d
  samtools faidx $X/asm/zebrafish_GRCz12tu.fna
fi
echo "GRCz12tu 서열 $(wc -l < $X/asm/zebrafish_GRCz12tu.fna.fai) · md5(zip) $(md5sum $X/asm/grcz12.zip | cut -c1-12)"
if [ ! -x $X/tools/modkit ]; then
  curl -sL -o $X/tools/modkit.tgz https://github.com/nanoporetech/modkit/releases/download/v0.6.4/modkit_v0.6.4_u16_x86_64.tar.gz
  tar -xzf $X/tools/modkit.tgz -C $X/tools && find $X/tools -name modkit -type f -exec mv {} $X/tools/modkit \; 2>/dev/null
  chmod +x $X/tools/modkit
fi
echo "modkit $($X/tools/modkit --version 2>&1 | head -1)"
bash $X/zf_extract.sh $X/zf/grcz11_5S_regions.bed $X/zf/grcz12_5S_arrays.bed
