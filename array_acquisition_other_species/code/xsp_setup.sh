#!/bin/bash
# Other species 5S arrays and methylation: prepares the work node
#   The assembly, the Rfam 5S model and modkit are downloaded under xsp. Existing folders are not touched.
set -u
X=$HOME/xsp
mkdir -p $X/asm $X/tools $X/rfam $X/log
cd $X || exit 1
log() { echo "[$(date -u +%H:%M:%S)] $*" >> $X/log/setup.log; }
log "start"
# modkit
if [ ! -x $X/tools/modkit ]; then
  curl -sL -o $X/tools/modkit.tgz https://github.com/nanoporetech/modkit/releases/download/v0.6.4/modkit_v0.6.4_u16_x86_64.tar.gz
  tar -xzf $X/tools/modkit.tgz -C $X/tools && find $X/tools -name modkit -type f -exec mv {} $X/tools/modkit \; 2>/dev/null
  chmod +x $X/tools/modkit 2>/dev/null
  log "modkit $($X/tools/modkit --version 2>&1 | head -1)"
fi
# Rfam 5S (RF00001)
curl -sL -o $X/rfam/RF00001.cm https://rfam.org/family/RF00001/cm
log "RF00001.cm $(grep -c '^INFERNAL' $X/rfam/RF00001.cm) models, $(wc -c < $X/rfam/RF00001.cm) bytes"
# the NCBI assembly, through the datasets API
getncbi() {
  a=$1; nm=$2
  [ -s $X/asm/$nm.fna ] && { log "$nm already present"; return; }
  if [ "$a" = "GCA_049306965.1" ] && [ -s /tmp/grcz12.zip ]; then
    cp /tmp/grcz12.zip $X/asm/$nm.zip
  else
    curl -sL --retry 3 -o $X/asm/$nm.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$a/download?include_annotation_type=GENOME_FASTA"
  fi
  unzip -o -q $X/asm/$nm.zip -d $X/asm/$nm.d && cat $(find $X/asm/$nm.d -name "*.fna") > $X/asm/$nm.fna && rm -rf $X/asm/$nm.d
  samtools faidx $X/asm/$nm.fna
  log "$nm ($a) : sequences $(wc -l < $X/asm/$nm.fna.fai), $(awk '{s+=$2} END{printf "%.1f Mb", s/1e6}' $X/asm/$nm.fna.fai), md5 $(md5sum $X/asm/$nm.zip | cut -c1-12)"
}
getncbi GCA_049306965.1 zebrafish_GRCz12tu
getncbi GCA_051427915.1 zfinch_bTaeGut7_mat
getncbi GCA_051428105.1 zfinch_bTaeGut7_pat
getncbi GCA_022117705.1 maize_Mo17_T2T
getncbi GCA_028009825.2 arabidopsis_ColCC
# zebrafish GRCz11; the ENA BAM coordinate system uses the Ensembl names 1, 2, ... over 993 sequences, so the Ensembl toplevel is used
if [ ! -s $X/asm/zebrafish_GRCz11_ens.fa ]; then
  curl -sL --retry 3 -o $X/asm/zebrafish_GRCz11_ens.fa.gz https://ftp.ensembl.org/pub/release-113/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz
  gzip -dc $X/asm/zebrafish_GRCz11_ens.fa.gz > $X/asm/zebrafish_GRCz11_ens.fa && samtools faidx $X/asm/zebrafish_GRCz11_ens.fa
  log "GRCz11 Ensembl toplevel : sequences $(wc -l < $X/asm/zebrafish_GRCz11_ens.fa.fai)"
fi
log "done"
touch $X/log/SETUP_DONE
