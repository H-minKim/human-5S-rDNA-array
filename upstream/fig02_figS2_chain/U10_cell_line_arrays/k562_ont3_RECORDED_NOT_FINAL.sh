#!/bin/bash
# NOT the final version that ran. The node copy differs in download/check lines only
# (run URL from sdl_url.py; patch_ont3.py: tolerate missing md5; disk threshold).
# The line that makes the data, campaign_extract_v5.py (U1a) on each run, then concatenation into
# K562_ONT_ALL.hits_5S.fastq.gz, is not confirmed against the node copy.
# The UL2 assembly (last step here) is run by ul2_run.sh.
# K562 ONT, remaining 21 runs: collect only the 5S reads, then reassemble with --ul together with HiFi
#   purpose = more spanning reads to settle the array unit counts (assembly 139/129 vs reads ~112)
#   one run at a time (disk 47G), delete the .sra when done, piece names %02d (keeps the piece order)
#   the final assembly uses --ul (HiFi + ONT) instead of --ont (accuracy 98-99%, cannot separate haplotypes)
set -u -o pipefail
W=~/p11_encode/k562ont; mkdir -p $W/tmp $W/hits; cd $W
G5=~/p11_encode/v5; BB=$HOME/opt/bbmap/bbduk.sh
log(){ echo "[$(date +%H:%M)] $*" >> $W/RUN3.log; }
log "start : remaining run $(wc -l < ~/p11_encode/k562ont_todo.txt)"
for acc in $(cat ~/p11_encode/k562ont_todo.txt); do
  [ -e hits/$acc.hits.fastq.gz ] && continue
  u=$(curl -s "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=$acc&location=s3.us-east-1" | python3 -c "import sys,json; d=json.load(sys.stdin)['result'][0]['files'][0]; print(d['locations'][0]['link'], d['size'], d['md5'])")
  [ -z "$u" ] && { log "$acc address not received"; continue; }
  set -- $u; url=$1; sz=$2; md5=$3
  df_free=$(df --output=avail -BG ~ | tail -1 | tr -dc 0-9)
  if [ "$df_free" -lt 40 ]; then log "disk ${df_free}G, stopped"; break; fi
  log "$acc download $((sz/1000000000)) GB (free ${df_free}G)"
  N=16; P=$(( (sz + N - 1) / N ))
  for i in $(seq 0 $((N-1))); do s=$((i*P)); e=$(( s+P-1 )); [ $e -ge $sz ] && e=$((sz-1)); curl -sS --retry 5 -r $s-$e -o tmp/p_$(printf %02d $i) $url & done; wait
  cat tmp/p_* > tmp/$acc.sra && rm -f tmp/p_*
  got=$(md5sum < tmp/$acc.sra | cut -c1-32)
  if [ "$got" != "$md5" ]; then log "$acc md5 mismatch, once more"; rm -f tmp/$acc.sra
    curl -sS --retry 5 -o tmp/$acc.sra $url
    got=$(md5sum < tmp/$acc.sra | cut -c1-32)
    [ "$got" = "$md5" ] || { log "$acc md5 mismatch twice, skipped"; rm -f tmp/$acc.sra; continue; }
  fi
  fastq-dump --stdout tmp/$acc.sra 2> tmp/$acc.fqd.log \
    | python3 $G5/campaign_extract_v5.py --nts $G5/5S_unit_discseg.fa --gene $G5/5S_gene_query.fa --gene-db $G5/gene_db --bbduk $BB --out tmp/$acc --threads 10 --tmp tmp > tmp/$acc.v5.log 2>&1
  n=$(( $(wc -l < tmp/$acc.hits.fastq) / 4 ))
  gzip -c tmp/$acc.hits.fastq > hits/$acc.hits.fastq.gz; rm -f tmp/$acc.hits.fastq tmp/$acc.sra
  log "$acc 5S read $n · $(tail -1 tmp/$acc.v5.log)"
done
zcat hits/*.hits.fastq.gz | gzip -c > K562_ONT_ALL.hits_5S.fastq.gz
log "merged 5S read $(( $(zcat K562_ONT_ALL.hits_5S.fastq.gz | wc -l) / 4 )), median length $(zcat K562_ONT_ALL.hits_5S.fastq.gz | awk 'NR%4==2{print length($0)}' | sort -n | awk '{a[NR]=$1} END{print a[int(NR/2)]}'), longest $(zcat K562_ONT_ALL.hits_5S.fastq.gz | awk 'NR%4==2{if(length($0)>m)m=length($0)} END{print m}')"
U=~/p11_encode/k562ul2; mkdir -p $U; cd $U
~/hifiasm -o K562ul2 -t 14 --ul $W/K562_ONT_ALL.hits_5S.fastq.gz ~/p11_encode/k562hifi_dl/k562_region.fastq.gz > hifiasm.log 2>&1
for g in *p_ctg.gfa; do awk '/^S/{print ">"$2"\n"$3}' $g > ${g%.gfa}.fa; done
log "UL2 assembly done, $(ls $U/*p_ctg.fa 2>/dev/null | wc -l) fa"
log K562ONT3_DONE
