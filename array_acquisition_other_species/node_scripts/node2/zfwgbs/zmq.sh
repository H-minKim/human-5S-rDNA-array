#!/bin/bash
set -u
cd ~/zfwgbs
sed 's/def boot(per_read):/MQ=int(__import__("os").environ.get("MQ","0"))\ndef boot(per_read):/; s/if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:/if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.mapping_quality < MQ:/' bismark_parts_meth.py > bismark_parts_mq.py
for M in 0 30; do
  export MQ=$M
  rm -f ZFPARTS_mq${M}.tsv
  for X in "egg_deep 알_Jiang" "SRR800059_deep 알_Potok" "SRR800081_deep 근육성체_Potok" \
           "SRR800077 sphere_4hpf" "SRR837895 Germring_6hpf" "SRR12797061 장성체" "SRR837894 1k세포_MBT"; do
    set -- $X
    [ -f bam/$1.bam ] && B=bam/$1.bam || B=bam/$1.sorted.bam
    [ -f "$B" ] && python3 bismark_parts_mq.py "$B" "$2" zf_regions_parts.tsv ZFPARTS_mq${M}.tsv 1000 >/dev/null
  done
  echo "MQ>=$M 완료"
done
# read/kb 몰림도 같이
python3 - <<'PY'
import io,collections
TAB=chr(9)
bp={}
for M in ("0","30"):
    d={}
    for L in io.open("ZFPARTS_mq%s.tsv"%M,encoding="utf8"):
        f=L.rstrip(chr(10)).split(TAB)
        if f[0]=="sample": continue
        d[(f[0],f[1],f[2])]=f[3:]
    print("=== mapq >= %s"%M)
    print("시료\t부류\tread\tread/kb\tCpG\tfrac\tlo95\thi95")
    for s in ["알_Jiang","알_Potok","1k세포_MBT","sphere_4hpf","Germring_6hpf","근육성체_Potok","장성체"]:
        for k in ("코딩부","스페이서","바깥_flankL","바깥_flankR"):
            v=d.get((s,"chr18_somatic",k))
            if not v or v[0]=="-": continue
            b=float(v[6]); r=int(v[0])
            print("%s\t%s\t%d\t%.2f\t%s\t%s\t%s\t%s"%(s,k,r,1000.0*r/b,v[1],v[3],v[4],v[5]))
PY
