#!/bin/bash
set -u
cd ~/zfwgbs
ls ref/*.fa ref/*.fasta 2>/dev/null | head
python3 - <<'PY'
import io,collections,glob,pysam
TAB=chr(9)
f=[p for p in glob.glob("ref/*.fa")+glob.glob("ref/*.fasta") if "bisulfite" not in p]
print("참조:",f[:2])
fa=pysam.FastaFile(f[0])
parts=collections.defaultdict(list)
for L in io.open("zf_regions_parts.tsv",encoding="utf8"):
    x=L.rstrip(chr(10)).split(TAB)
    b,_,k=x[0].rpartition("__")
    parts[(b,k)].append((x[1],int(x[2]),int(x[3])))
print("array\t부류\tCpG\tbp\tCpG/kb")
for base in ["chr18_somatic","chr4_maternal_A"]:
    for k in ("코딩부","스페이서","바깥_flankL"):
        iv=parts.get((base,k))
        if not iv: continue
        n=bp=0
        for ch,s,e in iv:
            sq=fa.fetch(ch,s-1,e).upper()
            n+=sq.count("CG"); bp+=len(sq)
        print("%s\t%s\t%d\t%d\t%.2f"%(base,k,n,bp,1000.0*n/bp))
PY
