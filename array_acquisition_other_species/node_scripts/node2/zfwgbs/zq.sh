#!/bin/bash
set -u
cd ~/zfwgbs
python3 - <<'PY'
import pysam,io,collections,bisect
TAB=chr(9)
parts=collections.defaultdict(list)
for L in io.open("zf_regions_parts.tsv",encoding="utf8"):
    f=L.rstrip(chr(10)).split(TAB)
    b,_,k=f[0].rpartition("__")
    parts[(b,k)].append((f[1],int(f[2]),int(f[3])))
print("sample\tpart\tn_read\tmapq0%\tmapq>=20%\tmean_mapq\tCpG_per_kb")
for bam,lab in [("bam/SRR800059_deep.bam","알_Potok"),("bam/SRR800081_deep.bam","근육_Potok"),("bam/egg_deep.bam","알_Jiang")]:
    b=pysam.AlignmentFile(bam,"rb")
    for k in ("코딩부","스페이서"):
        iv=sorted((s,e) for _,s,e in parts[("chr18_somatic",k)])
        ch=parts[("chr18_somatic",k)][0][0]
        st=[s for s,e in iv]; en=[e for s,e in iv]
        bp=sum(e-s+1 for s,e in iv)
        def ins(p):
            i=bisect.bisect_right(st,p)-1
            return i>=0 and p<=en[i]
        qs=[]; nc=0
        for r in b.fetch(ch,iv[0][0]-1,iv[-1][1]):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate: continue
            try: xm=r.get_tag("XM")
            except KeyError: continue
            c=0
            for qp,rp in r.get_aligned_pairs(matches_only=True):
                if rp is None or qp is None or not ins(rp+1): continue
                if (xm[qp] if qp<len(xm) else ".") in ("z","Z"): c+=1
            if c: qs.append(r.mapping_quality); nc+=c
        if qs:
            print("%s\t%s\t%d\t%.1f\t%.1f\t%.1f\t%.2f"%(lab,k,len(qs),100.0*sum(1 for q in qs if q==0)/len(qs),100.0*sum(1 for q in qs if q>=20)/len(qs),sum(qs)/float(len(qs)),1000.0*nc/bp))
    b.close()
PY
