# -*- coding: utf-8 -*-
"""
S14 - validate the segments from S13 and split them more finely.

Problem found in S13
  k=3 of basis A is skewed at 49,097 : 1,818 : 881, so it was 'outlier picking', not segmentation.
  The length contribution in that state (intermediate type 0.97) comes from the intermediate type being 95% of the units, which is close to a tautology.

What is done here
  (1) Individual-specificity check: how many haplotypes each cluster comes from, and the largest % taken by one haplotype.
      A cluster made by a few individuals is a trait of those people, not a 'segment'.
  (2) Finer: raise to k=12 and, **for each cluster separately**, compare positional concentration with a null (shuffling the order within the array).
      Only clusters with positional structure are accepted as segments.
  (3) Recompute the length contribution with the accepted segments. Large residual clusters are grouped as 'BULK' and shown separately.
Output: S14_segments.tsv, S14_length.tsv
"""
import sys,io,os,collections
import numpy as np
from sklearn.cluster import KMeans
import v11_io as G

PRESET=sys.argv[1] if len(sys.argv)>1 else "primary"
OUT=os.path.join(G.OUT,PRESET); os.makedirs(OUT,exist_ok=True)
Z=np.load(os.path.join(OUT,"S13_labels.npy"))
LA,subA,LB,subB,ap,ad,Nn,hid=[Z[i].astype(np.int64) for i in range(8)]
rel=(ap-0.5)/Nn
NU=len(ap); NH=int(hid.max())+1
rng=np.random.default_rng(20260913)

print("=== (1) 개체 특이성 검증 — S13 의 주구획 ===")
print("  구획은 '여러 하플로타입에서 고루' 나와야 array 의 구획이다.")
print("  기준 군집  n유닛   기여 하플로타입  최대 1개가 차지한 몫  상위3 몫")
for tag,lab,names in (("A",LA,{0:"c0",1:"c1",2:"c2"}),("B",LB,{0:"c0",1:"c1",2:"c2"})):
    for c in sorted(set(lab.tolist())):
        m=lab==c
        cc=collections.Counter(hid[m].tolist())
        n=int(m.sum()); tot=float(n)
        top=cc.most_common(3)
        print("   %s %-4s %7d   %6d          %.3f               %.3f"%(
            tag,"cl%d"%c,n,len(cc),top[0][1]/tot,sum(t[1] for t in top)/tot))

# ---------- rebuild the features (same definitions as S13) ----------
U,A_,M=G.load(PRESET,verbose=False)
CAN=G.CAN; BLOCKS=G.BLOCKS
HK=sorted(U); HIi={k:i for i,k in enumerate(HK)}
units=[]
for k in HK:
    for r in U[k]: units.append((HIi[k],r[4]))
assert len(units)==NU
def parse(t):
    if t.startswith("DB"):
        b=int(t[2:].split(":")[0]); k=int(t[2:].split(":")[1])
        return ("del",(BLOCKS[b-1][0] if 1<=b<=len(BLOCKS) else None),"","",k,True)
    if t.startswith("IB"):
        b=int(t[2:].split(":")[0]); k=int(t[2:].split(":")[1])
        return ("ins",(BLOCKS[b-1][0] if 1<=b<=len(BLOCKS) else None),"","",k,True)
    if t[0]=="D":
        p,k=t[1:].split(":"); return ("del",int(p),"","",int(k),False)
    if t[0]=="I":
        q=t[1:].split(":"); return ("ins",int(q[0]),"","",int(q[1]),False)
    if t[0]=="R":
        p=int(t[1:-2]); return ("snv",p,t[-2],t[-1],1,True)
    p=int(t[:-2]); return ("snv",p,t[-2],t[-1],1,False)
def is_cpg(p,ref):
    if p is None: return False
    if ref=="C": return p<2240 and CAN[p]=="G"
    if ref=="G": return p>=2 and CAN[p-2]=="C"
    return False
def elem(p):
    if p is None: return None
    for n,a,b in (("AluY",1,318),("NTS1",319,1178),("GArich",1179,1551),("NTS2",1552,2121),("coding",2122,2240)):
        if a<=p<=b: return n
cnt=collections.Counter()
for u in units: cnt.update(u[1])
VA=sorted([v for v,c in cnt.items() if 0.005*NU<=c<=0.995*NU]); VAI={v:i for i,v in enumerate(VA)}
XA=np.zeros((NU,len(VA)),dtype=np.float32)
FB=["deamination_CpG","C>T_nonCpG","oxidation_C>A_G>T","transition_T>C_A>G","tv_C>G","tv_T>A","tv_T>G",
    "del_repeat","del_unique","del_ge10","ins_repeat","ins_unique","coding_snv_xU84C","U84C","coding_indel",
    "AluY","NTS1","GArich","NTS2","n_var"]
FBI={f:i for i,f in enumerate(FB)}
XB=np.zeros((NU,len(FB)),dtype=np.float32)
for i,(h,vs) in enumerate(units):
    for t in vs:
        j=VAI.get(t)
        if j is not None: XA[i,j]=1.0
        cl,p,r,a,ln,rep=parse(t)
        XB[i,FBI["n_var"]]+=1
        if cl=="snv":
            if (r,a) in (("C","T"),("G","A")): XB[i,FBI["deamination_CpG" if is_cpg(p,r) else "C>T_nonCpG"]]+=1
            elif (r,a) in (("C","A"),("G","T")): XB[i,FBI["oxidation_C>A_G>T"]]+=1
            elif (r,a) in (("T","C"),("A","G")): XB[i,FBI["transition_T>C_A>G"]]+=1
            elif (r,a) in (("C","G"),("G","C")): XB[i,FBI["tv_C>G"]]+=1
            elif (r,a) in (("T","A"),("A","T")): XB[i,FBI["tv_T>A"]]+=1
            elif (r,a) in (("T","G"),("A","C")): XB[i,FBI["tv_T>G"]]+=1
            if p is not None and p>=2122:
                XB[i,FBI["U84C" if p==2205 else "coding_snv_xU84C"]]+=1
        else:
            XB[i,FBI[("del_repeat" if rep else "del_unique") if cl=="del" else ("ins_repeat" if rep else "ins_unique")]]+=1
            if cl=="del" and ln>=10: XB[i,FBI["del_ge10"]]+=1
            if p is not None and p>=2122: XB[i,FBI["coding_indel"]]+=1
        e=elem(p)
        if e and e!="coding": XB[i,FBI[e]]+=1
def center_scale(X):
    s=np.zeros((NH,X.shape[1]),dtype=np.float64); c=np.zeros(NH)
    np.add.at(s,hid,X.astype(np.float64)); np.add.at(c,hid,1)
    m=s/np.maximum(c,1)[:,None]
    Xc=X-m[hid]; sd=Xc.std(0); sd[sd==0]=1
    return (Xc/sd).astype(np.float32)
ZA=center_scale(XA); ZB=center_scale(XB)

# ---------- (2) k=12 + null test per cluster ----------
NB=20
B0=np.minimum((rel*NB).astype(int),NB-1)
tot=np.bincount(B0,minlength=NB).astype(float)
order=np.argsort(hid,kind="stable"); bounds=np.searchsorted(hid[order],np.arange(NH+1))
def tvd(mask,Bv):
    h=np.bincount(Bv[mask],minlength=NB).astype(float)
    f=h/np.maximum(tot,1)
    if f.sum()==0: return np.nan
    f=f/f.sum()
    return float(np.abs(f-1.0/NB).sum()/2)
def null_perm(rep=200):
    out=[]
    for _ in range(rep):
        Bp=B0.copy()
        for h in range(NH):
            idx=order[bounds[h]:bounds[h+1]]
            if len(idx)>1: Bp[idx]=B0[rng.permutation(idx)]
        out.append(Bp)
    return out
PERM=null_perm(200)

g=io.open(os.path.join(OUT,"S14_segments.tsv"),"w",encoding="utf8",newline="\n")
g.write("basis\tcluster\tn_units\tn_haplotypes\ttop_hap_share\trel_mean\tap_median\tad_median\ttvd\ttvd_null_p975\tpositional\tlabel\ttop_features\n")
KEEP={}
for tag,Z_,names in (("A",ZA,VA),("B",ZB,FB)):
    K=12
    lab=KMeans(n_clusters=K,n_init=20,random_state=0).fit(Z_).labels_
    print("\n=== (2) 기준 %s · k=%d — 군집마다 따로 위치 검정 ==="%(tag,K))
    print("  군집   n유닛  하플로  최대몫  상대위치  근위중앙 원위중앙   TVD   귀무97.5%  판정   특징")
    rows=[]
    for c in range(K):
        m=lab==c
        n=int(m.sum())
        if n<30: continue
        cc=collections.Counter(hid[m].tolist())
        topshare=cc.most_common(1)[0][1]/float(n)
        t=tvd(m,B0)
        nulls=[tvd(m,Bp) for Bp in PERM]
        n975=float(np.nanpercentile(nulls,97.5))
        posit=(t>n975) and (len(cc)>=20) and (topshare<0.5)
        mu=Z_[m].mean(0); o=np.argsort(mu)[::-1][:4]
        feat=", ".join("%s%+.1f"%(names[j],mu[j]) for j in o)
        rows.append((c,n,len(cc),topshare,float(np.mean(rel[m])),float(np.median(ap[m])),float(np.median(ad[m])),t,n975,posit,feat))
    rows.sort(key=lambda r:r[4])
    lbl={}
    for c,n,nh,ts,rm,apm,adm,t,n975,posit,feat in rows:
        if not posit: name="BULK/개체특이"
        elif rm<0.35: name="근위구획"
        elif rm>0.65: name="원위구획"
        else: name="가운데구획"
        lbl[c]=name
        print("  %2d  %7d  %5d  %.3f   %.3f    %5.0f    %5.0f  %.3f   %.3f   %-12s %s"%(
            c,n,nh,ts,rm,apm,adm,t,n975,name+("" if posit else ""),feat[:60]))
        g.write("%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.1f\t%.1f\t%.4f\t%.4f\t%s\t%s\t%s\n"%(
            tag,c,n,nh,ts,rm,apm,adm,t,n975,"Y" if posit else "N",name,feat))
    KEEP[tag]=(lab,lbl)
g.close()

# ---------- (3) length contribution ----------
Nh=np.array([U[k][0][2] for k in HK],dtype=float)
print("\n=== (3) 길이 기여 — 채택된 구획 기준 (기울기 합 = 1) ===")
g=io.open(os.path.join(OUT,"S14_length.tsv"),"w",encoding="utf8",newline="\n")
g.write("basis\tgroup\tn_units\tmean_size\tslope\tslope_lo\tslope_hi\tintercept\tspearman\n")
for tag in ("A","B"):
    lab,lbl=KEEP[tag]
    groups=collections.defaultdict(list)
    for c,name in lbl.items(): groups[name].append(c)
    # clusters left out for having fewer than 30 units go to BULK
    print("\n  [기준 %s]"%tag)
    print("   구획          유닛수   평균크기   기울기 b [95%CI]        절편   rho")
    A2=np.vstack([np.ones_like(Nh),Nh]).T
    res=[]
    assigned=np.zeros(NU,dtype=bool)
    for name,cs in groups.items():
        m=np.isin(lab,cs); assigned|=m
        cntv=np.zeros(NH); np.add.at(cntv,hid[m],1)
        b=np.linalg.lstsq(A2,cntv,rcond=None)[0]
        bs=[np.linalg.lstsq(A2[ii],cntv[ii],rcond=None)[0][1] for ii in (rng.integers(0,NH,NH) for _ in range(400))]
        lo,hi=np.percentile(bs,[2.5,97.5])
        rx=np.argsort(np.argsort(Nh)).astype(float); ry=np.argsort(np.argsort(cntv)).astype(float)
        rx-=rx.mean(); ry-=ry.mean()
        sp=float((rx*ry).sum()/np.sqrt((rx*rx).sum()*(ry*ry).sum()))
        res.append((name,int(m.sum()),cntv.mean(),b[1],lo,hi,b[0],sp))
    if (~assigned).sum()>0:
        cntv=np.zeros(NH); np.add.at(cntv,hid[~assigned],1)
        b=np.linalg.lstsq(A2,cntv,rcond=None)[0]
        res.append(("(소군집 잔여)",int((~assigned).sum()),cntv.mean(),b[1],np.nan,np.nan,b[0],np.nan))
    res.sort(key=lambda r:-r[3])
    for name,nu_,mu,sl,lo,hi,ic,sp in res:
        print("   %-12s %7d  %7.1f   %+.3f [%+.3f,%+.3f]  %6.1f  %+.2f"%(name,nu_,mu,sl,lo,hi,ic,sp))
        g.write("%s\t%s\t%d\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n"%(tag,name,nu_,mu,sl,lo,hi,ic,sp))
g.close()
print("\n저장: %s"%OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('S14_length.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
