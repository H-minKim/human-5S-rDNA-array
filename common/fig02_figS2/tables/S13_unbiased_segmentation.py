# -*- coding: utf-8 -*-
"""
S13 - segments arrays **from unit sequence features alone** (an unbiased approach that uses no position information at all)

Two criteria
  A. individual variants: presence (0/1) vector of variants with a global unit frequency of 0.5-99.5%
  B. variant character: deamination (CpG C>T), oxidation (C>A,G>T), other transitions, by transversion kind,
                 deletions and insertions inside/outside repeat blocks, large deletions, coding region (U84C kept separate), by segment (AluY/NTS1/GA/NTS2) and so on

Mandatory preprocessing: **centering within haplotype**
   Units of the same array share the lineage background of that array (measured: between-haplotype variance share per variant, median 0.265,
   17 of 54 above 0.5). Without centering, the clusters split by 'individual' rather than by 'position'.
   After centering, each unit value = 'how much it differs from the mean unit of that array'.

Clustering: sklearn KMeans(n_init=20, random_state=0), same settings as S3
   main segments k=3 (diagnostics for k=2..6 reported alongside) -> each main segment split again with k=2,3 into subsegments

Validation: the clusters do not know position. So 'how much a cluster concentrates at particular positions'
   is compared with a null that **shuffles unit order** within each haplotype (200 rounds).

Length contribution: per haplotype, units per segment are counted and regressed on N (total unit count).
   slope = how many units that segment gains when N grows by 1 unit. The slopes sum to 1.
Output: S13_segments_<A|B>.tsv, S13_length_contribution.tsv
"""
import sys,io,os,collections
import numpy as np
from sklearn.cluster import KMeans
import v11_io as G

PRESET=sys.argv[1] if len(sys.argv)>1 else "primary"
OUT=os.path.join(G.OUT,PRESET); os.makedirs(OUT,exist_ok=True)
U,A_,M=G.load(PRESET)
CAN=G.CAN; BLOCKS=G.BLOCKS
rng=np.random.default_rng(20260913)

# ---------- unit list ----------
HK=sorted(U)
HI={k:i for i,k in enumerate(HK)}
units=[]
for k in HK:
    for r in U[k]:
        units.append((HI[k],r[0],r[1],r[2],r[4]))     # hap, ap (proximal distance), ad (distal distance), N, variant set
NU=len(units)
hid=np.array([u[0] for u in units],dtype=np.int32)
ap=np.array([u[1] for u in units],dtype=np.int32)
ad=np.array([u[2] for u in units],dtype=np.int32)
Nn=np.array([u[3] for u in units],dtype=np.int32)
rel=(ap-0.5)/Nn                                        # 0 = proximal end, 1 = distal end
print("유닛 %d · 하플로타입 %d"%(NU,len(HK)))

# ---------- variant parsing ----------
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

# ---------- feature A: individual variants ----------
cnt=collections.Counter()
for u in units: cnt.update(u[4])
VA=sorted([v for v,c in cnt.items() if 0.005*NU<=c<=0.995*NU])
VAI={v:i for i,v in enumerate(VA)}
XA=np.zeros((NU,len(VA)),dtype=np.float32)
for i,u in enumerate(units):
    for v in u[4]:
        j=VAI.get(v)
        if j is not None: XA[i,j]=1.0
print("특징 A (개별 변이) : %d 종"%len(VA))

# ---------- feature B: variant character ----------
FB=["deamination_CpG","C>T_nonCpG","oxidation_C>A_G>T","transition_T>C_A>G",
    "tv_C>G","tv_T>A","tv_T>G","del_repeat","del_unique","del_ge10",
    "ins_repeat","ins_unique","coding_snv_xU84C","U84C","coding_indel",
    "AluY","NTS1","GArich","NTS2","n_var"]
FBI={f:i for i,f in enumerate(FB)}
XB=np.zeros((NU,len(FB)),dtype=np.float32)
for i,u in enumerate(units):
    for t in u[4]:
        cl,p,r,a,ln,rep=parse(t)
        XB[i,FBI["n_var"]]+=1
        if cl=="snv":
            if (r,a) in (("C","T"),("G","A")):
                XB[i,FBI["deamination_CpG" if is_cpg(p,r) else "C>T_nonCpG"]]+=1
            elif (r,a) in (("C","A"),("G","T")): XB[i,FBI["oxidation_C>A_G>T"]]+=1
            elif (r,a) in (("T","C"),("A","G")): XB[i,FBI["transition_T>C_A>G"]]+=1
            elif (r,a) in (("C","G"),("G","C")): XB[i,FBI["tv_C>G"]]+=1
            elif (r,a) in (("T","A"),("A","T")): XB[i,FBI["tv_T>A"]]+=1
            elif (r,a) in (("T","G"),("A","C")): XB[i,FBI["tv_T>G"]]+=1
            if p is not None and p>=2122:
                if p==2205: XB[i,FBI["U84C"]]+=1
                else: XB[i,FBI["coding_snv_xU84C"]]+=1
        else:
            XB[i,FBI[("del_repeat" if rep else "del_unique") if cl=="del" else ("ins_repeat" if rep else "ins_unique")]]+=1
            if cl=="del" and ln>=10: XB[i,FBI["del_ge10"]]+=1
            if p is not None and p>=2122: XB[i,FBI["coding_indel"]]+=1
        e=elem(p)
        if e and e!="coding": XB[i,FBI[e]]+=1
print("특징 B (변이 성격) : %d 개"%len(FB))

# ---------- centering within haplotype plus standardization ----------
def center_scale(X):
    Xc=X.copy()
    s=np.zeros((len(HK),X.shape[1]),dtype=np.float64); c=np.zeros(len(HK))
    np.add.at(s,hid,X.astype(np.float64)); np.add.at(c,hid,1)
    m=s/np.maximum(c,1)[:,None]
    Xc=X-m[hid]
    sd=Xc.std(0); sd[sd==0]=1
    return (Xc/sd).astype(np.float32)
ZA=center_scale(XA); ZB=center_scale(XB)

# ---------- position profile ----------
def profile(lab,K):
    """Where each cluster sits (position was not used in clustering)"""
    out={}
    for c in range(K):
        m=lab==c
        if m.sum()==0: continue
        out[c]=dict(n=int(m.sum()),
                    ap=float(np.median(ap[m])),ad=float(np.median(ad[m])),
                    rel=float(np.mean(rel[m])),
                    f_prox8=float((ap[m]<=8).mean()),f_dist8=float((ad[m]<=8).mean()),
                    f_mid=float(((ap[m]>30)&(ad[m]>30)).mean()))
    return out
def concentration(lab,K):
    """Position concentration = maximum deviation of the distribution over 20 relative-position bins (mean per cluster)"""
    B=np.minimum((rel*20).astype(int),19)
    tot=np.bincount(B,minlength=20).astype(float)
    v=[]
    for c in range(K):
        m=lab==c
        if m.sum()<50: continue
        h=np.bincount(B[m],minlength=20).astype(float)
        f=h/np.maximum(tot,1); f=f/f.sum()
        v.append(float(np.abs(f-1.0/20).sum()/2))       # total variation distance
    return float(np.mean(v)) if v else float("nan")
def null_conc(lab,K,rep=200):
    """Null that shuffles unit order within each haplotype (features unchanged, only positions reassigned)"""
    order=np.argsort(hid,kind="stable")
    bounds=np.searchsorted(hid[order],np.arange(len(HK)+1))
    out=[]
    B0=np.minimum((rel*20).astype(int),19)
    tot=np.bincount(B0,minlength=20).astype(float)
    for _ in range(rep):
        Bp=B0.copy()
        for h in range(len(HK)):
            idx=order[bounds[h]:bounds[h+1]]
            if len(idx)>1: Bp[idx]=B0[rng.permutation(idx)]
        v=[]
        for c in range(K):
            m=lab==c
            if m.sum()<50: continue
            hh=np.bincount(Bp[m],minlength=20).astype(float)
            f=hh/np.maximum(tot,1); f=f/f.sum()
            v.append(float(np.abs(f-1.0/20).sum()/2))
        out.append(np.mean(v) if v else np.nan)
    return float(np.nanmean(out)),float(np.nanpercentile(out,97.5))

def run(Z,tag,name):
    print("\n"+"#"*70)
    print("### 기준 %s — %s"%(tag,name))
    print("#"*70)
    print("\n[진단] k 별 위치 집중도 (군집은 위치를 모른다 · 귀무는 배열 안 순서 섞기)")
    diag={}
    for K in (2,3,4,5,6,8):
        lab=KMeans(n_clusters=K,n_init=20,random_state=0).fit(Z).labels_
        cc=concentration(lab,K)
        diag[K]=(lab,cc)
        print("   k=%d  집중도 %.4f"%(K,cc))
    KMAIN=3
    lab=diag[KMAIN][0]
    nm,nhi=null_conc(lab,KMAIN)
    print("   주구획 k=%d : 집중도 %.4f · 귀무 평균 %.4f (97.5%% %.4f) → %s"%(
        KMAIN,diag[KMAIN][1],nm,nhi,"위치 구조 있음" if diag[KMAIN][1]>nhi else "귀무와 구별 안 됨"))
    pr=profile(lab,KMAIN)
    # name clusters by median relative position (0 = proximal)
    order=sorted(pr,key=lambda c:pr[c]["rel"])
    nameof={}
    for rank,c in enumerate(order):
        nameof[c]=["근위형","중간형","원위형"][rank] if KMAIN==3 else "군집%d"%c
    print("\n[주구획] 위치는 군집화에 쓰이지 않았다 — 아래는 결과로 나온 분포")
    print("   구획      n       상대위치평균  근위거리중앙 원위거리중앙  근위8내 원위8내 가운데")
    for c in order:
        p=pr[c]
        print("   %-8s %6d   %.3f        %5.0f       %5.0f      %.3f  %.3f  %.3f"%(
            nameof[c],p["n"],p["rel"],p["ap"],p["ad"],p["f_prox8"],p["f_dist8"],p["f_mid"]))
    # subsegments
    sub=np.full(NU,-1,dtype=np.int32); subname={}
    nxt=0
    print("\n[하위구획] 각 주구획을 다시 k=3 으로")
    for c in order:
        m=lab==c
        sl=KMeans(n_clusters=3,n_init=20,random_state=0).fit(Z[m]).labels_
        idx=np.where(m)[0]
        rr=[(np.mean(rel[idx[sl==s]]),s) for s in range(3) if (sl==s).sum()>0]
        rr.sort()
        for rank,(rv,s) in enumerate(rr):
            gi=nxt; nxt+=1
            sub[idx[sl==s]]=gi
            subname[gi]="%s-%d"%(nameof[c],rank+1)
            mm=idx[sl==s]
            print("   %-10s n=%6d  상대위치 %.3f  근위거리중앙 %5.0f  원위거리중앙 %5.0f  가운데비율 %.3f"%(
                subname[gi],len(mm),rv,np.median(ap[mm]),np.median(ad[mm]),
                float(((ap[mm]>30)&(ad[mm]>30)).mean())))
    # save
    g=io.open(os.path.join(OUT,"S13_segments_%s.tsv"%tag),"w",encoding="utf8",newline="\n")
    g.write("level\tsegment\tn_units\trel_mean\tap_median\tad_median\tfrac_prox8\tfrac_dist8\tfrac_mid\n")
    for c in order:
        p=pr[c]
        g.write("main\t%s\t%d\t%.4f\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\n"%(nameof[c],p["n"],p["rel"],p["ap"],p["ad"],p["f_prox8"],p["f_dist8"],p["f_mid"]))
    for gi in sorted(subname):
        m=sub==gi
        g.write("sub\t%s\t%d\t%.4f\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\n"%(subname[gi],int(m.sum()),float(np.mean(rel[m])),
                float(np.median(ap[m])),float(np.median(ad[m])),float((ap[m]<=8).mean()),float((ad[m]<=8).mean()),
                float(((ap[m]>30)&(ad[m]>30)).mean())))
    g.close()
    # feature interpretation: the most prominent features in each main segment
    print("\n[무엇이 그 구획을 만드나] 중심화한 값의 군집 평균 (표준편차 단위) 상위/하위 5")
    for c in order:
        m=lab==c
        mu=Z[m].mean(0)
        names=(VA if tag=="A" else FB)
        o=np.argsort(mu)
        top=[(names[j],mu[j]) for j in o[::-1][:5]]
        bot=[(names[j],mu[j]) for j in o[:5]]
        print("   %-8s ↑ %s"%(nameof[c],", ".join("%s %+.2f"%(n,v) for n,v in top)))
        print("            ↓ %s"%(", ".join("%s %+.2f"%(n,v) for n,v in bot)))
    return lab,order,nameof,sub,subname

LA,orderA,nameA,subA,subnameA=run(ZA,"A","개별 변이 (SNV+indel) 종합")
LB,orderB,nameB,subB,subnameB=run(ZB,"B","변이의 성격 (탈아미노화·산화·코딩부 등)")

# ---------- length contribution ----------
print("\n"+"#"*70)
print("### 하플로타입 간 유닛 수 차이는 '어느 구획이 길어져서' 인가")
print("#"*70)
Nh=np.array([U[k][0][2] for k in HK],dtype=float)
g=io.open(os.path.join(OUT,"S13_length_contribution.tsv"),"w",encoding="utf8",newline="\n")
g.write("basis\tlevel\tsegment\tslope\tslope_lo\tslope_hi\tintercept\tspearman\tmean_size\n")
def contrib(lab,names,basis,level):
    K=int(lab.max())+1
    print("\n  [%s · %s] 구획 크기 = a + b*N  (b 합 = 1)"%(basis,level))
    print("   구획        평균 크기   기울기 b [95%CI]        해석")
    tot=np.zeros(len(HK))
    rows=[]
    for c in range(K):
        cnt=np.zeros(len(HK))
        np.add.at(cnt,hid[lab==c],1)
        tot+=cnt
        A=np.vstack([np.ones_like(Nh),Nh]).T
        b=np.linalg.lstsq(A,cnt,rcond=None)[0]
        bs=[]
        for _ in range(400):
            ii=rng.integers(0,len(HK),len(HK))
            bs.append(np.linalg.lstsq(A[ii],cnt[ii],rcond=None)[0][1])
        lo,hi=np.percentile(bs,[2.5,97.5])
        rx=np.argsort(np.argsort(Nh)).astype(float); ry=np.argsort(np.argsort(cnt)).astype(float)
        rx-=rx.mean(); ry-=ry.mean()
        sp=float((rx*ry).sum()/np.sqrt((rx*rx).sum()*(ry*ry).sum()))
        rows.append((names[c] if c in names else str(c),cnt.mean(),b[1],lo,hi,b[0],sp))
    rows.sort(key=lambda r:-r[2])
    for nm,mu,sl,lo,hi,ic,sp in rows:
        tag=("길이 변이의 주된 자리" if sl>0.4 else ("거의 고정" if abs(sl)<0.1 else ""))
        print("   %-10s %8.1f   %+.3f [%+.3f,%+.3f]   절편 %6.1f · rho %+.2f  %s"%(nm,mu,sl,lo,hi,ic,sp,tag))
        g.write("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.3f\n"%(basis,level,nm,sl,lo,hi,ic,sp,mu))
contrib(LA,nameA,"A","주구획"); contrib(subA,subnameA,"A","하위구획")
contrib(LB,nameB,"B","주구획"); contrib(subB,subnameB,"B","하위구획")
g.close()
np.save(os.path.join(OUT,"S13_labels.npy"),np.vstack([LA,subA,LB,subB,ap,ad,Nn,hid]))
print("\n저장: %s"%OUT)
