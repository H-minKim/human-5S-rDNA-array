# -*- coding: utf-8 -*-
"""
S16 - test by measurement the claim that "the ends must not be drawn in relative position", and redraw the figure in absolute coordinates.

Claim (user)
  In section 13 the place that absorbs length variation came out as the 'middle'. Then even if a terminal segment has the same shape in absolute distance,
  a longer array **pushes it toward the end and narrows it in relative position**. So the ends must use absolute distance.

Test method (both with the same cluster labels)
  (1) For each length stratum (N40-80 / N81-110 / N111-140 / N141+) build the positional profile of each segment,
      and measure the spread between strata (standard deviation) in **absolute coordinates** (units from the junction) and in **relative coordinates** (0-1).
      The coordinate with the smaller spread is the 'conserved axis'.
  (2) How the mean position of a segment moves with N: is the absolute one (ad) fixed, and is the relative one (rel) pushed?
  (3) Prediction check: if the end is fixed in absolute terms then rel ≈ 1 − ad/N, so it should approach 1 as N grows.
Output: S16_abs_vs_rel.tsv, fig_v11_segments_abs_<preset>.png
"""
import sys,io,os,collections
import numpy as np
from sklearn.cluster import KMeans
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import v11_io as G

PRESET=sys.argv[1] if len(sys.argv)>1 else "primary"
OUT=os.path.join(G.OUT,PRESET)
U,A_,M=G.load(PRESET,verbose=False)
CAN=G.CAN; BLOCKS=G.BLOCKS
HK=sorted(U); HI={k:i for i,k in enumerate(HK)}; NH=len(HK)
units=[];ap=[];ad=[];Nn=[];hid=[]
for k in HK:
    for r in U[k]:
        units.append(r[4]); ap.append(r[0]); ad.append(r[1]); Nn.append(r[2]); hid.append(HI[k])
ap=np.array(ap);ad=np.array(ad);Nn=np.array(Nn,float);hid=np.array(hid)
NU=len(ap); rel=(ap-0.5)/Nn
rng=np.random.default_rng(20260913)

# ---------- same features and clusters as S15 ----------
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
    if t[0]=="R": return ("snv",int(t[1:-2]),t[-2],t[-1],1,True)
    return ("snv",int(t[:-2]),t[-2],t[-1],1,False)
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
for v in units: cnt.update(v)
VA=sorted([v for v,c in cnt.items() if 0.005*NU<=c<=0.995*NU]); VAI={v:i for i,v in enumerate(VA)}
FB=["deamination_CpG","C>T_nonCpG","oxidation_C>A_G>T","transition_T>C_A>G","tv_C>G","tv_T>A","tv_T>G",
    "del_repeat","del_unique","del_ge10","ins_repeat","ins_unique","coding_snv_xU84C","U84C","coding_indel",
    "AluY","NTS1","GArich","NTS2","n_var"]
FBI={f:i for i,f in enumerate(FB)}
XA=np.zeros((NU,len(VA)),dtype=np.float32); XB=np.zeros((NU,len(FB)),dtype=np.float32)
for i,vs in enumerate(units):
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
            if p is not None and p>=2122: XB[i,FBI["U84C" if p==2205 else "coding_snv_xU84C"]]+=1
        else:
            XB[i,FBI[("del_repeat" if rep else "del_unique") if cl=="del" else ("ins_repeat" if rep else "ins_unique")]]+=1
            if cl=="del" and ln>=10: XB[i,FBI["del_ge10"]]+=1
            if p is not None and p>=2122: XB[i,FBI["coding_indel"]]+=1
        e=elem(p)
        if e and e!="coding": XB[i,FBI[e]]+=1
def cs(X):
    s=np.zeros((NH,X.shape[1]),dtype=np.float64); c=np.zeros(NH)
    np.add.at(s,hid,X.astype(np.float64)); np.add.at(c,hid,1)
    Xc=X-(s/np.maximum(c,1)[:,None])[hid]; sd=Xc.std(0); sd[sd==0]=1
    return (Xc/sd).astype(np.float32)
LAB={"A":KMeans(n_clusters=12,n_init=20,random_state=0).fit(cs(XA)).labels_,
     "B":KMeans(n_clusters=12,n_init=20,random_state=0).fit(cs(XB)).labels_}
FEAT={"A":(cs(XA),VA),"B":(cs(XB),FB)}

STR=[(40,80,"N40-80"),(81,110,"N81-110"),(111,140,"N111-140"),(141,10**9,"N141+")]
g=io.open(os.path.join(OUT,"S16_abs_vs_rel.tsv"),"w",encoding="utf8",newline="\n")
g.write("basis\tcluster\tstratum\tn\tmean_ad\tmean_ap\tmean_rel\tpred_rel_from_ad\n")

print("=== (1) 좌표별 층간 흩어짐 — 작은 쪽이 '보존되는 축' ===")
print("   구획을 길이층으로 나눠 프로파일을 만들고 층 사이 표준편차를 평균한다.")
print("   (말단 구획만 대상 : 상대위치 평균이 0.35 미만 또는 0.65 초과)")
print("")
print("  기준 구획  상대위치  n유닛 | 절대좌표 흩어짐  상대좌표 흩어짐   어느 쪽이 보존?")
summ=[]
for tag in ("A","B"):
    lab=LAB[tag]
    for c in range(12):
        m=lab==c
        if m.sum()<200: continue
        rm=float(np.mean(rel[m]))
        if 0.35<=rm<=0.65: continue                      # terminal segments only
        side="dist" if rm>0.5 else "prox"
        dv=ad if side=="dist" else ap
        # absolute: at depth 1..30, the fraction of units at that depth that belong to this segment
        prof_abs=[];prof_rel=[]
        for lo,hi,lab_s in STR:
            sm=(Nn>=lo)&(Nn<=hi)
            if sm.sum()<500: prof_abs.append(None); prof_rel.append(None); continue
            pa=[]
            for d in range(1,31):
                den=((dv==d)&sm).sum()
                pa.append(((dv==d)&sm&m).sum()/den if den>=20 else np.nan)
            prof_abs.append(np.array(pa))
            pr=[]
            for b in range(20):
                bl,bh=b/20.0,(b+1)/20.0
                den=((rel>=bl)&(rel<bh)&sm).sum()
                pr.append(((rel>=bl)&(rel<bh)&sm&m).sum()/den if den>=20 else np.nan)
            prof_rel.append(np.array(pr))
        PA=np.array([p for p in prof_abs if p is not None])
        PR=np.array([p for p in prof_rel if p is not None])
        if len(PA)<3: continue
        # the spread is expressed 'relative to the mean' to remove the scale difference between the coordinates
        sa=np.nanmean(np.nanstd(PA,0))/max(np.nanmean(PA),1e-9)
        sr=np.nanmean(np.nanstd(PR,0))/max(np.nanmean(PR),1e-9)
        win="절대" if sa<sr else "상대"
        print("  %s c%-2d  %.3f  %6d |    %.3f          %.3f        %s"%(tag,c,rm,int(m.sum()),sa,sr,win))
        summ.append((tag,c,rm,sa,sr))
        for lo,hi,lab_s in STR:
            sm=(Nn>=lo)&(Nn<=hi)&m
            if sm.sum()<20: continue
            g.write("%s\t%d\t%s\t%d\t%.3f\t%.3f\t%.4f\t%.4f\n"%(tag,c,lab_s,int(sm.sum()),
                float(np.mean(ad[sm])),float(np.mean(ap[sm])),float(np.mean(rel[sm])),
                float(np.mean(1.0-ad[sm]/Nn[sm])) if rm>0.5 else float(np.mean(ap[sm]/Nn[sm]))))
g.close()
nabs=sum(1 for t,c,rm,sa,sr in summ if sa<sr)
print("")
print("  ⇒ 말단 구획 %d 개 중 %d 개에서 **절대 좌표**의 층간 흩어짐이 더 작다"%(len(summ),nabs))

print("")
print("=== (2) 구획의 평균 위치가 배열 길이에 따라 움직이는가 ===")
print("   절대(접합에서 유닛 수)가 고정이고 상대가 밀리면, 상대좌표는 길이에 오염된 축이다.")
for tag in ("A","B"):
    lab=LAB[tag]
    for c in range(12):
        m=lab==c
        rm=float(np.mean(rel[m])) if m.sum() else 9
        if m.sum()<200 or rm<=0.65: continue             # distal terminal type only
        print("")
        print("  [기준 %s · c%d · 원위 말단형 n=%d]"%(tag,c,int(m.sum())))
        print("    길이층        n    평균 원위거리(절대)   평균 상대위치   예측 1-ad/N")
        for lo,hi,lab_s in STR:
            sm=m&(Nn>=lo)&(Nn<=hi)
            if sm.sum()<20: continue
            print("    %-10s %5d       %6.2f            %.3f          %.3f"%(
                lab_s,int(sm.sum()),float(np.mean(ad[sm])),float(np.mean(rel[sm])),
                float(np.mean(1.0-ad[sm]/Nn[sm]))))

# ---------- (3) absolute-coordinate figure ----------
REG=[("P%d"%d,("ap",d)) for d in range(1,31)]+[("MID",None)]+[("D%d"%d,("ad",d)) for d in range(30,0,-1)]
def col_mask(spec):
    if spec is None: return (ap>30)&(ad>30)
    which,d=spec
    return (ap==d) if which=="ap" else (ad==d)
fig=plt.figure(figsize=(18,10)); gs=fig.add_gridspec(2,1,hspace=0.42)
CM=plt.get_cmap("turbo")
for gi,tag,title in ((0,"A","A  Segments from INDIVIDUAL VARIANTS — ABSOLUTE distance axis"),
                     (1,"B","B  Segments from VARIANT CHARACTER — ABSOLUTE distance axis")):
    ax=fig.add_subplot(gs[gi])
    lab=LAB[tag]
    Z_,names=FEAT[tag]
    cl=[]
    for c in range(12):
        m=lab==c
        if m.sum()<30: continue
        mu=Z_[m].mean(0); j=int(np.argmax(mu))
        cl.append((float(np.mean(rel[m])),c,names[j]))
    cl.sort()
    M=np.zeros((len(cl),len(REG)))
    for xi,(nm,spec) in enumerate(REG):
        mm=col_mask(spec); den=mm.sum()
        if den==0: continue
        for ci,(rm,c,f) in enumerate(cl):
            M[ci,xi]=((lab==c)&mm).sum()/den
    bottom=np.zeros(len(REG))
    xs=np.arange(len(REG))
    for ci,(rm,c,f) in enumerate(cl):
        ax.fill_between(xs,bottom,bottom+M[ci],color=CM(ci/float(max(len(cl)-1,1))),
                        alpha=0.92,lw=0,label="c%d  %s"%(c,f[:20]))
        bottom+=M[ci]
    ax.set_xlim(0,len(REG)-1); ax.set_ylim(0,1)
    ax.axvline(30,color="w",lw=1.5); ax.axvline(30,color="k",lw=0.6,ls=":")
    tk=[i for i,(nm,_) in enumerate(REG) if nm.startswith("P") and int(nm[1:])%3==1]+[30]+\
       [i for i,(nm,_) in enumerate(REG) if nm.startswith("D") and int(nm[1:])%3==1]
    ax.set_xticks(tk); ax.set_xticklabels([REG[i][0] for i in tk],fontsize=7.5,rotation=90)
    ax.set_xlabel("absolute distance from the junction  (P = units from proximal PJ · MID = both >30 · D = units from distal DJ)",fontsize=9)
    ax.set_ylabel("share of units at that distance",fontsize=9)
    ax.set_title(title+"\nno relative coordinate: terminal domains are fixed in ABSOLUTE units, so relative position would smear them",loc="left",fontsize=10)
    ax.legend(fontsize=6.5,loc="center left",bbox_to_anchor=(1.005,0.5),framealpha=0.95)
fig.suptitle("5S rDNA — sequence-only segmentation drawn on the ABSOLUTE junction-distance axis (middle kept in place)\n"
             "v11 · %s · KMeans(k=12, n_init=20, random_state=0) · position never used for clustering"%PRESET,fontsize=12,y=0.98)
p=os.path.join(G.OUT,"fig_v11_segments_abs_%s.png"%PRESET)
fig.savefig(p,dpi=150,bbox_inches="tight"); print("\nsaved:",p)
