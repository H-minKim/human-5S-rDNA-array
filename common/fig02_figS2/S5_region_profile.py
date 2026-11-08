# -*- coding: utf-8 -*-
# S5 - variant pattern and its character by region (proximal 1..30, middle, distal 30..1)
#   metric: count per unit (two-stage aggregation: mean within haplotype -> mean across haplotypes)
#   character: class (SNV/deletion/insertion x inside/outside repeat block), 12 SNV substitution types, CpG context, 5 element compartments within the unit
#   ratio = that region / middle (ratio of means), haplotype bootstrap CI
import io,os,sys,collections
import numpy as np
import v11_io as G
PRESET=sys.argv[1] if len(sys.argv)>1 else "primary"
OUT=os.path.join(G.OUT,PRESET); os.makedirs(OUT,exist_ok=True)
U,A,M=G.load(PRESET)
CAN=G.CAN
def parse(x):
    if x.startswith("DB"): return ("del",None,"","",int(x[2:].split(":")[1]),True)
    if x.startswith("IB"): return ("ins",None,"","",int(x[2:].split(":")[1]),True)
    if x[0]=="D":
        p,n=x[1:].split(":"); p=int(p); return ("del",p,CAN[p-1],"",int(n),False)
    if x[0]=="I":
        q=x[1:].split(":"); return ("ins",int(q[0]),"","",int(q[1]),False)
    if x[0]=="R":
        p=int(x[1:-2]); return ("snv",p,x[-2],x[-1],1,True)
    p=int(x[:-2]); return ("snv",p,x[-2],x[-1],1,False)
def elem(p):
    if p is None: return None
    for a,b,n in [(1,318,"AluY"),(319,1178,"NTS1"),(1179,1551,"GArich"),(1552,2121,"NTS2"),(2122,2240,"coding")]:
        if a<=p<=b: return n
def iscpg(p):
    if p is None: return False
    c=CAN[p-1]
    return (c=="C" and p<2240 and CAN[p]=="G") or (c=="G" and p>=2 and CAN[p-2]=="C")
SUBS=["C>T","G>A","C>A","G>T","C>G","G>C","T>C","A>G","T>G","A>C","T>A","A>T"]
CLS=["snv","snv_rep","del","del_rep","ins","ins_rep"]
MET=["n_var"]+CLS+["del_ge10","ins_ge10"]+["sub:"+s for s in SUBS]+["C>T|CpG","C>T|non"]+["elem:"+e for e in ("AluY","NTS1","GArich","NTS2","coding")]
MI={m:i for i,m in enumerate(MET)}
def counts(vs):
    v=np.zeros(len(MET))
    v[MI["n_var"]]=len(vs)
    for x in vs:
        cl,p,r,a,ln,rep=parse(x)
        v[MI[cl+("_rep" if rep else "")]]+=1
        if cl=="del" and ln>=10: v[MI["del_ge10"]]+=1
        if cl=="ins" and ln>=10: v[MI["ins_ge10"]]+=1
        if cl=="snv" and r in "ACGT" and a in "ACGT":
            k="sub:%s>%s"%(r,a)
            if k in MI: v[MI[k]]+=1
            if r=="C" and a=="T": v[MI["C>T|CpG" if iscpg(p) else "C>T|non"]]+=1
        e=elem(p)
        if e: v[MI["elem:"+e]]+=1
    return v
REG=["P%d"%d for d in range(1,31)]+["MID"]+["D%d"%d for d in range(30,0,-1)]
RI={r:i for i,r in enumerate(REG)}
HK=sorted(U)
S=np.zeros((len(HK),len(REG),len(MET))); C=np.zeros((len(HK),len(REG)))
for hi,k in enumerate(HK):
    for r in U[k]:
        v=counts(r[4])
        if r[0]<=30: S[hi,RI["P%d"%r[0]]]+=v; C[hi,RI["P%d"%r[0]]]+=1
        if r[1]<=30: S[hi,RI["D%d"%r[1]]]+=v; C[hi,RI["D%d"%r[1]]]+=1
        if r[0]>30 and r[1]>30: S[hi,RI["MID"]]+=v; C[hi,RI["MID"]]+=1
Cn=np.where(C>0,C,np.nan)
Rate=S/Cn[:,:,None]                      # mean within haplotype (stage 1)
def agg(idx=None):
    X=Rate if idx is None else Rate[idx]
    n=np.sum(~np.isnan(X[:,:,0]),0)
    return np.nanmean(X,0),n
mean,nhap=agg()
rng=np.random.default_rng(20260912); B=600
BS=np.empty((B,len(REG),len(MET)))
for b in range(B):
    BS[b]=agg(rng.integers(0,len(HK),len(HK)))[0]
lo=np.nanpercentile(BS,2.5,0); hi=np.nanpercentile(BS,97.5,0)
mid=mean[RI["MID"]]; midB=BS[:,RI["MID"],:]
ratio=mean/np.where(mid>0,mid,np.nan)
rBS=BS/np.where(midB>0,midB,np.nan)[:,None,:]
rlo=np.nanpercentile(rBS,2.5,0); rhi=np.nanpercentile(rBS,97.5,0)
g=io.open(os.path.join(OUT,"S5_region_profile.tsv"),"w",encoding="utf8",newline="\n")
g.write("region\tside\tunit_from_junction\tmetric\tn_haplotypes\tper_unit_mean\tci_lo\tci_hi\tratio_vs_mid\tratio_lo\tratio_hi\n")
for r in REG:
    i=RI[r]; sd=("MID" if r=="MID" else ("prox" if r[0]=="P" else "dist")); d=(-1 if r=="MID" else int(r[1:]))
    for m in MET:
        j=MI[m]
        g.write("%s\t%s\t%d\t%s\t%d\t%.5f\t%.5f\t%.5f\t%s\t%s\t%s\n"%(r,sd,d,m,nhap[i],mean[i,j],lo[i,j],hi[i,j],
            ("%.4f"%ratio[i,j]) if np.isfinite(ratio[i,j]) else "",("%.4f"%rlo[i,j]) if np.isfinite(rlo[i,j]) else "",
            ("%.4f"%rhi[i,j]) if np.isfinite(rhi[i,j]) else ""))
g.close()
print("저장: S5_region_profile.tsv  (영역 %d x 지표 %d)"%(len(REG),len(MET)))
def line(m,regs):
    j=MI[m]; return "  ".join("%5.2f"%mean[RI[r],j] for r in regs)
sel=["P1","P2","P3","P5","P10","P20","P30","MID","D30","D20","D10","D5","D3","D2","D1"]
print("")
print("=== 유닛당 개수 (영역 %s) ==="%" ".join(sel))
print("  지표          "+"  ".join("%5s"%s for s in sel))
for m in ["n_var","snv","snv_rep","del","del_rep","ins","ins_rep","del_ge10"]:
    print("  %-12s "%m+line(m,sel))
print("")
print("=== 가운데 대비 비 (CI 가 1 을 벗어나면 ) ===")
print("  지표          "+"  ".join("%5s"%s for s in sel))
for m in ["snv","snv_rep","del","del_rep","ins","ins_rep","del_ge10"]+["sub:"+s for s in SUBS]+["C>T|CpG","C>T|non"]+["elem:coding","elem:GArich","elem:AluY","elem:NTS1","elem:NTS2"]:
    j=MI[m]
    cells=[]
    for r in sel:
        i=RI[r]
        v=ratio[i,j]; st="" if (np.isfinite(rlo[i,j]) and (rlo[i,j]>1 or rhi[i,j]<1)) else " "
        cells.append(("%5.2f%s"%(v,st)) if np.isfinite(v) else "    - ")
    print("  %-12s "%m+" ".join(cells))
# ---- figure
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows=["snv","snv_rep","del","del_rep","ins","ins_rep","del_ge10"]+["sub:"+s for s in SUBS]+["C>T|CpG","C>T|non"]+["elem:AluY","elem:NTS1","elem:GArich","elem:NTS2","elem:coding"]
Mt=np.array([[np.log2(ratio[RI[r],MI[m]]) if np.isfinite(ratio[RI[r],MI[m]]) and ratio[RI[r],MI[m]]>0 else np.nan for r in REG] for m in rows])
fig=plt.figure(figsize=(17,10)); gs=fig.add_gridspec(2,1,height_ratios=[1,3.2],hspace=0.13)
ax=fig.add_subplot(gs[0])
ax.plot(range(len(REG)),[mean[RI[r],MI["n_var"]] for r in REG],"-o",color="#2E5E8E",ms=3)
ax.axvline(RI["MID"],color="k",lw=1)
ax.set_xlim(-0.5,len(REG)-0.5); ax.set_xticks([]); ax.set_ylabel("variants per unit",fontsize=9)
ax.set_title("A  Total variants per unit by region (two-stage mean)",loc="left",fontsize=10)
ax=fig.add_subplot(gs[1])
im=ax.imshow(Mt,aspect="auto",cmap="RdBu_r",vmin=-1.5,vmax=1.5)
ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows,fontsize=8)
tk=[RI["P%d"%d] for d in range(1,31,3)]+[RI["MID"]]+[RI["D%d"%d] for d in range(30,0,-3)]
ax.set_xticks(tk); ax.set_xticklabels(["P%d"%d for d in range(1,31,3)]+["MID"]+["D%d"%d for d in range(30,0,-3)],fontsize=8,rotation=90)
ax.axvline(RI["MID"],color="k",lw=1)
plt.colorbar(im,ax=ax,fraction=0.018,pad=0.01,label="log2( region / middle )")
ax.set_title("B  What kind of variant each region carries, relative to the middle\n"
             "rows: class (SNV/del/ins, _rep = inside a repeat block) · substitution type · CpG context · element of the 2,240 bp unit",
             loc="left",fontsize=10)
fig.suptitle("5S rDNA — variant profile by region (P1..P30 | middle | D30..D1) · v11 · genomic polarity · %s set"%PRESET,fontsize=12,y=0.95)
p=os.path.join(G.OUT,"fig_v11_region_profile_%s.png"%PRESET)
fig.savefig(p,dpi=150,bbox_inches="tight"); print("saved:",p)
