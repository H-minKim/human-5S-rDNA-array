# -*- coding: utf-8 -*-
"""
S20 - finds change points (a) without bias and (b) with a fixed number of change points.

Change-point method (exact definition)
  Order the units of a haplotype as proximal 1 -> N and build a d-dimensional feature vector x_i for each unit.
  Cost of segment [a,b) = sum_{i in [a,b)} ||x_i - mean(x_[a,b))||^2  (deviation from the segment mean)
  Splitting the whole at one point t lowers the cost. Pick the t with the largest decrease, then
  (a) unbiased: accept only if the decrease exceeds the penalty pen = PEN * d * log(N), and repeat until none is left
  (b) fixed K: ignore the penalty and repeat the split with the largest decrease K times
  Minimum segment length 3. Features: centering within haplotype -> global PCA with 5 components (noise reduction).
  Null: shuffle the unit order of the same array and apply the same procedure. If the order has no structure, no change points should appear.
Output: S20_cp_unbiased.tsv, S20_cp_fixedK.tsv, fig_v11_changepoints_<preset>.png
"""
import sys,io,os,collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import v11_io as G

PRESET=sys.argv[1] if len(sys.argv)>1 else "primary"
OUT=os.path.join(G.OUT,PRESET)
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"S16_abs_vs_rel.py"),encoding="utf8").read()
exec(compile(src.split("STR=[(40,80")[0],"head","exec"))
rng=np.random.default_rng(20260913)

def pca(X,k=5):
    Xc=X-X.mean(0); C=np.cov(Xc.T); w,v=np.linalg.eigh(C)
    P=Xc@v[:,np.argsort(w)[::-1][:k]]
    return P/np.maximum(P.std(0),1e-9)
ZB=pca(cs(XB),5); ZA=pca(cs(XA),5)

def best_split(S,min_size=3):
    m=len(S)
    if m<2*min_size: return None
    cs1=np.cumsum(S,0); cs2=np.cumsum(S*S,0); tot1=cs1[-1]; tot2=cs2[-1]
    mu=tot1/m; c0=float((tot2-2*mu*tot1+m*mu*mu).sum())
    best=None
    for t in range(min_size,m-min_size+1):
        l1=cs1[t-1]; l2=cs2[t-1]
        cl=float((l2-l1*l1/t).sum())
        r1=tot1-l1; r2=tot2-l2; rr=m-t
        cr=float((r2-r1*r1/rr).sum())
        gain=c0-(cl+cr)
        if best is None or gain>best[0]: best=(gain,t,cl,cr)
    return best
def binseg(X,pen=None,K=None,min_size=3):
    n=len(X); segs=[(0,n)]; cps=[]
    while True:
        if K is not None and len(cps)>=K: break
        cand=None
        for si,(a,b) in enumerate(segs):
            r=best_split(X[a:b],min_size)
            if r is None: continue
            if cand is None or r[0]>cand[0]: cand=(r[0],si,r[1])
        if cand is None: break
        if pen is not None and cand[0]<=pen: break
        g,si,t=cand
        a,b=segs.pop(si); segs.append((a,a+t)); segs.append((a+t,b)); segs.sort()
        cps.append(a+t)
    return sorted(cps)

ORD={h:np.where(hid==h)[0][np.argsort(ap[np.where(hid==h)[0]])] for h in range(NH)}
Nh=np.array([U[k][0][2] for k in HK],dtype=float)

# ---------- (a) unbiased ----------
print("=== (가) 비편향 변화점 (벌점 PEN*d*logN) ===")
UNB={}
g=io.open(os.path.join(OUT,"S20_cp_unbiased.tsv"),"w",encoding="utf8",newline="\n")
g.write("basis\tpen\thap\tN\tn_cp\tcps_from_prox\n")
for tag,Z in (("B",ZB),("A",ZA)):
    for PEN in (1.0,2.0):
        ncp=[];Nv=[];pos=[];nullc=[]
        for h in range(NH):
            idx=ORD[h]; n=len(idx)
            if n<40: continue
            Xi=Z[idx]; pen=PEN*Xi.shape[1]*np.log(n)
            c=binseg(Xi,pen=pen)
            ncp.append(len(c)); Nv.append(Nh[h]); pos.extend([(x,n) for x in c])
            nullc.append(len(binseg(Xi[rng.permutation(n)],pen=pen)))
            g.write("%s\t%.1f\t%s_%s\t%d\t%d\t%s\n"%(tag,PEN,HK[h][0],HK[h][1],n,len(c),",".join(map(str,c))))
        ncp=np.array(ncp,float); Nv=np.array(Nv,float)
        A2=np.vstack([np.ones_like(Nv),Nv]).T
        b=np.linalg.lstsq(A2,ncp,rcond=None)[0]
        bs=[np.linalg.lstsq(A2[ii],ncp[ii],rcond=None)[0][1] for ii in (rng.integers(0,len(Nv),len(Nv)) for _ in range(400))]
        lo,hi=np.percentile(bs,[2.5,97.5])
        print("  [%s PEN=%.1f] 변화점 개수 평균 %.2f (중앙 %.0f · 최대 %d) · 귀무 %.2f"%(
            tag,PEN,ncp.mean(),np.median(ncp),int(ncp.max()),float(np.mean(nullc))))
        print("     개수 vs N : 기울기 %+.4f [%+.4f,%+.4f] · 절편 %.2f  (N 100 늘 때 %+.2f 개)"%(
            b[1],lo,hi,b[0],b[1]*100))
        # mean count per length stratum
        s=[]
        for lo_,hi_,lab_ in ((40,80,"N40-80"),(81,110,"N81-110"),(111,140,"N111-140"),(141,10**9,"N141+")):
            m=(Nv>=lo_)&(Nv<=hi_)
            if m.sum()>=10: s.append("%s %.2f"%(lab_,ncp[m].mean()))
        print("     길이층별 평균 : "+" · ".join(s))
        UNB[(tag,PEN)]=(ncp,Nv,pos)
g.close()

# ---------- (b) fixed K ----------
print("")
print("=== (나) 변화점 개수를 고정 (벌점 없이 K 번 분할) ===")
FIX={}
g=io.open(os.path.join(OUT,"S20_cp_fixedK.tsv"),"w",encoding="utf8",newline="\n")
g.write("basis\tK\thap\tN\tcp_index\tfrom_prox\tfrom_dist\trel\n")
for tag,Z in (("B",ZB),("A",ZA)):
    for K in (1,2,3):
        rows=[]
        for h in range(NH):
            idx=ORD[h]; n=len(idx)
            if n<40: continue
            c=binseg(Z[idx],K=K)
            if len(c)<K: continue
            for j,x in enumerate(sorted(c)):
                rows.append((Nh[h],j,x,n-x,x/float(n)))
                g.write("%s\t%d\t%s_%s\t%d\t%d\t%d\t%d\t%.4f\n"%(tag,K,HK[h][0],HK[h][1],n,j+1,x,n-x,x/float(n)))
        R=np.array(rows,float)
        FIX[(tag,K)]=R
        print("")
        print("  [%s · K=%d] 하플로타입 %d"%(tag,K,len(set(R[:,0].tolist())) if len(R) else 0))
        print("   경계   근위에서(중앙)  원위에서(중앙)  상대(중앙)  |  근위거리~N 기울기   상대~N 기울기(x100)")
        for j in range(K):
            m=R[:,1]==j
            if m.sum()<30: continue
            Nv=R[m,0]; fp=R[m,2]; rl=R[m,4]
            A2=np.vstack([np.ones_like(Nv),Nv]).T
            b1=np.linalg.lstsq(A2,fp,rcond=None)[0][1]
            b2=np.linalg.lstsq(A2,rl,rcond=None)[0][1]*100
            print("    %d 번   %8.0f       %8.0f       %.3f     |   %+.3f            %+.3f"%(
                j+1,np.median(fp),np.median(R[m,3]),np.median(rl),b1,b2))
g.close()

# ---------- figure ----------
fig=plt.figure(figsize=(16,11)); gs=fig.add_gridspec(2,2,hspace=0.33,wspace=0.25)
ncp,Nv,pos=UNB[("B",1.0)]
ax=fig.add_subplot(gs[0,0])
jit=rng.normal(0,0.08,len(Nv))
ax.plot(Nv,ncp+jit,".",ms=3,alpha=0.45,color="#2E5E8E")
bins=[(40,70),(71,90),(91,110),(111,130),(131,160),(161,260)]
bx=[];by=[];be=[]
for lo,hi in bins:
    m=(Nv>=lo)&(Nv<=hi)
    if m.sum()<10: continue
    bx.append(np.median(Nv[m])); by.append(ncp[m].mean()); be.append(ncp[m].std()/np.sqrt(m.sum()))
ax.errorbar(bx,by,yerr=be,fmt="o-",color="#B03A2E",ms=6,lw=2,capsize=3,label="mean per length bin")
ax.set_xlabel("array length N (units)",fontsize=9); ax.set_ylabel("number of change-points found",fontsize=9)
ax.set_title("A  Unbiased: how many change-points does an array get?\nBIC-style penalty, basis = variant character",loc="left",fontsize=10)
ax.legend(fontsize=8)
ax=fig.add_subplot(gs[0,1])
P=np.array(pos,float)
ax.hist(P[:,0],bins=np.arange(0,121,3),alpha=0.75,color="#2E5E8E",label="distance from proximal PJ")
ax.hist(P[:,1]-P[:,0],bins=np.arange(0,121,3),alpha=0.55,color="#9C3B32",label="distance from distal DJ")
ax.set_xlabel("change-point position (units from that junction)",fontsize=9)
ax.set_ylabel("count",fontsize=9)
ax.set_title("B  Where do the change-points fall? (absolute)\nno pile-up near either junction = boundaries are not junction-anchored",loc="left",fontsize=10)
ax.legend(fontsize=8)
ax=fig.add_subplot(gs[1,0])
ax.hist(P[:,0]/P[:,1],bins=np.arange(0,1.01,0.025),color="#4F7942",alpha=0.8)
ax.axvline(0.5,color="k",ls=":",lw=1)
ax.set_xlabel("change-point position (fraction of array length)",fontsize=9)
ax.set_ylabel("count",fontsize=9)
ax.set_title("C  Same change-points in relative coordinate\nbroad, centred — consistent with a gradient, not a fixed domain edge",loc="left",fontsize=10)
ax=fig.add_subplot(gs[1,1])
for K,c,mk in ((1,"#2E5E8E","o"),(2,"#B03A2E","s"),(3,"#4F7942","^")):
    R=FIX[("B",K)]
    for j in range(K):
        m=R[:,1]==j
        if m.sum()<30: continue
        Nv2=R[m,0]; fp=R[m,2]
        o=np.argsort(Nv2)
        w=21
        sm=np.convolve(fp[o],np.ones(w)/w,mode="valid")
        ax.plot(Nv2[o][w//2:len(sm)+w//2],sm,"-",color=c,lw=1.6,
                label=("K=%d"%K) if j==0 else None)
ax.plot([40,230],[20,115],"k:",lw=1,label="proportional (N/2)")
ax.set_xlabel("array length N",fontsize=9); ax.set_ylabel("boundary position (units from proximal PJ)",fontsize=9)
ax.set_title("D  Fixed number of change-points: where do the boundaries sit?\nrunning mean; boundaries track N (relative), they do not stay at a fixed distance",loc="left",fontsize=10)
ax.legend(fontsize=8)
fig.suptitle("5S rDNA — change-point segmentation of each array (unit sequence only)\n"
             "v11 · %s · within-haplotype centred · PCA 5 · binary segmentation"%PRESET,fontsize=12,y=0.975)
p=os.path.join(G.OUT,"fig_v11_changepoints_%s.png"%PRESET)
fig.savefig(p,dpi=150,bbox_inches="tight"); print("\nsaved:",p)
