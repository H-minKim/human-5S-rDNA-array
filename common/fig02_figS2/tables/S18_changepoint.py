# -*- coding: utf-8 -*-
"""
S18: segments **each array separately** (change-point detection) and looks at which segment grows as arrays get longer.

Why global clustering is not enough
  S13 to S17 clustered units globally and then looked at positions, so they never asked directly 'where is the boundary of this array'.
  Here the units of each haplotype are laid out in proximal -> distal order and **change points are searched within that array.**

Algorithm
  - Features: 20 variant characters (B) / 160 individual variants (A). Centered within haplotype, then 5 components by global PCA.
  - Segmentation: binary segmentation. Cost = sum of squared deviations within segments. A cut is accepted only if the cost reduction
           exceeds the penalty pen = PEN x d x log(N) (BIC form). Minimum segment length 3.
  - Null: shuffle the unit order of the same array and run the same procedure. If the order has no structure, no segments should appear.
  - PEN is run at several values to check whether the conclusion depends on the penalty.

Decision
  Boundary positions are recorded as absolute distances: from the proximal junction to the first boundary (= proximal segment length),
  from the distal junction to the last boundary (= distal segment length), and the sum of inner segment lengths.
  -> Each is regressed on N. If the ends are fixed the slope is 0, if everything is proportional the slope equals the share of that segment.
Output: S18_changepoint.tsv
"""
import sys,io,os,collections
import numpy as np
import v11_io as G

PRESET=sys.argv[1] if len(sys.argv)>1 else "primary"
OUT=os.path.join(G.OUT,PRESET)
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"S16_abs_vs_rel.py"),encoding="utf8").read()
exec(compile(src.split("STR=[(40,80")[0],"head","exec"))   # units, ap, ad, Nn, hid, XA, XB, cs(), HK, U

rng=np.random.default_rng(20260913)

def pca(X,k=5):
    Xc=X-X.mean(0)
    C=np.cov(Xc.T)
    w,v=np.linalg.eigh(C)
    idx=np.argsort(w)[::-1][:k]
    P=Xc@v[:,idx]
    return P/np.maximum(P.std(0),1e-9)

def binseg(X,pen,min_size=3):
    """Binary segmentation. Returns: sorted inner boundary positions (1-based, 'number of units before the boundary')"""
    n=len(X)
    def cost(a,b):
        if b-a<=0: return 0.0
        S=X[a:b]; mu=S.mean(0)
        return float(((S-mu)**2).sum())
    segs=[(0,n,cost(0,n))]
    cps=[]
    while True:
        best=None
        for si,(a,b,c0) in enumerate(segs):
            if b-a<2*min_size: continue
            # fast, using cumulative sums
            S=X[a:b]
            cs1=np.cumsum(S,0); cs2=np.cumsum(S*S,0)
            tot1=cs1[-1]; tot2=cs2[-1]; m=b-a
            bestg=None
            for t in range(min_size,m-min_size+1):
                l1=cs1[t-1]; l2=cs2[t-1]
                cl=float((l2-l1*l1/t).sum())
                r1=tot1-l1; r2=tot2-l2; rr=m-t
                cr=float((r2-r1*r1/rr).sum())
                g=c0-(cl+cr)
                if bestg is None or g>bestg[0]: bestg=(g,t,cl,cr)
            if bestg and (best is None or bestg[0]>best[0]):
                best=(bestg[0],si,bestg[1],bestg[2],bestg[3])
        if best is None or best[0]<=pen: break
        g,si,t,cl,cr=best
        a,b,c0=segs.pop(si)
        segs.append((a,a+t,cl)); segs.append((a+t,b,cr))
        cps.append(a+t)
        segs.sort()
    return sorted(cps)

ORD={}
for h in range(NH):
    idx=np.where(hid==h)[0]
    ORD[h]=idx[np.argsort(ap[idx])]          # proximal 1 -> N

Nh=np.array([U[k][0][2] for k in HK],dtype=float)
g=io.open(os.path.join(OUT,"S18_changepoint.tsv"),"w",encoding="utf8",newline="\n")
g.write("basis\tpen\thap\tN\tn_seg\tprox_seg_len\tdist_seg_len\tinner_len\tboundaries\n")

for tag,X in (("B",XB),("A",XA)):
    Z=pca(cs(X),5)
    for PEN in (1.0,2.0,4.0):
        rec=[];nullrec=[]
        for h in range(NH):
            idx=ORD[h]; n=len(idx)
            if n<40: continue
            Xi=Z[idx]
            pen=PEN*Xi.shape[1]*np.log(n)
            cps=binseg(Xi,pen)
            nseg=len(cps)+1
            pl=cps[0] if cps else n
            dl=n-cps[-1] if cps else n
            inner=n-pl-dl if cps and len(cps)>=2 else 0
            rec.append((Nh[h],nseg,pl,dl,inner))
            g.write("%s\t%.1f\t%s_%s\t%d\t%d\t%d\t%d\t%d\t%s\n"%(tag,PEN,HK[h][0],HK[h][1],n,nseg,pl,dl,inner,
                    ",".join(str(c) for c in cps)))
            # null: shuffle the order
            Xp=Xi[rng.permutation(n)]
            nullrec.append(len(binseg(Xp,pen))+1)
        R=np.array(rec,float)
        if len(R)<50: continue
        Nv=R[:,0]
        print("\n"+"="*92)
        print("=== 기준 %s · 벌점 PEN=%.1f · 하플로타입 %d ==="%(tag,PEN,len(R)))
        print("  구간 수 : 실제 평균 %.2f (중앙 %.0f) · 순서 섞기 귀무 평균 %.2f"%(
            R[:,1].mean(),np.median(R[:,1]),float(np.mean(nullrec))))
        if np.mean(nullrec)>=R[:,1].mean()*0.9:
            print("  ⚠ 귀무와 비슷하다 = 순서 구조를 잡은 게 아니라 잡음을 자른 것. 이 벌점은 못 쓴다.")
        A2=np.vstack([np.ones_like(Nv),Nv]).T
        print("")
        print("  무엇이 N 과 함께 늘어나는가 (기울기 · 비례라면 그 구간의 몫과 같아야 한다)")
        print("   항목              평균     기울기 b [95%CI]        절편    비례기대   판정")
        for ci,nm in ((1,"구간 수"),(2,"근위 첫 구간 길이"),(3,"원위 마지막 구간 길이"),(4,"내부 구간 합")):
            y=R[:,ci]
            b=np.linalg.lstsq(A2,y,rcond=None)[0]
            bs=[np.linalg.lstsq(A2[ii],y[ii],rcond=None)[0][1] for ii in (rng.integers(0,len(R),len(R)) for _ in range(400))]
            lo,hi=np.percentile(bs,[2.5,97.5])
            exp=(y.mean()/Nv.mean()) if ci!=1 else float("nan")
            if ci==1: jud="—"
            elif hi<exp*0.5: jud="고정에 가까움"
            elif lo>exp*1.5: jud="비례보다 빠름"
            else: jud="비례 성장"
            print("   %-16s %7.1f   %+.4f [%+.4f,%+.4f]  %6.1f   %s   %s"%(
                nm,y.mean(),b[1],lo,hi,b[0],("%.4f"%exp) if exp==exp else "   -  ",jud))
g.close()
print("\n저장: %s"%os.path.join(OUT,"S18_changepoint.tsv"))
