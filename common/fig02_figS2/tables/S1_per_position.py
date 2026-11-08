# -*- coding: utf-8 -*-
"""
S1 - Result 1-1 (coding part variants), 1-2 (all variants, SNV / indel separately), by unit position
Reference build : v11 (canonical cmp4d_v10 + 03_maps/g8canonv10/CONSENSUS_MAP_G8_full.tsv, gate v11_io)
Compartments    : proximal1 .. proximal30, middle (more than 30 units from both junctions), distal30 .. distal1   (user specified)

Carried over identically from the old code (verified)
  metric definitions   binvar.py lines 8-30
      n_var  = number of unit tokens
      n_snv  = snv + snv_rep            n_del = del + del_rep        n_ins = ins + ins_rep
      n_cod = number of tokens with position >= 2122 (blocks merged into DB/IB have no position and are not counted as coding)
      n_cod_x= n_cod excluding '2205TC'(hs5S U84C)
      n_del_ge10 = deletions (including DB) with length >= 10   (del10 of fig18b.py / binvar4.py)
  two-stage aggregation  binvar2.py profile() : ① mean within haplotype ② mean across haplotypes, bootstrap 800, seed 4242
                    min_hap = 150 (all) / 40 (length strata)
  length strata    binvar.py STR fixed bins N40-80 / N81-110 / N111+  (N<40 is dropped, the count is reported)
  significance     position-dependence statistic T=|mean(1-8)-mean(rest)| of binvar4.py / shiftdemo.py, exhaustive circular shifts
                    (the segment-SSE version of binvar2.py failed the positive control and was discarded)
                    30 positions -> exact p lower bound = 1/30

Newly added (not in the old code, marked)
  n_cod_snv, n_cod_indel     : coding part split into SNV / indel (requested to be shown separately, like 1-2)
  single-position specific variants : the zone-specific call of P2 (paired haplotype difference + haplotype bootstrap 95% CI)
                                 narrowed from the 'terminal 1-8 bundle' to 'a single position'.
                                 d_h = 1[the unit at that position carries the variant] - (carrier fraction among the middle units of that haplotype)
                                 bootstrap uses multinomial weights (equivalent to row resampling with replacement), 2000 times, seed 20260830
                                 Because this is multiple testing of 60 positions x hundreds of variants, BH q is reported too.
                                   Counting with only the old criterion (CI excludes 0) mixes in chance positives.
"""
import sys,io,os,re,collections
import numpy as np
import v11_io as G

PRESET=sys.argv[1] if len(sys.argv)>1 else "primary"
OUT=os.path.join(G.OUT,PRESET); os.makedirs(OUT,exist_ok=True)
U,A,M=G.load(PRESET)
MAXD=30; MID_GATE=30
STR=[(40,80,"N40-80"),(81,110,"N81-110"),(111,400,"N111+")]
MIN_ZS=30

CAN=G.CAN
def ispos(x):                                   # as in binvar.py lines 8-12
    if x.startswith("DB") or x.startswith("IB"): return None
    if x[0] in "DI": return int(x[1:].split(":")[0])
    if x[0]=="R": return int(x[1:-2])
    return int(x[:-2])
def cls(x):                                     # as in binvar.py lines 13-19
    if x.startswith("DB"): return "del_rep"
    if x.startswith("IB"): return "ins_rep"
    if x[0]=="D": return "del"
    if x[0]=="I": return "ins"
    if x[0]=="R": return "snv_rep"
    return "snv"
def dlen(x):                                    # as in fig18b.py dlen
    if x.startswith("DB") or x.startswith("IB"): return int(x[2:].split(":")[1])
    if x[0] in "DI": return int(x[1:].split(":")[1])
    return 0
def element(p):
    if p is None: return ""
    for a,b,n in [(1,318,"AluY"),(319,1178,"NTS1"),(1179,1551,"GArich"),(1552,2121,"NTS2"),(2122,2240,"coding")]:
        if a<=p<=b: return n
    return ""
def iscpg(p):
    if p is None: return False
    c=CAN[p-1]
    return (c=="C" and p<2240 and CAN[p]=="G") or (c=="G" and p>=2 and CAN[p-2]=="C")

def metrics(vs):
    c=collections.Counter(cls(x) for x in vs)
    cod=[x for x in vs if (ispos(x) or 0)>=2122]
    return dict(
        n_var=float(len(vs)),
        n_snv=float(c["snv"]+c["snv_rep"]),
        n_del=float(c["del"]+c["del_rep"]),
        n_ins=float(c["ins"]+c["ins_rep"]),
        n_indel=float(c["del"]+c["del_rep"]+c["ins"]+c["ins_rep"]),
        n_del_ge10=float(sum(1 for x in vs if x[0]=="D" and dlen(x)>=10)),
        n_cod=float(len(cod)),
        n_cod_x=float(sum(1 for x in cod if x!="2205TC")),
        n_cod_snv=float(sum(1 for x in cod if cls(x) in ("snv","snv_rep"))),
        n_cod_indel=float(sum(1 for x in cod if cls(x) in ("del","ins"))),
        u84c=float("2205TC" in vs))
METRICS=["n_var","n_snv","n_del","n_ins","n_indel","n_del_ge10",
         "n_cod","n_cod_x","n_cod_snv","n_cod_indel","u84c"]
ROWS={}
for k,rs in U.items():
    ROWS[k]=[(r[0],r[1],r[2],metrics(r[4])) for r in rs]

rng=np.random.default_rng(4242)                 # binvar2.py
def profile(col,side,maxp=MAXD,min_hap=150,nlo=0,nhi=10**9,nb=800):
    out=[]
    for d in range(1,maxp+1):
        g=[]
        for k,rs in ROWS.items():
            if not (nlo<=rs[0][2]<=nhi): continue
            v=[r[3][col] for r in rs if (r[0] if side=="prox" else r[1])==d]
            if v: g.append(np.mean(v))
        if len(g)<min_hap: continue
        g=np.array(g,float)
        b=[g[rng.integers(0,len(g),len(g))].mean() for _ in range(nb)]
        out.append((d,len(g),g.mean(),g.std(ddof=1),np.percentile(b,2.5),np.percentile(b,97.5)))
    return out
def mid_level(col,nlo=0,nhi=10**9,nb=800):
    g=[]
    for k,rs in ROWS.items():
        if not (nlo<=rs[0][2]<=nhi): continue
        v=[r[3][col] for r in rs if r[0]>MID_GATE and r[1]>MID_GATE]
        if v: g.append(np.mean(v))
    if not g: return (np.nan,0,np.nan,np.nan)
    g=np.array(g,float)
    b=[g[rng.integers(0,len(g),len(g))].mean() for _ in range(nb)]
    return (g.mean(),len(g),np.percentile(b,2.5),np.percentile(b,97.5))
def end_ratio(col,side,nlo=0,nhi=10**9,nb=2000):
    """Terminal 1-8 / middle, paired by haplotype, ratio of means, haplotype bootstrap"""
    e=[];m=[]
    for k,rs in ROWS.items():
        if not (nlo<=rs[0][2]<=nhi): continue
        a=[r[3][col] for r in rs if (r[0] if side=="prox" else r[1])<=8]
        c=[r[3][col] for r in rs if r[0]>MID_GATE and r[1]>MID_GATE]
        if a and c: e.append(np.mean(a)); m.append(np.mean(c))
    e=np.array(e); m=np.array(m)
    if len(e)<40 or m.mean()<=0: return (np.nan,np.nan,np.nan,len(e))
    bs=[]
    for _ in range(nb):
        j=rng.integers(0,len(e),len(e))
        if m[j].mean()>0: bs.append(e[j].mean()/m[j].mean())
    return (e.mean()/m.mean(),np.percentile(bs,2.5),np.percentile(bs,97.5),len(e))
def shift_test(y):
    y=np.asarray(y,float); n=len(y)
    T=lambda v: abs(v[:8].mean()-v[8:].mean())
    t0=T(y); ge=sum(1 for s in range(1,n) if T(np.roll(y,s))>=t0)
    return t0,(ge+1)/float(n),n

# ================== 1-1, 1-2 per-position profile ==================
Ns=np.array([rs[0][2] for rs in ROWS.values()])
print("\n배열 길이: 슬롯 %d · N<40 %d 슬롯은 고정 길이층에서 빠진다"%(len(Ns),int((Ns<40).sum())))
g=io.open(os.path.join(OUT,"S1_per_position_profile.tsv"),"w",encoding="utf8",newline="\n")
g.write("stratum\tmetric\tside\tunit_from_junction\tkb_from_junction\tn_haplotypes\tmean\tsd_across_haplotypes\tci_lo\tci_hi\n")
SUM={}
for lo,hi,lab in [(0,10**9,"all")]+STR:
    for col in METRICS:
        for side in ("prox","dist"):
            pr=profile(col,side,MAXD,150 if lab=="all" else 40,lo,hi)
            SUM[(lab,col,side)]=pr
            for d,n,mn,sd,cl,ch in pr:
                g.write("%s\t%s\t%s\t%d\t%.2f\t%d\t%.5f\t%.5f\t%.5f\t%.5f\n"%(lab,col,side,d,d*2.24,n,mn,sd,cl,ch))
        mm=mid_level(col,lo,hi)
        SUM[(lab,col,"MID")]=mm
        g.write("%s\t%s\tMID\t-1\t-1\t%d\t%.5f\t\t%.5f\t%.5f\n"%(lab,col,mm[1],mm[0],mm[2],mm[3]))
g.close()

print("\n=== 1-1 · 1-2 말단 vs 가운데 (분모 = 그 자리의 유닛 1 개 · 2단계 집계) ===")
print("  %d 개체 · %d 슬롯 · %d 유닛 (preset=%s)"%(M["n_samples"],M["n_slots"],M["n_units"],PRESET))
print("  지표          근위1-8  가운데   원위1-8 | 근위/가운데 [95%CI]      원위/가운데 [95%CI]      | 순환이동 p 근위 원위")
RT=io.open(os.path.join(OUT,"S1_end_vs_mid.tsv"),"w",encoding="utf8",newline="\n")
RT.write("stratum\tmetric\tside\tend1_8_over_mid\tci_lo\tci_hi\tn_hap\tshift_T\tshift_p_exact\tp_floor\n")
for lo,hi,lab in [(0,10**9,"all")]+STR:
    if lab!="all": print("  --- 길이층 %s"%lab)
    for col in METRICS:
        pp=SUM[(lab,col,"prox")]; dd=SUM[(lab,col,"dist")]; mm=SUM[(lab,col,"MID")]
        if len(pp)<MAXD or len(dd)<MAXD:
            if lab=="all": print("  %-12s 자리가 %d/%d 만 채워져 순환이동 생략"%(col,len(pp),len(dd)))
            continue
        a=np.mean([x[2] for x in pp[:8]]); b=np.mean([x[2] for x in dd[:8]])
        rp=end_ratio(col,"prox",lo,hi); rd=end_ratio(col,"dist",lo,hi)
        tp,pP,n=shift_test([x[2] for x in pp]); td,pD,_=shift_test([x[2] for x in dd])
        for side,r,t,p in (("prox",rp,tp,pP),("dist",rd,td,pD)):
            RT.write("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%d\t%.5f\t%.4f\t%.4f\n"%(lab,col,side,r[0],r[1],r[2],r[3],t,p,1.0/n))
        if lab=="all" or col in ("n_var","n_del_ge10","n_cod_x"):
            print("  %-12s %7.3f %7.3f %8.3f | %6.3f [%5.3f,%5.3f] %6.3f [%5.3f,%5.3f] | %.4f %.4f"%(
                col,a,mm[0],b,rp[0],rp[1],rp[2],rd[0],rd[1],rd[2],pP,pD))
RT.close()
print("  정확 p 하한 = 1/%d = %.4f"%(MAXD,1.0/MAXD))

print("\n=== 자리 하나하나 — 요청 구획 (전체) ===")
print("  자리   n_var 근위 / 원위     n_snv 근위 / 원위    n_indel 근위 / 원위   coding(U84C제외) 근위 / 원위   U84C 보유비율 근위 / 원위")
P={c:{x[0]:x[2] for x in SUM[("all",c,"prox")]} for c in METRICS}
D={c:{x[0]:x[2] for x in SUM[("all",c,"dist")]} for c in METRICS}
for d in range(1,MAXD+1):
    if d not in P["n_var"] or d not in D["n_var"]: continue
    print("  %-4d %7.3f / %6.3f   %6.3f / %6.3f   %6.3f / %6.3f     %6.4f / %6.4f           %6.4f / %6.4f"%(
        d,P["n_var"][d],D["n_var"][d],P["n_snv"][d],D["n_snv"][d],P["n_indel"][d],D["n_indel"][d],
        P["n_cod_x"][d],D["n_cod_x"][d],P["u84c"][d],D["u84c"][d]))
print("  가운데 %6.3f          %6.3f           %6.3f            %6.4f                  %6.4f"%(
    SUM[("all","n_var","MID")][0],SUM[("all","n_snv","MID")][0],SUM[("all","n_indel","MID")][0],
    SUM[("all","n_cod_x","MID")][0],SUM[("all","u84c","MID")][0]))

# ================== position-specific variants at single-position resolution ==================
NU=collections.Counter()
for rs in U.values():
    for r in rs: NU.update(r[4])
VX=sorted([x for x,c in NU.items() if c>=MIN_ZS],key=lambda x:-NU[x])
XI={x:i for i,x in enumerate(VX)}
print("\n=== 자리 하나 단위 특이 변이 — 대상 변이 %d 종 (보유 유닛 >= %d, P2 기준) ==="%(len(VX),MIN_ZS))
# per-haplotype carrier fraction in the middle
HK=[k for k,rs in U.items() if sum(1 for r in rs if r[0]>MID_GATE and r[1]>MID_GATE)>=2]   # P2: MID >= 2 units
FM={}
for k in HK:
    mids=[r[4] for r in U[k] if r[0]>MID_GATE and r[1]>MID_GATE]
    v=np.zeros(len(VX),np.float32)
    for s in mids:
        for x in s:
            j=XI.get(x)
            if j is not None: v[j]+=1
    FM[k]=v/len(mids)
brng=np.random.default_rng(20260830); B=2000
out=[]
for side in ("prox","dist"):
    for d in range(1,MAXD+1):
        hs=[];A=[]
        for k in HK:
            us=[r for r in U[k] if (r[0] if side=="prox" else r[1])==d]
            if not us: continue
            a=np.zeros(len(VX),np.float32)
            for x in us[0][4]:
                j=XI.get(x)
                if j is not None: a[j]=1
            hs.append(k); A.append(a)
        n=len(hs)
        if n<100: continue
        A=np.array(A); Fm=np.array([FM[k] for k in hs])
        Dm=A-Fm
        W=brng.multinomial(n,np.full(n,1.0/n),size=B).astype(np.float32)
        BS=(W@Dm)/n
        mean=Dm.mean(0); lo=np.percentile(BS,2.5,axis=0); hi=np.percentile(BS,97.5,axis=0)
        pneg=(np.sum(BS<=0,axis=0)+1)/(B+1.0); ppos=(np.sum(BS>=0,axis=0)+1)/(B+1.0)
        pb=np.minimum(1.0,2*np.minimum(pneg,ppos))
        fz=A.mean(0); fm=Fm.mean(0)
        for j,x in enumerate(VX):
            if fz[j]==0 and fm[j]==0: continue
            out.append([side,d,x,n,float(fz[j]),float(fm[j]),float(mean[j]),float(lo[j]),float(hi[j]),float(pb[j])])
# BH
pv=np.array([o[9] for o in out]); order=np.argsort(pv); m=len(pv)
q=np.empty(m); run=1.0
for rank,i in enumerate(order[::-1]):
    rk=m-rank
    run=min(run,pv[i]*m/rk); q[i]=run
for o,qq in zip(out,q): o.append(float(qq))
g=io.open(os.path.join(OUT,"S1_per_position_specific_variants.tsv"),"w",encoding="utf8",newline="\n")
g.write("side\tunit_from_junction\tvariant\tclass\tpos\tblock_span\telement\tcpg\tn_hap\t"
        "freq_at_position\tfreq_mid_hapmean\tdelta\tci_lo\tci_hi\tp_boot\tq_bh\tdirection\tci_excl0\n")
CNT=collections.Counter(); TOP=collections.defaultdict(list)
for side,d,x,n,fz,fm,dl,lo,hi,pb,qq in out:
    ex=(lo>0 or hi<0)
    if not (ex or qq<0.05): continue
    p=ispos(x); c=cls(x)
    span=""
    if x.startswith(("DB","IB")):
        b=int(x[2:].split(":")[0]); span="%d-%d"%G.BLOCKS[b-1]
    dr="up" if dl>0 else "down"
    g.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.5f\t%.5f\t%+.5f\t%+.5f\t%+.5f\t%.4g\t%.4g\t%s\t%s\n"%(
        side,d,x,c,("" if p is None else p),span,element(p),("Y" if iscpg(p) else ""),n,fz,fm,dl,lo,hi,pb,qq,dr,"Y" if ex else "N"))
    if ex: CNT[(side,d,dr,"ci")]+=1
    if qq<0.05:
        CNT[(side,d,dr,"q")]+=1; TOP[(side,d)].append((dl,x,c,fz,fm,qq))
g.close()
print("  검정 %d 건 (자리 x 변이) · CI 가 0 을 벗어난 것 %d · q<0.05 %d"%(
    len(out),sum(1 for o in out if o[7]>0 or o[8]<0),int((q<0.05).sum())))
print("  옛 기준(CI 만)으로 기대되는 우연 양성 ≈ 5%% x %d = %d"%(len(out),int(0.05*len(out))))
print("\n  자리    근위 농축/결핍 (q<0.05)   원위 농축/결핍 (q<0.05)")
for d in range(1,MAXD+1):
    print("  %-4d        %4d / %-4d               %4d / %-4d"%(d,CNT[("prox",d,"up","q")],CNT[("prox",d,"down","q")],
          CNT[("dist",d,"up","q")],CNT[("dist",d,"down","q")]))
print("\n  자리별 q<0.05 농축 상위 (근위1-3 · 원위1-3)")
for side in ("prox","dist"):
    for d in (1,2,3):
        t=sorted([z for z in TOP[(side,d)] if z[0]>0],key=lambda z:-z[0])[:8]
        print("   %s%d : %s"%("근위" if side=="prox" else "원위",d,
              ", ".join("%s(%s %.3f vs %.3f)"%(x,c,fz,fm) for dl,x,c,fz,fm,qq in t) or "-"))
print("\n저장: %s"%OUT)
