# -*- coding: utf-8 -*-
# Compares three zonings and builds and checks a combined grid.
#   (A)  the six methylation zones: P1 proximal 1-3, P2 proximal 4-10, P3 proximal 11-25 (soft), INT, D2 distal 11-22, D1 distal 1-10
#   (B)  the sequence zones: proximal 1-3, proximal 4-17, the interior (proximal 18+ and distal 9+), and distal 1-8
# (C) here : methylation proximal 9-11, distal 10-11 / sequence proximal 4-5, distal 7-8, distal 2 alone
import io,os,warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
OUT=__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"; T=OUT+"/02_tables"
LOG=io.open(OUT+"/04_logs/V34_ZONING.txt","w",encoding="utf8",newline="\n")
rng=np.random.RandomState(20260915)
def P(s=""):
    try: print(s)
    except UnicodeEncodeError: print(s.encode("ascii","replace").decode())
    LOG.write(s+"\n"); LOG.flush()
D=pd.read_parquet(T+"/L2_HAP_UNIT.parquet"); D=D[D.n_read>=5].copy()
D["sl"]=D["sample"].astype(str)+"_"+D.slot.astype(str)
ap=D.ap.values.astype(int); ad=D.ad.values.astype(int)
near_p=ap<=ad
P("="*104); P("boundary comparison of the three zonings and the combined grid check"); P("="*104); P()

# ---- the third proximal boundary: 25 or 17, the one described as soft
P("[1] the third proximal boundary : 25.5 (soft), 17, and 17|18 from methylation here")
P()
P("   length stratum  f_hi : prox14 prox16 prox18 prox20 prox25 prox30  interior")
for lo,hi in ((30,59),(60,99),(100,139),(140,300)):
    s=D[(D.N>=lo)&(D.N<=hi)&near_p]
    v=[s.loc[s.ap==k,"f_hi"].mean() for k in (14,16,18,20,25,30)]
    it=D[(D.N>=lo)&(D.N<=hi)&(D.ap>30)&(D.ad>30)]["f_hi"].mean()
    P("   N %3d-%3d              "%(lo,hi)+"  ".join("%.3f"%x for x in v)+"   %.3f"%it)
P()
P("   position where the interior value is reached (98% or more of interior), per length stratum")
for lo,hi in ((30,59),(60,99),(100,139),(140,300)):
    s=D[(D.N>=lo)&(D.N<=hi)&near_p]
    it=D[(D.N>=lo)&(D.N<=hi)&(D.ap>30)&(D.ad>30)]["f_hi"].mean()
    hit=None
    for k in range(4,31):
        v=s.loc[s.ap==k,"f_hi"].mean()
        if v==v and it==it and v>=0.98*it: hit=k; break
    P("   N %3d-%3d : interior value reached at proximal %s (interior %.3f)"%(lo,hi,hit,it))
P()

# ---- the combined grid
def mzone(a,d,np_):
    if np_:
        if a<=3: return "M1 근1-3"
        if a<=10: return "M2 근4-10"
        if a<=17: return "M3 근11-17"
        return "MI 내부"
    if d<=10: return "D1 원1-10"
    if d<=22: return "D2 원11-22"
    return "MI 내부"
def szone(a,d,np_):
    if np_:
        if a<=3: return "S1 근1-3"
        if a<=17: return "S2 근4-17"
        return "SI 내부"
    if d<=2: return "S0 원1-2"
    if d<=8: return "SD 원3-8"
    return "SI 내부"
D["MZ"]=[mzone(a,d,n) for a,d,n in zip(ap,ad,near_p)]
D["SZ"]=[szone(a,d,n) for a,d,n in zip(ap,ad,near_p)]
D["Z"]=D.MZ+" / "+D.SZ
P("[2] cross table of methylation zone x sequence zone (number of position-haplotypes)")
P()
ct=pd.crosstab(D.MZ,D.SZ)
P(ct.to_string())
P()
P("  -> filled cells %d = candidate integrated zones"%int((ct.values>0).sum()))
P()
FM=["f_hi","b_lo_mean","p_switch","blk_len_mean","r_lo_gene_sp"]
FS=["n_var_rare","n_var_cpg","n_del","unit_len"]
P("[3] properties per integrated zone (cells with n>=300 only)")
P()
P("  zone                          n       "+" ".join("%-9s"%c for c in FM)+" | "+" ".join("%-9s"%c for c in FS))
keep=[]
for z,s in D.groupby("Z"):
    if len(s)<300: continue
    keep.append(z)
    P("  %-28s %6s  "%(z,format(len(s),","))
      +" ".join("%-9.4f"%s[c].mean() for c in FM)+" | "+" ".join("%-9.3f"%s[c].mean() for c in FS))
P()
P("[4] do the two axes give independent information: methylation axis only / sequence axis only / both")
P()
for tgt,lab in (("f_hi","메틸 표적 f_hi"),("n_var_rare","서열 표적 희귀변이 수")):
    y=D[tgt].values.astype(float); ok=np.isfinite(y)
    tot=((y[ok]-y[ok].mean())**2).sum()
    def ev(col):
        m=D[col].values
        pred=pd.Series(y[ok]).groupby(pd.Series(m[ok]).values).transform("mean").values
        return 1-((y[ok]-pred)**2).sum()/tot
    P("   %-22s methylation axis %.4f, sequence axis %.4f, both %.4f  (added %+.4f / %+.4f)"
      %(lab,ev("MZ"),ev("SZ"),ev("Z"),ev("Z")-ev("MZ"),ev("Z")-ev("SZ")))
P()
P("[5] cells created by boundary mismatch: distal 9-10 (methylation says 'off domain', sequence says 'interior')")
P()
a=D[(~near_p)&(D.ad>=9)&(D.ad<=10)]; b=D[(~near_p)&(D.ad>=3)&(D.ad<=8)]
c=D[(~near_p)&(D.ad>=11)&(D.ad<=22)]
P("   range       n      f_hi     rare var   CpG var   del     unit len")
for nm,s in (("원3-8",b),("원9-10",a),("원11-22",c)):
    P("   %-10s %6s  %.4f   %.3f     %.3f    %.3f   %.1f"
      %(nm,format(len(s),","),s.f_hi.mean(),s.n_var_rare.mean(),s.n_var_cpg.mean(),
        s.n_del.mean(),s.unit_len.mean()))
P()
P("   -> by methylation distal 9-10 sides with distal 3-8 (f_hi 0.03 vs 0.005), but")
P("     by sequence it is closer to distal 11-22 (rare variants, deletions) = the two axes really measure different things.")
P()
D[["sample","slot","ap","ad","N","MZ","SZ","Z"]].to_csv(
    os.path.join(T,"ZONE_INTEGRATED.tsv"),sep="\t",index=False)
P("  saved : 02_tables/ZONE_INTEGRATED.tsv")
LOG.close(); print("done")
