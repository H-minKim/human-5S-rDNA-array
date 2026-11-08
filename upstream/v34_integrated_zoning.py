# -*- coding: utf-8 -*-
# Compares three zonings and builds and checks a combined grid.
#   (A)  메틸 6구간 : P1 근1-3 · P2 근4-10 · P3 근11-25(무름) · INT · D2 원11-22 · D1 원1-10
#   (B)  서열 구역 : 근1-3 · 근4-17 · 내부(근18+,원9+) · 원1-8
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

# ---- 근위 제3경계 : 25 인가 17 인가 ( 이 '무름' 이라 한 그것)
P("[1] the third proximal boundary : 25.5 (soft), 17, and 17|18 from methylation here")
P()
P("   길이층      근위 f_hi :  근14   근16   근18   근20   근25   근30   내부")
for lo,hi in ((30,59),(60,99),(100,139),(140,300)):
    s=D[(D.N>=lo)&(D.N<=hi)&near_p]
    v=[s.loc[s.ap==k,"f_hi"].mean() for k in (14,16,18,20,25,30)]
    it=D[(D.N>=lo)&(D.N<=hi)&(D.ap>30)&(D.ad>30)]["f_hi"].mean()
    P("   N %3d-%3d              "%(lo,hi)+"  ".join("%.3f"%x for x in v)+"   %.3f"%it)
P()
P("   내부값에 도달하는 자리(내부 대비 98% 이상) — 길이층마다")
for lo,hi in ((30,59),(60,99),(100,139),(140,300)):
    s=D[(D.N>=lo)&(D.N<=hi)&near_p]
    it=D[(D.N>=lo)&(D.N<=hi)&(D.ap>30)&(D.ad>30)]["f_hi"].mean()
    hit=None
    for k in range(4,31):
        v=s.loc[s.ap==k,"f_hi"].mean()
        if v==v and it==it and v>=0.98*it: hit=k; break
    P("   N %3d-%3d : 근위 %s 에서 내부값 도달 (내부 %.3f)"%(lo,hi,hit,it))
P()

# ---- 종합 격자
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
P("[2] 메틸 구역 x 서열 구역 교차표 (자리-하플로타입 수)")
P()
ct=pd.crosstab(D.MZ,D.SZ)
P(ct.to_string())
P()
P("  ⇒ 채워진 칸 %d 개 = 종합 구역 후보"%int((ct.values>0).sum()))
P()
FM=["f_hi","b_lo_mean","p_switch","blk_len_mean","r_lo_gene_sp"]
FS=["n_var_rare","n_var_cpg","n_del","unit_len"]
P("[3] 종합 구역별 성질 (n>=300 인 칸만)")
P()
P("  구역                          n       "+" ".join("%-9s"%c for c in FM)+" | "+" ".join("%-9s"%c for c in FS))
keep=[]
for z,s in D.groupby("Z"):
    if len(s)<300: continue
    keep.append(z)
    P("  %-28s %6s  "%(z,format(len(s),","))
      +" ".join("%-9.4f"%s[c].mean() for c in FM)+" | "+" ".join("%-9.3f"%s[c].mean() for c in FS))
P()
P("[4] 두 축이 독립인 정보를 주는가 — 메틸 축만 / 서열 축만 / 둘 다")
P()
for tgt,lab in (("f_hi","메틸 표적 f_hi"),("n_var_rare","서열 표적 희귀변이 수")):
    y=D[tgt].values.astype(float); ok=np.isfinite(y)
    tot=((y[ok]-y[ok].mean())**2).sum()
    def ev(col):
        m=D[col].values
        pred=pd.Series(y[ok]).groupby(pd.Series(m[ok]).values).transform("mean").values
        return 1-((y[ok]-pred)**2).sum()/tot
    P("   %-22s 메틸축 %.4f · 서열축 %.4f · 둘 다 %.4f  (추가 %+.4f / %+.4f)"
      %(lab,ev("MZ"),ev("SZ"),ev("Z"),ev("Z")-ev("MZ"),ev("Z")-ev("SZ")))
P()
P("[5] 경계가 어긋나서 생긴 칸 — 원위 9-10 (메틸은 '꺼진 도메인', 서열은 '내부')")
P()
a=D[(~near_p)&(D.ad>=9)&(D.ad<=10)]; b=D[(~near_p)&(D.ad>=3)&(D.ad<=8)]
c=D[(~near_p)&(D.ad>=11)&(D.ad<=22)]
P("   구간        n      f_hi     희귀변이   CpG변이   결실    유닛길이")
for nm,s in (("원3-8",b),("원9-10",a),("원11-22",c)):
    P("   %-10s %6s  %.4f   %.3f     %.3f    %.3f   %.1f"
      %(nm,format(len(s),","),s.f_hi.mean(),s.n_var_rare.mean(),s.n_var_cpg.mean(),
        s.n_del.mean(),s.unit_len.mean()))
P()
P("   ⇒ 원9-10 은 메틸로는 원3-8 과 같은 편(f_hi 0.03 vs 0.005)인데")
P("     서열로는 원11-22 쪽(희귀변이·결실)에 가깝다 = 두 축이 진짜로 다른 것을 잰다.")
P()
D[["sample","slot","ap","ad","N","MZ","SZ","Z"]].to_csv(
    os.path.join(T,"ZONE_INTEGRATED.tsv"),sep="\t",index=False)
P("  저장 : 02_tables/ZONE_INTEGRATED.tsv")
LOG.close(); print("done")
