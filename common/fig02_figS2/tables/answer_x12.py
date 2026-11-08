# -*- coding: utf-8 -*-
# Answers to questions - reproduce the majority-type path of figx9_12_prep.py exactly, then measure 9,10,11,12
import io, os, sys, collections
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
T   = (__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/02_tables/")
MAJ = (__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/03_figures/latest_2026-10-01_majority_reference/L2_MAJREF_variant_cols_only.parquet")
OUTD = (__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/03_figures/00_variant_position_answer_2026-10-11")
os.makedirs(OUTD, exist_ok=True)
CUT, NTOP, MINU = 0.90, 40, 2

def var_windows(D):
    return [("근1_3",(D.ap<=3)),("근4_7",(D.ap>=4)&(D.ap<=7)),("근8_12",(D.ap>=8)&(D.ap<=12)),
            ("근13_25",(D.ap>=13)&(D.ap<=25)),("중앙",(D.ap>25)&(D.ad>25)),
            ("원13_25",(D.ad>=13)&(D.ad<=25)),("원8_12",(D.ad>=8)&(D.ad<=12)),
            ("원4_7",(D.ad>=4)&(D.ad<=7)),("원1_3",(D.ad<=3))]

def build_var(D, label):
    D = D.copy(); D["sl"] = D["sample"].astype(str)+"_"+D.slot.astype(str)
    slots = sorted(D.sl.unique())
    tok = D.varlist.fillna("-").astype(str).str.split(",")
    cnt = collections.Counter()
    for tl in tok:
        for t in tl:
            if t and t != "-": cnt[t] += 1
    cand = [t for t,c in cnt.items() if 0.02*len(D) <= c <= 0.98*len(D)]
    TOKM = pd.DataFrame({t: tok.map(lambda L,t=t: 1.0 if t in L else 0.0).values for t in cand}, index=D.index)
    sd = TOKM.groupby(D.sl.values).mean().std()
    TOP = sd.sort_values(ascending=False).head(NTOP).index.tolist()
    out = pd.DataFrame(index=slots)
    for wn, m in var_windows(D):
        sub = TOKM[m.values]; key = D.sl.values[m.values]
        g = sub.groupby(key); keep = g.size()[g.size()>=MINU].index
        mu = g.mean().reindex(slots); mu.loc[~mu.index.isin(keep), :] = np.nan
        for t in TOP: out["%s|VAR_%s" % (wn,t)] = mu[t].values
    print("   [%s] 후보 토큰 %d → 쓰는 토큰 %d · 열 %d · 하플로 %d" % (label,len(cand),len(TOP),out.shape[1],out.shape[0]))
    return out, TOP

L = pd.read_parquet(T+"L2_HAP_UNIT.parquet", columns=["sample","slot","array_pos","ap","ad","N","n_read","varlist"])
L = L[L.n_read>=5].copy()
M = pd.read_parquet(MAJ, columns=["sample","slot","array_pos","ap","ad","N","varlist"])
M = M.merge(L[["sample","slot","array_pos","n_read"]], on=["sample","slot","array_pos"], how="inner")
print("MAJREF 자리 %s · 하플로 %d" % (format(len(M),","), (M["sample"].astype(str)+"_"+M.slot.astype(str)).nunique()))
V, TOP = build_var(M, "다수형")

haps = list(V.index)
Xb = (V.fillna(0).values > 0).astype(np.int8)
tot = Xb.sum(1).astype(float)
inter = (Xb.astype(np.int32) @ Xb.astype(np.int32).T).astype(float)
uni = tot[:,None]+tot[None,:]-inter
J = np.where(uni>0, inter/uni, np.nan)
par = list(range(len(haps)))
def find(a):
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
iu = np.triu_indices(len(haps),1)
for a,b in zip(iu[0][J[iu]>=CUT], iu[1][J[iu]>=CUT]):
    ra,rb=find(a),find(b)
    if ra!=rb: par[rb]=ra
lab = np.array([find(i) for i in range(len(haps))])
sz = collections.Counter(lab)
big = [k for k,v in sorted(sz.items(), key=lambda x:-x[1]) if v>=3]
gid = {k:i+1 for i,k in enumerate(big)}
print("\n관문 — 발행 표와 맞는가")
print("   무리(>=3) %d · 무리 안 하플로 %d · 최대 무리 %d" % (len(big), sum(sz[k] for k in big), max(sz[k] for k in big)))
G = pd.read_csv((__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/03_figures/paper_Fig/FIGX9_12_GROUPS.tsv"), sep="\t")
GM = G[G.기준=="다수형"]
print("   발행 FIGX9_12_GROUPS 다수형 : 무리 %d · 하플로합 %d · 최대 %d" % (len(GM), GM.하플로.sum(), GM.하플로.max()))
rep_sizes = sorted([sz[k] for k in big], reverse=True)
pub_sizes = sorted(GM.하플로.tolist(), reverse=True)
print("   무리 크기 수열 일치 : %s" % (rep_sizes==pub_sizes))

# ---------------- Q9 : within-cluster pairwise Jaccard (chaining diagnostic)
rows=[]
for k in big:
    ii=[i for i in range(len(haps)) if lab[i]==k]
    s=J[np.ix_(ii,ii)]; pv=s[np.triu_indices(len(ii),1)]
    rows.append(dict(무리=gid[k], 하플로=len(ii), 쌍=len(pv),
                     자카드최소=round(float(pv.min()),4), 자카드중앙=round(float(np.median(pv)),4),
                     자카드평균=round(float(pv.mean()),4),
                     문턱이상_몫=round(float((pv>=CUT).mean()),4)))
CH=pd.DataFrame(rows)
CH.to_csv(os.path.join(OUTD,"ans9_within_cluster_pairwise_jaccard.tsv"), sep="\t", index=False)
print("\n[9] 무리 안 쌍별 자카드 (단일연결 사슬 진단)")
print(CH.head(13).to_string(index=False))

# ---------------- Q11 : member list
mem=[]
for i,h in enumerate(haps):
    mem.append(dict(무리=gid.get(lab[i],0), hap=h, sample=h.rsplit("_",1)[0], slot=h.rsplit("_",1)[1],
                    무리크기=int(sz[lab[i]]), 보유토큰수=int(Xb[i].sum())))
MEMD=pd.DataFrame(mem).sort_values(["group","hap"])
MEMD.to_csv(os.path.join(OUTD,"ans11_cluster_haplotype_list.tsv"), sep="\t", index=False)
print("\n[11] 명단 저장 — 무리에 든 하플로 %d · 무리 0(안 든 것) %d"
      % (int((MEMD.무리>0).sum()), int((MEMD.무리==0).sum())))

# ---------------- Q12 : two slots of the same individual
smp=pd.Series([h.rsplit("_",1)[0] for h in haps])
cc=smp.value_counts()
two=[s for s,c in cc.items() if c==2]
bothin=same=0
for s in two:
    ii=[i for i in range(len(haps)) if haps[i].rsplit("_",1)[0]==s]
    if all(lab[i] in big for i in ii):
        bothin+=1
        if lab[ii[0]]==lab[ii[1]]: same+=1
print("\n[12] 개체 %d (두 슬롯 %d · 한 슬롯 %d)" % (cc.size,len(two),int((cc==1).sum())))
print("   두 슬롯이 둘 다 무리에 든 개체 %d · 그 중 같은 무리 %d (%.0f%%)"
      % (bothin, same, 100*same/max(bothin,1)))

# ---------------- Q10 : distinguishing tokens
VF=V.fillna(0); hasz=(VF>0)
drows=[]
inb=np.array([lab[i] in big for i in range(len(haps))])
for k in big:
    m=np.array([lab[i]==k for i in range(len(haps))])
    o=~m
    for c in VF.columns:
        v=hasz[c].values
        pin=v[m].mean(); pout=v[o].mean()
        if pin>=0.80 and (pin-pout)>=0.30:
            z,t=c.split("|VAR_")
            drows.append(dict(무리=gid[k], 구역=z, 토큰=t,
                              무리안_보유율=round(float(pin),4), 무리밖_보유율=round(float(pout),4),
                              차=round(float(pin-pout),4)))
DT=pd.DataFrame(drows).sort_values(["group","difference"],ascending=[True,False])
DT.to_csv(os.path.join(OUTD,"ans10_cluster_distinguishing_tokens.tsv"), sep="\t", index=False)
print("\n[10] 가르는 토큰 (무리안>=0.80 & 차>=0.30) %d 줄 · 무리 %d / %d"
      % (len(DT), DT.무리.nunique(), len(big)))
print(DT.groupby("group").size().head(13).to_string())
io.open(os.path.join(OUTD,"ans_tokens_used_40.txt"),"w",encoding="utf8").write("\n".join(TOP))
print("\n[저장] %s" % OUTD)
