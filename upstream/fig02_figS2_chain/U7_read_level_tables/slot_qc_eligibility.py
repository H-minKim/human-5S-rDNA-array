# -*- coding: utf-8 -*-
# Transcribed from the run record: the code ran as a shell heredoc or a one-line command, so no script file existed on disk.
# Writes SLOT_QC.tsv in the 02 tables folder (590 HPRC-ONT slots, eligible 510), read by v13_L1_read_unit.py.
# Eligible = N>=30, no map break, covered positions/N >= 0.8, reads per position >= 5, gate pass >= 0.85, reads >= 10.
# Machine paths are given by RDNA_DATA_ROOT / RDNA_WORK_DIR (no default); the computing lines are unchanged.
import numpy as np, pandas as pd, os, io
OUT=__import__("os").environ["RDNA_DATA_ROOT"]+"/v11_readvar"; V=__import__("os").environ["RDNA_DATA_ROOT"]+"/for_paper/08_versions/v11"
A=pd.read_parquet(os.path.join(OUT,"02_tables","ru_v11.parquet"))
A=A[A.data_source.astype(str)=="HPRC_ONT_submissions"].copy()
A["slotS"]=A.slot.astype(str)
A["gate"]=((A.frac_canon>=0.9)&(A.n_canon>=150)).astype(int)
S=pd.read_csv(os.path.join(V,"03_maps","g8canonv10","CONSENSUS_SUM_G8_full.tsv"),sep="\t")
S["slotS"]=S.slot.astype(str)
g=A.groupby(["sample","slotS"],observed=True)
H=g.agg(N=("N","first"),nru=("b","size"),nread=("read_id","nunique"),
        gate=("gate","mean"),npos=("ap_dj","nunique")).reset_index()
H=H.merge(S[["sample","slotS","N_units","n_break","revcomp"]],on=["sample","slotS"],how="left")
H["cov_frac"]=H.npos/H.N.clip(lower=1)
H["depth"]=H.nru/H.npos.clip(lower=1)
print("슬롯 %d"%len(H))
print()
for c,lab in (("N","array 길이 N"),("nread","read 수"),("npos","덮인 자리 수"),
              ("cov_frac","자리 피복 비율"),("depth","자리당 read 깊이"),
              ("gate","관문 통과율"),("n_break","지도 끊김")):
    v=H[c].astype(float)
    print("  %-16s 중앙 %8.3f · 사분위 %8.3f-%8.3f · 최소 %8.3f · 최대 %8.3f"%(lab,v.median(),v.quantile(.25),v.quantile(.75),v.min(),v.max()))
print()
print("  후보 관문별로 남는 슬롯")
for nm,m in (("전체",np.ones(len(H),bool)),
             ("N>=30",H.N>=30),
             ("끊김 0",H.n_break.fillna(9)==0),
             ("자리 피복>=0.8",H.cov_frac>=0.8),
             ("자리당 깊이>=5",H.depth>=5),
             ("관문 통과율>=0.85",H.gate>=0.85),
             ("read>=10",H.nread>=10)):
    print("    %-18s %3d"%(nm,int(np.sum(m))))
ok=(H.N>=30)&(H.n_break.fillna(9)==0)&(H.cov_frac>=0.8)&(H.depth>=5)&(H.gate>=0.85)&(H.nread>=10)
print()
print("  전부 만족 : %d 슬롯 (개체 %d)"%(int(ok.sum()),H[ok]["sample"].nunique()))
print("     그중 같은 개체 두 슬롯 다 통과 : %d 개체"%int((H[ok].groupby("sample").size()==2).sum()))
H.assign(eligible=ok.astype(int)).to_csv(os.path.join(OUT,"02_tables","SLOT_QC.tsv"),sep="\t",index=False,float_format="%.4f")
print("  저장: 02_tables\\SLOT_QC.tsv")
