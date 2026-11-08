#!/bin/bash
# Transcribed from the run record: the code ran as a shell heredoc or a one-line command, so no script file existed on disk.
# Writes the cohort tier list TIERS_H0H3_SLOTS_2026-10-05.tsv (2,161 rows: H0 604, H1 557, H2 510, H3 490).
# Machine paths are given by RDNA_DATA_ROOT / RDNA_WORK_DIR (no default); the computing lines are unchanged.
# step 1 (record line 19136): copy the five slot lists and the dataset ledger
D="$RDNA_WORK_DIR/coh"; mkdir -p "$D" && cp "$RDNA_DATA_ROOT/for_paper/00_handover/handover_v11_data_ledger_2026-10-04/COHORT_SLOTS_ledger_20261005.tsv" "$D/ledger.tsv" && cp "$RDNA_DATA_ROOT/for_paper/00_paper_figures_distal_layout/tables/COHORT_SLOTS_variant_position_20261005.tsv" "$D/varpos.tsv" && cp "$RDNA_DATA_ROOT/v11_readvar/00_handover_control/COHORT_SLOTS_5S_methylation_20261005.tsv" "$D/meth5s.tsv" && cp "$RDNA_DATA_ROOT/5S_methylation_model/cmd_20261005/out/COHORT_SLOTS_haplotype_effect_20261005.tsv" "$D/haploinf.tsv" && cp "$RDNA_DATA_ROOT/5S_array_methylation_handover/02_output_tables/cohort_list_2026-10-05/COHORT_SLOTS_methylation_blocks_20261005.tsv" "$D/methblock.tsv" && cp "$RDNA_DATA_ROOT/for_paper/00_handover/handover_v11_data_ledger_2026-10-04/V11F4_DATASET_hap_flags.tsv" "$D/flags.tsv" && ls -la "$D" && md5sum "$D"/*.tsv

# step 2 (record line 19141): load every cohort as a set of sample|slot -> coh.pkl
cd "$RDNA_WORK_DIR/coh" && PYTHONIOENCODING=utf-8 python - <<'PYEOF'
import pandas as pd, pickle
F={"대장":"ledger.tsv","변이위치":"varpos.tsv","5S메틸":"meth5s.tsv","하플로영향":"haploinf.tsv","메틸블록":"methblock.tsv"}
C={}
for s,p in F.items():
    t=pd.read_csv(p,sep="\t",dtype=str).fillna(""); t.columns=[c.strip() for c in t.columns]
    n=0
    for c,g in t.groupby("cohort"):
        g=g[g["slot"].isin(["1","2"])]
        if len(g)==0: continue
        C[(s,c)]=set(g["sample"]+"|"+g["slot"]); n+=1
    print(f"{s:8s} 코호트 {n:2d} · 행 {len(t):5d} · 열 {list(t.columns)}")
pickle.dump(C,open("coh.pkl","wb"))
print(f"\n적재 {len(C)} 코호트")
for (s,c),v in sorted(C.items(), key=lambda x:(x[0][0],-len(x[1]))):
    ns=len({k.split('|')[0] for k in v})
    print(f"  {s:8s} {c:36s} {len(v):4d} / {ns:4d}")
PYEOF

# step 3 (record line 19728): H0 = HPRC, H1 = variant-position list cohort primary_564 & HPRC & not 1kGP/GIAB ONT, H2 = ledger cohort E_QC-eligible_510 (SLOT_QC eligible), H3 = H1 & H2
cd "$RDNA_WORK_DIR/coh" && PYTHONIOENCODING=utf-8 python - <<'PYEOF'
import pickle, pandas as pd
C=pickle.load(open("coh.pkl","rb"))
f=pd.read_csv("flags.tsv",sep="\t",dtype=str).fillna(""); f["key"]=f["sample"]+"|"+f["slot"]
K=lambda m:set(f.loc[m,"key"]); PRIM=C[("변이위치","primary_564")]; ELI=C[("대장","E_QC적격_510")]
hc=f["cohort"].eq("HPRC"); bad=f["ont_data_source"].isin(["1kGP_ONT","GIAB_ONT","GIAB_SOMATIC_ONT"])
H={"H0":K(hc),"H1":K(f["key"].isin(PRIM)&hc&~bad),"H2":ELI}; H["H3"]=H["H1"]&H["H2"]
idx=f.set_index("key")
rows=[{"tier":t,"sample":k.split("|")[0],"slot":k.split("|")[1],
       "N_units":idx.at[k,"N_units"],"ont_data_source":idx.at[k,"ont_data_source"],
       "canon_source":idx.at[k,"canon_source"],"qc_eligible":idx.at[k,"qc_eligible"]}
      for t in ["H0","H1","H2","H3"] for k in sorted(H[t])]
d=pd.DataFrame(rows); d.to_csv("TIERS.tsv",sep="\t",index=False,encoding="utf-8",lineterminator="\n")
print(f"행 {len(d)}")
for t in ["H0","H1","H2","H3"]:
    g=d[d["tier"].eq(t)]; p=g.groupby("sample")["slot"].nunique()
    print(f"   {t} 슬롯 {len(g):4d} · 개체 {g['sample'].nunique():4d} · 짝 {int((p==2).sum()):4d}")
PYEOF
cp TIERS.tsv "$RDNA_DATA_ROOT/rDNA_registry/registry/TIERS_H0H3_SLOTS_2026-10-05.tsv" && md5sum TIERS.tsv && cd "$RDNA_DATA_ROOT/rDNA_registry/registry" && md5sum TIERS_H0H3_SLOTS_2026-10-05.tsv
