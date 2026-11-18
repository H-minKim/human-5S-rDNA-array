# -*- coding: utf-8 -*-
# Step 1 of FASTA_REF_MAP_2026-10-07b.tsv: picks the 22 empty rows to fill (length match + sample name in file name, unique)
# Machine paths are given by RDNA_DATA_ROOT / RDNA_WORK_DIR (no default).
import sys, pandas as pd, re
sys.stdout.reconfigure(encoding='utf-8',errors='replace')
S=__import__('os').environ['RDNA_DATA_ROOT']+'/v11_readvar/'; FIG=S+'03_figures/paper_Fig/'
M=pd.read_csv(S+'02_tables/FASTA_REF_MAP.tsv',sep='\t',dtype=str)
M['fasta_len']=pd.to_numeric(M.fasta_len,errors='coerce')
OFF=pd.read_csv(__import__('os').environ['RDNA_DATA_ROOT']+'/for_paper/08_versions/v11/13_junction50k/mapping_table/FLANK_OFFSET.tsv',sep='\t',dtype=str)
OFF['ref_len']=pd.to_numeric(OFF.ref_len,errors='coerce')
OFF['key']=list(zip(OFF['sample'].astype(str),OFF.slot.astype(str)))
C504=pd.read_csv(FIG+'COHORT_FIGX5_BETA500_SLOTS.tsv',sep='\t',dtype=str)
C478=pd.read_csv(FIG+'COHORT_FIGX5_COMMON478.tsv',sep='\t',dtype=str)
miss=sorted(set(zip(C504['sample'],C504.slot))-set(zip(C478['sample'],C478.slot)))
nn=M[M.ref_name.notna()&(M.ref_name!='')]; taken=set(nn.ref_name)
E=M[M.ref_name.isna()|(M.ref_name=='')].copy()
print('missing slots %d, empty rows in the mapping table %d'%(len(miss),len(E)))
print()
print('%-13s %9s  %-34s %9s %-7s %s'%('slot','목표len','길이일치 후보(빈행)','len','name?','판정'))
ok=[];amb=[];none=[];skip=[]
for sm,sl in miss:
    r=OFF[OFF.key==(sm,sl)]
    if not len(r): skip.append((sm,sl,'FLANK_OFFSET 없음')); continue
    rn=r.iloc[0].ref_name; tgt=r.iloc[0].ref_len
    if rn in taken: skip.append((sm,sl,'이미 매핑됨')); continue
    hit=E[E.fasta_len==tgt]
    nm=hit[hit.fasta.str.contains(sm,regex=False)]
    if len(nm)==1:
        ok.append((sm,sl,rn,nm.iloc[0].fasta,int(tgt)))
        print('%-13s %9.0f  %-34s %9.0f %-7s match'%(sm+'|'+sl,tgt,nm.iloc[0].fasta[:34],tgt,'예'))
    elif len(nm)>1:
        amb.append((sm,sl,list(nm.fasta)))
        print('%-13s %9.0f  %-34s %9.0f %-7s candidates %d'%(sm+'|'+sl,tgt,'multiple',tgt,'예',len(nm)))
    else:
        none.append((sm,sl,int(tgt),len(hit)))
        print('%-13s %9.0f  %-34s %9s %-7s none (length-only matches %d)'%(sm+'|'+sl,tgt,'—','—','—',len(hit)))
print()
print('unique candidate matching both name and length : %d / ambiguous %d / none %d / skipped %d'%(len(ok),len(amb),len(none),len(skip)))
for a in skip: print('   skipped',a)
import collections
print()
print('duplicate length check : among the 40 empty rows, lengths that occur twice or more =',[k for k,v in collections.Counter(E.fasta_len.dropna()).items() if v>1])
pd.DataFrame(ok,columns=['sample','slot','ref_name','fasta','fasta_len']).to_csv('repro_all/REFMAP_FILL_PLAN_V2.tsv',sep='\t',index=False)
print('[saved] repro_all/REFMAP_FILL_PLAN_V2.tsv')
