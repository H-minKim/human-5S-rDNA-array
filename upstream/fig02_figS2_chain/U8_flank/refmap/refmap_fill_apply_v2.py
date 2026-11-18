# -*- coding: utf-8 -*-
# Step 2: fill FASTA_REF_MAP (22 rows)
#   criteria (not relying on names: slots in file names cannot be trusted)
#     ⑴ `fasta_len` == `FLANK_OFFSET.ref_len`  <- the very safeguard that v145 uses
#     ⑵ among the candidates of that length, the one with the individual name in its file name must be unique
#     ⑶ it must differ from the ref_len of the partner slot (if equal it is ambiguous, measured 0 cases)
#   the input map is not overwritten. Output goes under a `_b` name
import io, os, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, "repro_all")
T2 = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/02_tables/"
LOG = io.open(os.path.join(OUTD, "REFMAP_APPLY_V2.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
def m5(p): return hashlib.md5(io.open(p, "rb").read()).hexdigest()

SRC = T2 + "FASTA_REF_MAP.tsv"
OLD = T2 + "FASTA_REF_MAP_2026-10-07.tsv"
DST = T2 + "FASTA_REF_MAP_2026-10-07b.tsv"
M = pd.read_csv(SRC, sep="\t", dtype=str)
e0 = int((M.ref_name.isna() | (M.ref_name == "")).sum())
P("# (2) fill FASTA_REF_MAP (22 rows)")
P()
P("   original `FASTA_REF_MAP.tsv` md5 %s, %d rows, empty %d" % (m5(SRC)[:12], len(M), e0))
P("   FASTA_REF_MAP_2026-10-07.tsv (3 rows) md5 %s" % m5(OLD)[:12])
P()
PLAN = pd.read_csv(os.path.join(OUTD, "REFMAP_FILL_PLAN_V2.tsv"), sep="\t", dtype=str)
P("   plan %d rows" % len(PLAN))
P()
P("   %-13s %-35s %10s  %s" % ("slot", "fasta", "길이", "넣을 ref_name"))
n = 0
for _, r in PLAN.iterrows():
    i = M.index[M.fasta == r.fasta]
    assert len(i) == 1, "fasta appears %d times : %s" % (len(i), r.fasta)
    i = i[0]
    cur = M.at[i, "ref_name"]
    assert not (isinstance(cur, str) and cur.strip()), "already filled %s" % r.fasta
    M.at[i, "ref_name"] = r.ref_name
    M.at[i, "method"] = "길이일치_2026-10-07b"
    n += 1
    P("   %-13s %-35s %10s  %s" % (r["sample"] + "|" + r.slot, r.fasta[:35], r.fasta_len, r.ref_name))
P()
M.to_csv(DST, sep="\t", index=False)
nn = M[M.ref_name.notna() & (M.ref_name != "")]
P("   [new table] FASTA_REF_MAP_2026-10-07b.tsv, md5 %s" % m5(DST)[:12])
P("   filled rows %d, empty rows %d -> %d" % (n, e0, int((M.ref_name.isna() | (M.ref_name == "")).sum())))
P("   ref_name with 2 or more fasta : %d (must be 0)"
  % int((nn.groupby("ref_name").size() > 1).sum()))
P("   original unchanged md5 %s, FASTA_REF_MAP_2026-10-07.tsv unchanged md5 %s" % (m5(SRC)[:12], m5(OLD)[:12]))
LOG.close()
