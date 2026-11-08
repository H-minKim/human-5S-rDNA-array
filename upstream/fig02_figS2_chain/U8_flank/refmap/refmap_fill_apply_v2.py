# -*- coding: utf-8 -*-
# ② redo: fill FASTA_REF_MAP (22 rows)
#   why redo: the "3 fillable" of 10-07 was a bug in my own matching code.
#     `norm()` stripped only the `c2seq_`, `hiarr_`, `arr_` prefixes and the `_l0` suffix, so
#     the forms `_f2`, `_f3`, `_f4`, `_f5`, `_f6`, `_sw2`, `_sub250`, `_p1`, `_p2` were counted as 0 candidates.
#     Moreover it required `slot_in_name == slot`, but slots in file names cannot be trusted (v143 trap).
#   fixed criteria (not relying on names)
#     ⑴ `fasta_len` == `FLANK_OFFSET.ref_len`  <- the very safeguard that v145 uses
#     ⑵ among the candidates of that length, the one with the individual name in its file name must be unique
#     ⑶ it must differ from the ref_len of the partner slot (if equal it is ambiguous, measured 0 cases)
#   neither the original nor the 10-07 version (3 rows) is overwritten. Output goes under a `_b` name
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
P("# ② 다시 — FASTA_REF_MAP 메우기 (22 줄)")
P()
P("   원본 `FASTA_REF_MAP.tsv` md5 %s · %d 행 · 빈 %d" % (m5(SRC)[:12], len(M), e0))
P("   10-07 판(3 줄) md5 %s  ← 이 판은 버그 산출. 보존만 하고 안 쓴다" % m5(OLD)[:12])
P()
PLAN = pd.read_csv(os.path.join(OUTD, "REFMAP_FILL_PLAN_V2.tsv"), sep="\t", dtype=str)
P("   계획 %d 줄" % len(PLAN))
P()
P("   %-13s %-35s %10s  %s" % ("슬롯", "fasta", "길이", "넣을 ref_name"))
n = 0
for _, r in PLAN.iterrows():
    i = M.index[M.fasta == r.fasta]
    assert len(i) == 1, "fasta 가 %d 번 : %s" % (len(i), r.fasta)
    i = i[0]
    cur = M.at[i, "ref_name"]
    assert not (isinstance(cur, str) and cur.strip()), "이미 차 있음 %s" % r.fasta
    M.at[i, "ref_name"] = r.ref_name
    M.at[i, "method"] = "길이일치_2026-10-07b"
    n += 1
    P("   %-13s %-35s %10s  %s" % (r["sample"] + "|" + r.slot, r.fasta[:35], r.fasta_len, r.ref_name))
P()
M.to_csv(DST, sep="\t", index=False)
nn = M[M.ref_name.notna() & (M.ref_name != "")]
P("   [새 표] FASTA_REF_MAP_2026-10-07b.tsv · md5 %s" % m5(DST)[:12])
P("   채운 줄 %d · 빈 행 %d → %d" % (n, e0, int((M.ref_name.isna() | (M.ref_name == "")).sum())))
P("   한 ref_name 에 fasta 2 개 이상 : %d (0 이어야 한다)"
  % int((nn.groupby("ref_name").size() > 1).sum()))
P("   원본 불변 md5 %s · 10-07 판 불변 md5 %s" % (m5(SRC)[:12], m5(OLD)[:12]))
LOG.close()
