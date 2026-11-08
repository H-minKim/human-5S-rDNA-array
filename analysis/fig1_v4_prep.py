# -*- coding: utf-8 -*-
# Fig 1 (v4) 패널 A 자료 — HPRC ONT 만 · flank 와 array 를 하나의 bp 축에 끊김 없이
#  입력 : 02_tables/JUNCTION_FINE_METH.tsv (v142 · 하플로 x 쪽 x bp 창) · 02_tables/L1_READ_UNIT.parquet(가운데 · 유닛 축)
#         02_tables/JUNCTION_SEQ_ANNOT.tsv (v143 · 같은 창의 CpG 밀도 · Alu 유사 비율)
#  출력 : FIG1_V4_A_bp.tsv (쪽 x 창 x 메틸/서열) · FIG1_V4_A_mid.tsv (가운데 정규화 12 칸) · FIG1_V4_UNITS.tsv
import io, os, hashlib
import numpy as np, pandas as pd
S = os.environ["DATA_ROOT"].rstrip("/") + "/"   # 환경변수로 열었다 — 기본값은 전과 같다
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))   # 
LOG = io.open(os.path.join(HERE, "FIG1_V4_PREP.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
INP = [S + "02_tables/JUNCTION_FINE_METH.tsv", S + "02_tables/L1_READ_UNIT.parquet", S + "02_tables/JUNCTION_SEQ_ANNOT.tsv"]
pd.DataFrame([(p, md5(p), os.path.getsize(p)) for p in INP], columns=["path", "md5", "bytes"]).to_csv(
    os.path.join(HERE, "FIG1_V4_INPUTS.tsv"), sep="\t", index=False)
P("[입력]")
for p in INP: P("   %-30s %s %12d" % (os.path.basename(p), md5(p), os.path.getsize(p)))

D = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "read_id", "data_source", "ap", "ad", "N", "b_call", "gate"])
D = D[(D.gate == 1) & (D.data_source.astype(str) == "HPRC_ONT_submissions")].copy()
D["hap"] = D["sample"] + "|" + D.slot.astype(str)
HP = set(zip(D["sample"].astype(str), D.slot.astype(str)))
P()
P("[HPRC] 하플로 %d · 유닛 %d" % (len(HP), len(D)))

J = pd.read_csv(S + "02_tables/JUNCTION_FINE_METH.tsv", sep="\t", dtype={"slot": str})
J = J[J["sample"].astype(str) != "HG008-T"].copy()   # flank 도 HPRC 자료만 (공통 규칙 §3 · HG008-T 는 암·PacBio 혼합)
J["key"] = list(zip(J["sample"].astype(str), J.slot.astype(str)))
J = J[J.key.isin(HP)]
P("[bp 창] HPRC 하플로 %d · 행 %d" % (J.groupby(["sample", "slot"]).ngroups, len(J)))
J["beta"] = J.sum_p / J.n_call
agg = J.groupby(["side", "b"]).apply(
    lambda x: pd.Series({"beta": x.sum_p.sum() / x.n_call.sum(),
                         "q1": x.beta.quantile(.25), "q3": x.beta.quantile(.75),
                         "n_call": x.n_call.sum(), "n_hap": x.groupby(["sample", "slot"]).ngroups})).reset_index()
try:
    A = pd.read_csv(S + "02_tables/JUNCTION_SEQ_ANNOT.tsv", sep="\t", dtype={"slot": str})
    A["key"] = list(zip(A["sample"].astype(str), A.slot.astype(str)))
    A = A[A.key.isin(HP)]
    sq = A.groupby(["side", "b"]).agg(bp=("bp", "sum"), cg=("n_cpg", "sum"), al=("n_alu", "sum")).reset_index()
    sq["cpg_per_kb"] = 1000 * sq.cg / sq.bp
    sq["alu_frac"] = sq.al / sq.bp
    agg = agg.merge(sq[["side", "b", "cpg_per_kb", "alu_frac", "bp"]], on=["side", "b"], how="left")
    P("[서열] 주석 붙음 · 하플로 %d" % A.groupby(["sample", "slot"]).ngroups)
except FileNotFoundError:
    P("[서열] 아직 없음 (v143 실행 중)")
agg.to_csv(os.path.join(HERE, "FIG1_V4_A_bp.tsv"), sep="\t", index=False, float_format="%.5f")
P("   창 %d (하플로 50 개 이상인 창 %d)" % (len(agg), (agg.n_hap >= 50).sum()))

# 가운데 — 유닛 축(길이 정규화 12 칸 · 양 끝 20 유닛 밖)
NED, MIDB = 20, 12
d = D[D.N >= 2 * NED + MIDB].copy()
m = (d.ap > NED) & (d.ad > NED)
mid = d[m].copy()
u = (mid.ap - NED - 1) / np.maximum(mid.N - 2 * NED - 1, 1)
mid["k"] = np.clip((u * MIDB).astype(int), 0, MIDB - 1)
g = mid.groupby("k").b_call.agg(["median", "mean", "size", lambda v: v.quantile(.25), lambda v: v.quantile(.75)])
g.columns = ["med", "mean", "n", "q1", "q3"]
g.reset_index().to_csv(os.path.join(HERE, "FIG1_V4_A_mid.tsv"), sep="\t", index=False, float_format="%.4f")
P("[가운데] N>=%d 하플로의 안쪽 유닛 %d · 12 칸 중앙 β %s"
  % (2 * NED + MIDB, len(mid), " ".join("%.2f" % v for v in g["med"])))

# 유닛 축 값(끝 20 유닛) — bp 축 위에 유닛 눈금을 얹기 위한 참고값
un = []
for side, col in (("근위", "ap"), ("원위", "ad")):
    v = D[D[col] <= 20].groupby(col).b_call.agg(["median", "size"])
    for i, r in v.iterrows(): un.append((side, int(i), r["median"], int(r["size"])))
pd.DataFrame(un, columns=["side", "idx", "med", "n"]).to_csv(
    os.path.join(HERE, "FIG1_V4_UNITS.tsv"), sep="\t", index=False, float_format="%.4f")
P("[유닛] 끝 20 유닛 중앙값 저장")
LOG.close()
