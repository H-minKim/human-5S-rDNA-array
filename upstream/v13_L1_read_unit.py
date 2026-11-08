# -*- coding: utf-8 -*-
# L1 = read x 유닛 통합표.  요청 1-10 + 추가 (A) read 수준 통제 변수.
#
#  경로는 전부 forward slash (heredoc 백슬래시 뭉개짐 함정 회피)
#  beta 두 정의를 모두 들고 간다 : b_call(문턱 ML>=128) · b_prob(확률 평균, 옛 판/논문 정의)
#  문턱은 실측 골짜기 0.508 (b_prob 기준)
#  극성 : array_pos 1 = 게놈 원위(DJ).  ap = 근위 순번(ap_pj) · ad = 원위 순번(ap_dj)
import io, os
import numpy as np
import pandas as pd
from collections import Counter

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    이 글이 읽는 자료의 뿌리는 README.md 의
    "Running this code" 절에 적혀 있다."""
    import os as _o, sys as _s
    _x = _o.environ.get(_v)
    if not _x:
        _s.exit(
            "environment variable %s is not set; it is required by %s"
            % (_v, _o.path.basename(__file__))
            + "\n(a module stops the script that imports it)."
            + "\nThe deposited copy has no default; see README.md"
            + " section \"Running this code\"."
        )
    return _x


OUT = _need("RDNA_V11_ROOT")
T = OUT + "/02_tables"
LOGD = OUT + "/04_logs"
REF = OUT + "/05_reference"
MAP = os.path.join(_need("RDNA_VERSION_ROOT"), "08_versions", "v11", "03_maps", "g8canonv10", "CONSENSUS_MAP_G8_full.tsv")
POP = _need("RDNA_POP_PANEL")

THR = 0.508
LO = 0.35
HI = 0.65
COMMON = 0.5

LOG = io.open(os.path.join(LOGD, "V13_L1.txt"), "w", encoding="utf8", newline="\n")


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


P("=" * 100)
P("L1 = read x 유닛 통합표")
P("=" * 100)
P()

# ---------------------------------------------------------------- 0. 적격 슬롯
Q = pd.read_csv(os.path.join(T, "SLOT_QC.tsv"), sep="\t")
Q["slot"] = Q.slotS.astype(str)
E = Q[Q.eligible == 1][["sample", "slot", "N", "nread", "gate", "cov_frac",
                        "depth", "n_break", "revcomp"]]
E = E.rename(columns={"gate": "gate_rate", "nread": "n_read_slot", "N": "N_qc"})
P("  적격 슬롯 %d / %d · 개체 %d" % (len(E), len(Q), E["sample"].nunique()))
P("  두 슬롯 다 적격인 개체 %d" % int((E.groupby("sample").size() == 2).sum()))
P()

# ---------------------------------------------------------------- 1. read x 유닛
A = pd.read_parquet(os.path.join(T, "ru_v11.parquet"))
P("  ru_v11 전체 %s 행" % format(len(A), ","))
A = A[A.data_source.astype(str) == "HPRC_ONT_submissions"].copy()
A["slot"] = A.slot.astype(str)
A = A.merge(E, on=["sample", "slot"], how="inner")
P("  적격 슬롯으로 거른 뒤 %s 행 · read %s"
  % (format(len(A), ","), format(A.read_id.nunique(), ",")))
assert (A.ap_pj == A.N - A.ap_dj + 1).all(), "극성 축 불일치"
P("  극성 검산 ap_pj == N-ap_dj+1 : 전건 일치")
P()

A["array_pos"] = A.ap.astype("int32")        # 지도 조인 키 (원본 array_pos)
A["ad"] = A.ap_dj.astype("int32")            # 원위 순번
A["ap"] = A.ap_pj.astype("int32")            # 근위 순번
A = A[(A.ap >= 1) & (A.ad >= 1)].copy()

# --- 요청 2·3 : 자리
A["side"] = np.where(A.ap.values <= A.ad.values, "prox", "dist")
A["ord_near"] = np.minimum(A.ap.values, A.ad.values).astype("int32")
A["pos_label"] = np.where(A.side.values == "prox",
                          ["근위%d" % v for v in A.ord_near.values],
                          ["원위%d" % v for v in A.ord_near.values])
A["rel"] = ((A.ap - 1) / np.maximum(A.N - 1, 1)).astype("float32")

# --- 요청 6-9 : 메틸
A = A.rename(columns={"b": "b_call"})
A["st"] = (A.b_call > THR).astype("int8")          # 하드콜 (문턱 THR 그대로)
A["st3"] = np.where(A.b_call < LO, "L", np.where(A.b_call > HI, "H", "M"))  # 하드콜
A["r_gene_sp"] = np.where(A.b_sp.values > 0,
                          A.b_gene.values / np.maximum(A.b_sp.values, 1e-9), np.nan)
A["d_gene_sp"] = (A.b_gene - A.b_sp).astype("float32")
A["gate"] = ((A.frac_canon >= 0.9) & (A.n_canon >= 150)).astype("int8")

# --- (A) read 안 자리 : 말단 품질 통제
A = A.sort_values(["read_id", "array_pos"], kind="mergesort").reset_index(drop=True)
g = A.groupby("read_id", sort=False)
A["n_unit_read"] = g["array_pos"].transform("size").astype("int32")
A["unit_in_read"] = (g.cumcount() + 1).astype("int32")
A["dist_to_read_end"] = (A.n_unit_read - A.unit_in_read + 1).astype("int32")
A["dist_to_read_edge"] = np.minimum(A.unit_in_read, A.dist_to_read_end).astype("int32")

# --- (A) 블록 : 같은 read 안 연속 자리에서 같은 상태가 이어지는 덩어리
rid = A.read_id.values
apv = A.array_pos.values
stv = A.st.values
newseg = np.ones(len(A), dtype=bool)
newseg[1:] = (rid[1:] != rid[:-1]) | (apv[1:] != apv[:-1] + 1) | (stv[1:] != stv[:-1])
A["blk"] = (np.cumsum(newseg) - 1).astype("int32")
bg = A.groupby("blk", sort=False)
A["blk_len"] = bg["array_pos"].transform("size").astype("int32")
A["blk_idx"] = (bg.cumcount() + 1).astype("int32")
A["blk_state"] = A.st
blk_a0 = bg["array_pos"].transform("min")
blk_a1 = bg["array_pos"].transform("max")
lo_read = bg["unit_in_read"].transform("min") == 1
hi_read = bg["dist_to_read_end"].transform("min") == 1
A["blk_trunc"] = ((lo_read & (blk_a0 != 1)) | (hi_read & (blk_a1 != A.N))).astype("int8")

# --- 요청 10 : 양옆 (자리 번호가 연속일 때만 이웃)
key = A[["read_id", "array_pos", "b_prob", "b_call", "st", "b_gene", "b_sp"]].copy()
for shift, tag in ((1, "prev"), (-1, "next")):
    k = key.copy()
    k["array_pos"] = k.array_pos + shift
    k = k.rename(columns={"b_prob": "b_prob_" + tag, "b_call": "b_call_" + tag,
                          "st": "st_" + tag, "b_gene": "b_gene_" + tag,
                          "b_sp": "b_sp_" + tag})
    A = A.merge(k, on=["read_id", "array_pos"], how="left")
A["agree_prev"] = np.where(A.st_prev.notna(), (A.st == A.st_prev).astype(float), np.nan)
A["agree_next"] = np.where(A.st_next.notna(), (A.st == A.st_next).astype(float), np.nan)
P("  양옆 관측된 비율 : 이전 %.3f · 다음 %.3f"
  % (A.st_prev.notna().mean(), A.st_next.notna().mean()))
P()

# ---------------------------------------------------------------- 2. 지도 : 변이
MP = pd.read_csv(MAP, sep="\t", dtype={"sample": str, "slot": str})
P("  지도 전체 %s 행 · seg_type %s"
  % (format(len(MP), ","), MP.seg_type.astype(str).value_counts().to_dict()))
MP = MP[MP.array_pos > 0].copy()
MP = MP[["sample", "slot", "array_pos", "unit_len", "n_var", "n_var_gene",
         "varlist", "chain_ok", "gap_to_next"]]
MP["array_pos"] = MP.array_pos.astype("int32")
P("  온전 유닛(array_pos>0) %s 행" % format(len(MP), ","))
A = A.merge(MP, on=["sample", "slot", "array_pos"], how="left")
P("  지도가 붙은 read-유닛 %.4f" % A.varlist.notna().mean())
P()

# ---------------------------------------------------------------- 3. 변이 그룹
VG = pd.read_csv(os.path.join(REF, "V11_VARIANT_GROUPS.tsv"), sep="\t")
VG["tok"] = VG.g8_token_example.astype(str)
VG = VG.drop_duplicates("tok")
P("  변이 그룹표 %s 토큰" % format(len(VG), ","))

cnt = Counter()
for s in MP.varlist.dropna().astype(str).values:
    if s == "-" or not s:
        continue
    for t in s.split(","):
        cnt[t] += 1
TOK = pd.DataFrame({"tok": list(cnt.keys()), "n_unit_seen": list(cnt.values())})
TOK = TOK.merge(VG, on="tok", how="left")
TOK["matched"] = TOK.variant.notna().astype(int)
P("  지도 varlist 토큰 %s · 그룹표와 맞은 것 %s (%.4f)"
  % (format(len(TOK), ","), format(int(TOK.matched.sum()), ","), TOK.matched.mean()))
P("  맞은 토큰이 덮는 관측 비율 %.4f"
  % (TOK.loc[TOK.matched == 1, "n_unit_seen"].sum() / TOK.n_unit_seen.sum()))
TOK.sort_values("n_unit_seen", ascending=False).to_csv(
    os.path.join(T, "VAR_GROUP_MAP.tsv"), sep="\t", index=False)
P("  저장 : 02_tables/VAR_GROUP_MAP.tsv")
P()

CL = dict(zip(TOK.tok, TOK["class"].astype(str)))
EL = dict(zip(TOK.tok, TOK["element"].fillna("NA").astype(str)))
CPG = dict(zip(TOK.tok, (TOK.is_cpg.astype(str) == "Y").astype(int)))
COD = dict(zip(TOK.tok, (TOK.is_coding.astype(str) == "Y").astype(int)))
UFQ = dict(zip(TOK.tok, TOK.unit_freq.astype(float)))
ELEMS = ["AluY", "NTS1", "GArich", "NTS2", "coding", "NA"]

U = MP[["sample", "slot", "array_pos", "varlist"]].copy()
rec = []
for s in U.varlist.values:
    if not isinstance(s, str) or s in ("-", ""):
        rec.append((0, 0, 0, 0, 0, 0, 0, 0, [0] * 6, "", 0, np.nan))
        continue
    tk = s.split(",")
    ncom = nrare = nun = 0
    nsnv = nins = ndel = 0
    ncpg = ncod = 0
    ec = [0] * 6
    sig = []
    fmin = 1.0
    for t in tk:
        f = UFQ.get(t, np.nan)
        if not (f == f):
            nun += 1
        elif f >= COMMON:
            ncom += 1
        else:
            nrare += 1
            sig.append(t)
            fmin = min(fmin, f)
        c = CL.get(t, "")
        if c == "snv":
            nsnv += 1
        elif c == "ins":
            nins += 1
        elif c == "del":
            ndel += 1
        ncpg += CPG.get(t, 0)
        ncod += COD.get(t, 0)
        e = EL.get(t, "NA")
        ec[ELEMS.index(e)] += 1 if e in ELEMS else 0
        if e not in ELEMS:
            ec[5] += 1
    rec.append((ncom, nrare, nun, nsnv, nins, ndel, ncpg, ncod, ec,
                ";".join(sorted(sig)), len(sig), (fmin if sig else np.nan)))

U["n_var_common"] = [r[0] for r in rec]
U["n_var_rare"] = [r[1] for r in rec]
U["n_var_unmatched"] = [r[2] for r in rec]
U["n_snv"] = [r[3] for r in rec]
U["n_ins"] = [r[4] for r in rec]
U["n_del"] = [r[5] for r in rec]
U["n_var_cpg"] = [r[6] for r in rec]
U["n_var_coding"] = [r[7] for r in rec]
for i, e in enumerate(ELEMS):
    U["nv_" + e] = [r[8][i] for r in rec]
U["vsig"] = [r[9] for r in rec]
U["vsig_n"] = [r[10] for r in rec]
U["vsig_minfreq"] = [r[11] for r in rec]
U = U.drop(columns=["varlist"])
A = A.merge(U, on=["sample", "slot", "array_pos"], how="left")
P("  유닛 변이 요약 : 흔한 변이 중앙 %.0f · 희귀 변이 중앙 %.0f · 희귀 있는 유닛 %.3f"
  % (U.n_var_common.median(), U.n_var_rare.median(), (U.n_var_rare > 0).mean()))
P("  CpG 관련 변이 있는 유닛 %.3f · 코딩부 변이 있는 유닛 %.3f"
  % ((U.n_var_cpg > 0).mean(), (U.n_var_coding > 0).mean()))
P()

# ---------------------------------------------------------------- 4. 집단
PP = pd.read_csv(POP, sep="\t")[["sample", "pop", "super_pop"]].dropna(subset=["sample"])
PP = PP.drop_duplicates("sample")
A = A.merge(PP, on="sample", how="left")
P("  집단 라벨이 붙은 개체 %d / %d"
  % (A.dropna(subset=["super_pop"])["sample"].nunique(), A["sample"].nunique()))
P()

# ---------------------------------------------------------------- 5. 저장
COLS = (["sample", "slot", "read_id", "data_source", "cfg", "kind", "strand",
         "array_pos", "ap", "ad", "side", "ord_near", "pos_label", "rel", "N",
         "b_call", "b_prob", "st", "st3",
         "b_alu", "b_nts1", "b_ga", "b_nts2", "b_ups", "b_gene", "b_sp",
         "r_gene_sp", "d_gene_sp",
         "n_canon", "frac_canon", "gate", "n_raw",
         "n_alu", "n_nts1", "n_ga", "n_nts2", "n_ups", "n_gene", "n_sp",
         "bA", "bB", "nA", "nB",
         "unit_in_read", "dist_to_read_end", "dist_to_read_edge", "n_unit_read",
         "blk", "blk_len", "blk_idx", "blk_state", "blk_trunc",
         "b_prob_prev", "b_prob_next", "st_prev", "st_next",
         "agree_prev", "agree_next",
         "b_gene_prev", "b_gene_next", "b_sp_prev", "b_sp_next",
         "unit_len", "n_var", "n_var_gene", "varlist",
         "n_var_common", "n_var_rare", "n_var_unmatched",
         "n_snv", "n_ins", "n_del", "n_var_cpg", "n_var_coding"]
        + ["nv_" + e for e in ELEMS]
        + ["vsig", "vsig_n", "vsig_minfreq", "chain_ok", "gap_to_next",
           "N_qc", "n_read_slot", "gate_rate", "cov_frac", "depth", "n_break",
           "revcomp", "pop", "super_pop"])
COLS = [c for c in COLS if c in A.columns]
L1 = A[COLS].copy()
L1.to_parquet(os.path.join(T, "L1_READ_UNIT.parquet"), index=False)
P("  저장 : 02_tables/L1_READ_UNIT.parquet  %s 행 x %d 열"
  % (format(len(L1), ","), L1.shape[1]))
P("     read %s · 슬롯 %d · 개체 %d"
  % (format(L1.read_id.nunique(), ","),
     L1.groupby(["sample", "slot"], observed=True).ngroups, L1["sample"].nunique()))
P()
P("  열 목록")
for i in range(0, len(COLS), 6):
    P("    " + " ".join("%-20s" % c for c in COLS[i:i + 6]))
LOG.close()
print("done")
