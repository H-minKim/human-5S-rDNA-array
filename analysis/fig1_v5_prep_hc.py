# -*- coding: utf-8 -*-
# Fig 1 (v5) 패널 A 자료 — array 의 시작을 '첫 온전 유닛' 으로 옮긴다
#   지적 : 양 끝의 온전치 못한(부분) 유닛은 array 에 속하지 않는다.
#   지도(G8)에서 부분 유닛 길이 : PARTIAL_PROX = array_pos 1 앞 = 원위 쪽 (중앙 1,635 bp)
#                                  PARTIAL_DIST = 마지막 유닛 뒤 = 근위 쪽 (중앙 1,928 bp)
#     (근거 : L1 의 array_pos == ad 이므로 지도틀 array_pos 1 = 원위 끝 · 이름 PROX/DIST 는 방향 교정 이전 명칭)
#   새 좌표 u : 첫 온전 유닛의 경계에서의 거리(bp) · u>0 바깥(부분 유닛 → flank) · u<0 array 안(온전 유닛)
#   출력 FIG1_V5_A_bp.tsv · FIG1_V5_CHECK.txt
#
#   하드콜 판 (`fig1_v5_prep_hc.py`)
#     원본은 확률 평균(`sum_p / n_call`)으로 창 β 를 냈다. 그런데 같은 그림의 유닛 축은
#     하드콜 비율 이어서 두 선이 다른 자 위에 있었다 (원위에서 체계적으로 +0.103).
#     ⇒ 여기서는 `JUNCTION_FINE_METH_HC.tsv` 의 `sum_h`(ML≥128 호출 수) 를 읽어
#       `beta_hc` 를 함께 낸다. 옛 열(`beta` = 확률 평균)도 그대로 남긴다.
#     출력 : `FIG1_V5_A_bp_HC.tsv` (옛 파일을 덮지 않는다)
import io, os, hashlib
import numpy as np, pandas as pd
# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
DATA_ROOT = os.environ["DATA_ROOT"]
V_ROOT = os.environ["V_ROOT"]
S = DATA_ROOT + "/"
SUMF = V_ROOT + "/03_maps/g8canonv10/CONSENSUS_SUM_G8_full.tsv"
# 표를 읽고 쓰는 자리 — 저장소에서는 표가 코드 옆에 없으므로 환경변수로 받는다 (2026-10-21)
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
LOG = io.open(os.path.join(HERE, "FIG1_V5_PREP_HC.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
UNIT = 2240.0

D = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "data_source", "ap", "ad", "N", "b_call", "gate"])
D = D[(D.gate == 1) & (D.data_source.astype(str) == "HPRC_ONT_submissions")].copy()
HP = set(zip(D["sample"].astype(str), D.slot.astype(str)))
u1_prox = D[D.ap == 1].b_call.median(); u1_dist = D[D.ad == 1].b_call.median()
u2_prox = D[D.ap == 2].b_call.median(); u2_dist = D[D.ad == 2].b_call.median()
P("[유닛 축 기준값(L1)] 근위 첫 온전 유닛 β %.3f (둘째 %.3f) · 원위 첫 온전 유닛 β %.3f (둘째 %.3f)"
  % (u1_prox, u2_prox, u1_dist, u2_dist))

SUM = pd.read_csv(SUMF, sep="\t", dtype={"slot": str})
SUM["off_원위"] = SUM.off_prox          # 지도 시작쪽 부분 유닛 = 원위
SUM["off_근위"] = SUM.off_dist
P("[부분 유닛 길이] 원위 중앙 %d bp (10–90%% %d–%d) · 근위 중앙 %d bp (10–90%% %d–%d)"
  % (SUM.off_원위.median(), SUM.off_원위.quantile(.1), SUM.off_원위.quantile(.9),
     SUM.off_근위.median(), SUM.off_근위.quantile(.1), SUM.off_근위.quantile(.9)))

J = pd.read_csv(S + "02_tables/JUNCTION_FINE_METH_HC.tsv", sep="\t", dtype={"slot": str})
J = J[J["sample"].astype(str) != "HG008-T"].copy()   # flank 도 HPRC 자료만 (공통 규칙 §3 · HG008-T 는 암·PacBio 혼합)
J = J[[t in HP for t in zip(J["sample"].astype(str), J.slot.astype(str))]].copy()
A = pd.read_csv(S + "02_tables/JUNCTION_SEQ_ANNOT.tsv", sep="\t", dtype={"slot": str})
A = A[[t in HP for t in zip(A["sample"].astype(str), A.slot.astype(str))]].copy()
off = {}
for r in SUM.itertuples(index=False):
    off[(r.sample, str(r.slot), "원위")] = float(r.off_원위)
    off[(r.sample, str(r.slot), "근위")] = float(r.off_근위)
def newu(df):
    o = np.array([off.get((s, sl, sd), np.nan) for s, sl, sd in
                  zip(df["sample"].astype(str), df.slot.astype(str), df.side)])
    b = df.b.values.astype(float)
    # b<=0 : array 안 · |b| 가 접합부에서의 거리 → u = o - |b|   (u>0 이면 아직 부분 유닛 안)
    # b>0  : flank                             → u = b + o
    return np.where(b <= 0, o - np.abs(b), b + o)
for T in (J, A):
    T["u"] = newu(T)
J = J.dropna(subset=["u"]); A = A.dropna(subset=["u"])
P("[재좌표] 메틸 행 %d · 서열 행 %d · 하플로 %d" % (len(J), len(A), J.groupby(["sample", "slot"]).ngroups))

def rb(u):
    a = abs(u)
    if a < 2000: w = 50 if u < 0 else 100
    elif a < 10000: w = 200 if u < 0 else 500
    elif a <= 52000: w = 1000 if u < 0 else 5000
    else: return np.nan
    lo = (a // w) * w
    return lo if u > 0 else -lo
J["u2"] = [rb(x) for x in J.u]; A["u2"] = [rb(x) for x in A.u]
J = J.dropna(subset=["u2"]); A = A.dropna(subset=["u2"])
J["beta"] = J.sum_p / J.n_call
J["beta_hc"] = J.sum_h / J.n_call          # 하드콜 비율
g = J.groupby(["side", "u2"]).apply(lambda x: pd.Series(
    {"beta": x.sum_p.sum() / x.n_call.sum(),
     "beta_hc": x.sum_h.sum() / x.n_call.sum(),          # 같은 창을 하드콜 자로
     "q1": x.beta.quantile(.25), "q3": x.beta.quantile(.75),
     "q1_hc": x.beta_hc.quantile(.25), "q3_hc": x.beta_hc.quantile(.75),
     "n_call": x.n_call.sum(), "n_hap": x.groupby(["sample", "slot"]).ngroups})).reset_index()
s = A.groupby(["side", "u2"]).agg(bp=("bp", "sum"), cg=("n_cpg", "sum"), al=("n_alu", "sum")).reset_index()
s["cpg_per_kb"] = 1000 * s.cg / s.bp; s["alu_frac"] = s.al / s.bp
g = g.merge(s[["side", "u2", "cpg_per_kb", "alu_frac"]], on=["side", "u2"], how="left").rename(columns={"u2": "u"})
g.to_csv(os.path.join(HERE, "FIG1_V5_A_bp_HC.tsv"), sep="\t", index=False, float_format="%.5f")
P("[저장] FIG1_V5_A_bp_HC.tsv · 창 %d" % len(g))

# ---- 검산 : 새 좌표에서 첫 온전 유닛 구간의 β 가 유닛 축 값과 맞는가 ----
P()
P("[검산] 새 좌표 u 로 잰 '첫 온전 유닛' 구간 β vs L1 유닛 축 값")
def band(side, lo, hi):
    v = J[(J.side == side) & (-J.u2 >= lo) & (-J.u2 < hi)]      # u<0 = array 안
    return v.sum_p.sum() / v.n_call.sum(), v.n_call.sum()
for side, ref1, ref2 in (("근위", u1_prox, u2_prox), ("원위", u1_dist, u2_dist)):
    b1, n1 = band(side, 0, UNIT); b2, n2 = band(side, UNIT, 2 * UNIT)
    P("   %s : 첫 유닛 구간 β %.3f (L1 %.3f) · 둘째 유닛 구간 β %.3f (L1 %.3f)" % (side, b1, ref1, b2, ref2))
    # 옛 좌표(부분 유닛 포함)로 재면 얼마였나
    v = J[(J.side == side) & (J.b <= 0) & (-J.b < UNIT)]
    P("      (참고) 부분 유닛을 array 로 친 옛 좌표의 같은 구간 β %.3f ← 이래서 근위 첫 유닛이 고메틸로 보였다"
      % (v.sum_p.sum() / v.n_call.sum()))
LOG.close()
