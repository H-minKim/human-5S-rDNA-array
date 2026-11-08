# -*- coding: utf-8 -*-
# 보충 (A)의 CpG 선을 창 폭 1,000 bp 하나로 다시 모은다 2026-10-20
#
#   왜
#     지금 그림의 CpG 선은 창 폭이 자리마다 다르다 — 0~2 kb 는 50 bp · 2~8 kb 는 200 bp ·
#     8 kb 밖은 1,000 bp. 그래서 접합부 쪽만 출렁임이 커 보인다.
#     그것은 서열의 차이가 아니라 해상도의 차이다 — 정본 유닛을 20 벌 이어 붙인
#     완전히 똑같은 가짜 array 를 같은 창으로 재면 출렁임이 50 bp 200 · 200 bp 95 · 1 kb 15 로
#     똑같이 재현된다 (실측 2026-10-20). 밀도 자체(중앙값)는 구간마다 같다.
#   ⇒ 선 하나를 한 가지 창으로 그리면 그 가짜 구조가 사라진다. 패널이 할 말
#     ("array 안이 flank 의 두 배 · 근위와 원위가 거의 대칭")은 1 kb 로도 그대로 전해진다.
#
#   새로 계산하는 값은 없다 — `JUNCTION_SEQ_ANNOT.tsv` 의 bp · n_cpg · n_alu 를 1 kb 로 다시 더할 뿐이다.
#     좌표 u(가장 가까운 온전 유닛 끝에서의 거리)는 `fig1_v5_prep.py` 와 같은 식 을 쓴다.
#   한계 한 줄 : 원본 칸(50·200 bp)은 시작 자리가 드는 1 kb 칸 에 통째로 들어간다.
#     칸 경계에 걸친 것은 쪼개지 않는다 — 지금 표가 쓰는 것과 같은 어림 이다(bp 합은 보존된다).
import io, os, sys
import numpy as np, pandas as pd

# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
DATA_ROOT = os.environ["DATA_ROOT"]
V_ROOT = os.environ["V_ROOT"]
S = DATA_ROOT + "/"
SUMF = V_ROOT + "/03_maps/g8canonv10/CONSENSUS_SUM_G8_full.tsv"
# 표를 읽고 쓰는 자리 — 저장소에서는 표가 코드 옆에 없으므로 환경변수로 받는다 (2026-10-21)
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
W = 1000.0
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

D = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "data_source", "gate"])
D = D[(D.gate == 1) & (D.data_source.astype(str) == "HPRC_ONT_submissions")]
HP = set(zip(D["sample"].astype(str), D.slot.astype(str)))
print("[코호트] 하플로타입 %d" % len(HP))

SUM = pd.read_csv(SUMF, sep="\t", dtype={"slot": str})
SUM["off_원위"] = SUM.off_prox
SUM["off_근위"] = SUM.off_dist
off = {}
for r in SUM.itertuples(index=False):
    off[(r.sample, str(r.slot), "원위")] = float(r.off_원위)
    off[(r.sample, str(r.slot), "근위")] = float(r.off_근위)

A = pd.read_csv(S + "02_tables/JUNCTION_SEQ_ANNOT_V145.tsv", sep="\t", dtype={"slot": str})
#   입력을 V145 로 바꿨다
#     왜 : v143 은 flank fasta 를 파일 이름의 슬롯 으로 붙여 50.1 % 의 하플로타입에
#       다른 하플로의 서열 을 썼다. 원본 FASTA 로 전수(겹치는 하플로타입 457) 대조 :
#       한 칸도 안 다른 하플로타입 — v143 263/457 · V145 457/457 (실측 2026-10-29).
#     그림에 어떻게 나타났나 : 근위 array 안쪽으로 갈수록 CpG 가 오르는 가짜 기울기 —
#       끝(-43~-40 kb) − 가운데(-40~-5 kb) +5.30 → +1.51 CpG/kb (칸 사이 sd 3.15 의 0.5 배)
#       근위 flank 37.77 → 15.85 · 원위 flank 16.16 → 15.30 ⇒ 두 flank 비대칭 2.3 배도 가짜였다
#     분모가 바뀐다 : 하플로타입 470 → 487 (V145 가 더 많이 담고 있다).
A = A[[t in HP for t in zip(A["sample"].astype(str), A.slot.astype(str))]].copy()
o = np.array([off.get((s, sl, sd), np.nan) for s, sl, sd in
              zip(A["sample"].astype(str), A.slot.astype(str), A.side)])
b = A.b.values.astype(float)
A["u"] = np.where(b <= 0, o - np.abs(b), b + o)      # fig1_v5_prep 와 같은 식
A = A.dropna(subset=["u"])
print("[재좌표] 서열 행 %s · 하플로타입 %d" % (format(len(A), ","), A.groupby(["sample", "slot"]).ngroups))

# ---- 창 폭 하나(1 kb)로 다시 담는다 ----
a = np.abs(A.u.values)
lo = (a // W) * W
A["u1k"] = np.where(A.u.values > 0, lo, -lo)
A = A[a <= 52000]
g = (A.groupby(["side", "u1k"])
       .agg(bp=("bp", "sum"), n_cpg=("n_cpg", "sum"), n_alu=("n_alu", "sum"),
            n_hap=("sample", "size")).reset_index())
g["cpg_per_kb"] = 1000.0 * g.n_cpg / g.bp
g["alu_frac"] = g.n_alu / g.bp
g = g.rename(columns={"u1k": "u"}).sort_values(["side", "u"])
OUT = os.path.join(HERE, "FIG3_SUPP_A_SEQ1K.tsv")
g.to_csv(OUT, sep="\t", index=False, float_format="%.5f")

# ---- 스스로 검산 ----
old = pd.read_csv(os.path.join(HERE, "FIG1_V5_A_bp.tsv"), sep="\t")
print()
# 이 검산의 뜻이 바뀌었다 — 입력이 V145 가 되었기 때문이다.
#   맞대는 `FIG1_V5_A_bp.tsv` 는 v143 으로 만든 표 다. 창(1 kb)은 같으므로
#   남는 차이는 오직 입력 표의 차이 이고, 그것이 바로 고친 양 이다.
#   0 이라야 하는 것이 아니다 — 0 이면 오히려 V145 가 안 먹은 것이다.
print("검산 ⑴ — 깊은 쪽(1 kb 창)에서 v143 표와 얼마나 달라졌나")
n_same = n_diff = 0
for side in ("근위", "원위"):
    o2 = old[(old.side == side) & (old.u < -8000) & old.cpg_per_kb.notna()]
    n2 = g[(g.side == side)].set_index("u")
    for r in o2.itertuples(index=False):
        if r.u in n2.index:
            d = abs(float(n2.loc[r.u, "cpg_per_kb"]) - r.cpg_per_kb)
            if d < 0.05:
                n_same += 1
            else:
                n_diff += 1
                if n_diff <= 3:
                    print("   다름 %s u %.0f : 옛 %.2f ↔ 새 %.2f" % (side, r.u, r.cpg_per_kb, n2.loc[r.u, "cpg_per_kb"]))
print("   v143 과 같은 칸 %d · 달라진 칸 %d — 창이 같으므로 이 차이 = 고친 양 이다"
      % (n_same, n_diff))
if n_diff == 0:
    sys.exit("멈춘다 — 한 칸도 안 달라졌다. V145 가 안 읽힌 것이다.")
print()
print("검산 ⑵ — 구간마다 출렁임(이웃 칸 차의 중앙값)")
for side in ("근위", "원위"):
    for lo_, hi_, nm in ((0, 8000, "0~8 kb (전에는 50·200 bp 창)"), (8000, 43000, "8~43 kb (전에도 1 kb 창)")):
        for T, tag in ((old, "옛"), (g, "새")):
            d = T[(T.side == side) & (T.u < 0) & T.cpg_per_kb.notna()].copy()
            d["dist"] = -d.u
            q = d[(d.dist >= lo_) & (d.dist < hi_)].sort_values("dist")
            if len(q) > 2:
                j = np.median(np.abs(np.diff(q.cpg_per_kb.values)))
                print("   %-4s %-28s %s 톱니 %6.2f · 중앙 %6.2f · 칸 %d"
                      % (side, nm, tag, j, q.cpg_per_kb.median(), len(q)))
print()
print("낸 표 %s · 칸 %d (근위 %d · 원위 %d)"
      % (os.path.basename(OUT), len(g), (g.side == "근위").sum(), (g.side == "원위").sum()))
