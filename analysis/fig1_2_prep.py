# -*- coding: utf-8 -*-
# Fig1-2 자료 — 주황선을 검정선과 같은 자(Σ메틸 / Σ호출) 로 다시 낸다 2026-10-08
#   사용자 선택 ⒜ : 두 선을 같은 통계로 맞춘다. 결과(값)는 새 표 로만 내고 옛 표는 그대로 둔다.
#   관문 ① : 옛 표(중앙값·단순평균)를 그대로 재현 한 뒤에만 새 칸을 더한다.
#   같은 여과 : `gate == 1` · `HPRC_ONT_submissions` (fig1_v4_prep 와 같다)
#   가중치 = `n_canon` (그 read-유닛의 CpG 호출 수).
#   정정 — 이 줄에 있던 "검정선의 `sum_p / n_call` 과 같은 뜻이다" 는 틀렸다.
#     가중치는 같지만 평균하는 대상이 다르다 : 여기(유닛 축)는 하드콜 `b_call` 을 평균하고,
#     검정선(`fig1_v5_prep.py` · `FIG1_V5_A_bp.tsv`)은 확률 합 `sum_p` 를 평균한다.
#     가리는 법 : `b_call × n_canon` 은 정수(하드콜) · `sum_p` 는 정수가 아니다(확률 합).
#     크기 : 같은 read×유닛에서 원위 유닛1 하드콜 0.028 ↔ 확률평균 0.149(+0.121) ·
#              근위 유닛10 0.724 ↔ 0.720(−0.004). β 가 0 에 가까운 자리에서만 크게 벌어진다.
#     ⇒ 두 선을 한 그림에 겹치면 원위에서 체계적으로 어긋난다 (실측 차 중앙 +0.103).
#       `Fig 3 (A)` 의 legend 가 이 사실을 적고 있다.
import io, os, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd

HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))   # 
S = os.environ["DATA_ROOT"].rstrip("/") + "/"   # 환경변수로 열었다 — 기본값은 전과 같다
md5 = lambda p: hashlib.md5(io.open(p, "rb").read()).hexdigest()
NED, MIDB = 20, 12

D = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "ap", "ad", "N", "b_call", "gate",
                             "data_source", "n_canon"])
D = D[(D.gate == 1) & (D.data_source.astype(str) == "HPRC_ONT_submissions")].copy()
print("read-유닛 %d · 하플로 %d" % (len(D), D.groupby(["sample", "slot"]).ngroups))

# ---------- 가운데 12 칸 ----------
d = D[D.N >= 2 * NED + MIDB].copy()
mid = d[(d.ap > NED) & (d.ad > NED)].copy()
u = (mid.ap - NED - 1) / np.maximum(mid.N - 2 * NED - 1, 1)
mid["k"] = np.clip((u * MIDB).astype(int), 0, MIDB - 1)
gm = mid.groupby("k").apply(lambda v: pd.Series({
    "med": v.b_call.median(), "mean_simple": v.b_call.mean(),
    "mean_pooled": (v.b_call * v.n_canon).sum() / v.n_canon.sum(),
    "n": len(v), "n_call": int(v.n_canon.sum())}), include_groups=False).reset_index()

OLDM = pd.read_csv(HERE + "/FIG1_V4_A_mid.tsv", sep="\t")
m = OLDM.merge(gm[["k", "med", "mean_simple"]], on="k", suffixes=("_old", "_new"))
d1 = (m.med_old - m.med_new.round(4)).abs().max()
d2 = (m["mean"] - m.mean_simple.round(4)).abs().max()
print()
print("관문 ① 가운데 12 칸 : 중앙값 최대차 %.6f · 단순평균 최대차 %.6f" % (d1, d2))
if max(d1, d2) > 1e-9:
    sys.exit("재현 실패 — 멈춘다")
print("   ⇒ 그대로 재현된다")
gm.to_csv(HERE + "/FIG1_2_MID.tsv", sep="\t", index=False, float_format="%.5f")

# ---------- 끝 20 유닛 ----------
rows = []
for side, col in (("근위", "ap"), ("원위", "ad")):
    g = D[D[col] <= NED]
    for i, v in g.groupby(col):
        w = v.n_canon.astype(float)
        rows.append(dict(side=side, idx=int(i), med=float(v.b_call.median()),
                         mean_simple=float(v.b_call.mean()),
                         mean_pooled=float((v.b_call * w).sum() / w.sum()),
                         n=int(len(v)), n_call=int(w.sum())))
gu = pd.DataFrame(rows)
OLDU = pd.read_csv(HERE + "/FIG1_V4_UNITS.tsv", sep="\t")
mu = OLDU.merge(gu[["side", "idx", "med"]], on=["side", "idx"], suffixes=("_old", "_new"))
d3 = (mu.med_old - mu.med_new.round(4)).abs().max()
print("관문 ① 끝 20 유닛 : 중앙값 최대차 %.6f" % d3)
if d3 > 1e-9:
    sys.exit("재현 실패 — 멈춘다")
print("   ⇒ 그대로 재현된다")
gu.to_csv(HERE + "/FIG1_2_UNITS.tsv", sep="\t", index=False, float_format="%.5f")
print()
print("저장 : FIG1_2_UNITS.tsv (md5 %s) · FIG1_2_MID.tsv (md5 %s)"
      % (md5(HERE + "/FIG1_2_UNITS.tsv")[:12], md5(HERE + "/FIG1_2_MID.tsv")[:12]))
print()
print("가운데 12 칸 : 중앙 %s" % " ".join("%.2f" % v for v in gm.med))
print("             가중평균 %s" % " ".join("%.2f" % v for v in gm.mean_pooled))
