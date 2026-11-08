# -*- coding: utf-8 -*-
# Supp3_2 의 원자료 표를 낸다 2026-10-10
#   까닭 : `supp3_2_draw.py` 는 중간표 없이 `L1_READ_UNIT` 에서 바로 그린다.
#          ⇒ 원자료 xlsx 에 넣을 표가 없다 — 그림이 보이는 것을 그대로 표로 낸다.
#   그림을 다시 그리지 않는다 — 같은 여과·같은 정렬로 표만 낸다.
#   검산 : 분자 수 52 · 덮인 칸 1,200 · 자리 83 (그림 출력과 같아야 한다)
import io, os, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
DATA_ROOT = os.environ["DATA_ROOT"]

S = DATA_ROOT + "/"
# 표를 읽고 쓰는 자리 — 저장소에서는 표가 코드 옆에 없으므로 환경변수로 받는다 (2026-10-21)
# 예치표의 열 이름과 값을 영문으로 낸다
#   쓰기 바로 전 한 자리 에서만 바꾼다 — 이 글 안의 내부 이름은 그대로 둔다.
_KO2EN_COL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": {"구역": "zone", "유닛안자리": "pos_in_unit",
                                    "유전자기준": "pos_gene_ref", "호출수": "n_call",
                                    "메틸호출수": "n_meth_call"},
    "FIG3_B_POSBETA.tsv": {"층": "stratum", "축": "axis", "자리": "pos"},
    "FIG3_C_LOWCOUNT.tsv": {"컷": "cut", "층": "stratum", "N중앙": "N_median",
                            "저메틸수": "n_low_meth", "저메틸수_sd": "n_low_meth_sd",
                            "비례했다면": "expected_if_proportional",
                            "저메틸몫": "low_meth_fraction"},
    "SUPP3_2_POSSD.tsv": {"분자수": "n_molecules", "평균_b_call": "mean_b_call",
                          "분자사이_sd": "sd_between_molecules"},
    "SUPP3_B_HAP3.tsv": {"하플로타입": "haplotype", "길이층": "length_class",
                         "유닛수": "n_units", "분자순번": "molecule_index",
                         "유닛자리": "unit_position"},
}
_KO2EN_VAL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": ("구역", {"근위": "proximal", "내부": "interior",
                                            "원위": "distal", "전체": "all"}),
    "FIG3_B_POSBETA.tsv": ("층", {"전체": "all"}),
    "FIG3_C_LOWCOUNT.tsv": ("컷", {"0.508(여기서 잴 골)": "0.508 (trough measured here)"}),
    "SUPP3_B_HAP3.tsv": ("길이층", {"짧은 array": "short array", "중간 길이": "medium",
                                   "긴 array": "long array"}),
}


def _to_english(n, df):
    b = os.path.basename(n)
    v = _KO2EN_VAL.get(b)
    if v and v[0] in df.columns:
        df = df.copy()
        df[v[0]] = df[v[0]].map(lambda x: v[1].get(x, x))
    c = _KO2EN_COL.get(b)
    if c:
        df = df.rename(columns=c)
    return df

HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
HAPS, HAPSLOT = "NA19338", 2

L = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "read_id", "ap", "ad", "N", "gate", "b_call"])
L = L[(L.gate == 1) & (L["sample"].astype(str) == HAPS) & (L.slot.astype(int) == HAPSLOT)].copy()
NN = int(L.N.iloc[0])
order = list(L.groupby("read_id").b_call.mean().sort_values().index)
pos = {r: i for i, r in enumerate(order)}

# ---- A·B : 분자 x 유닛 (긴 꼴) ----
A = L[["read_id", "ap", "ad", "b_call"]].copy()
A["분자순번"] = A.read_id.map(pos)          # 평균 메틸 오름차순 (그림의 세로 순서)
A["하플로타입"] = "%s|%d" % (HAPS, HAPSLOT)
A["array_길이_유닛"] = NN
A = A.sort_values(["분자순번", "ap"])[
    ["하플로타입", "array_길이_유닛", "read_id", "분자순번", "ap", "ad", "b_call"]]

# ---- C : 자리별 분자 사이 흩어짐 ----
g = L.groupby("ap").b_call.agg(["mean", "std", "count"]).reset_index()
g.columns = ["ap", "평균_b_call", "분자사이_sd", "분자수"]
g["ad"] = NN - g.ap + 1
C = g[g.분자수 >= 5][["ap", "ad", "분자수", "평균_b_call", "분자사이_sd"]]

for df, n in ((A, "SUPP3_2_MOLECULES.tsv"), (C, "SUPP3_2_POSSD.tsv")):
    p = os.path.join(HERE, n)
    _to_english(n, df).to_csv(p, sep="\t", index=False, float_format="%.5f", lineterminator="\n")
    print("저장 %-28s 행 %6d · md5 %s"
          % (n, len(df), hashlib.md5(io.open(p, "rb").read()).hexdigest()[:12]))

print()
print("검산 — 그림 출력과 맞나")
print("   분자(read) %d        (그림 : 52)" % A.read_id.nunique())
print("   덮인 칸 %d           (그림 : 1200)" % len(A))
print("   array 유닛 %d        (그림 : 83)" % NN)
print("   read>=5 인 자리 %d" % len(C))
for nm, m in (("근위 1-10", C.ap <= 10), ("가운데 21-55", (C.ap > 20) & (C.ap <= 55)),
              ("끝 71-83", C.ap > 70)):
    s = C[m]
    print("   %-12s 평균 β %.3f · 분자사이 sd 평균 %.3f (최대 %.3f)"
          % (nm, s.평균_b_call.mean(), s.분자사이_sd.mean(), s.분자사이_sd.max()))
print("   (그림 legend : 0.093 / 0.059 / 0.026)")
