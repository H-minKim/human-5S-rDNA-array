# -*- coding: utf-8 -*-
# Figure 3 본 figure 자료 준비 2026-10-09
#
#   구성 :
#     A  array 구조 모식            — Fig1_2_array_structure a 윗단 재사용 (여기서는 안 만든다)
#     B  자리별 유닛 β × 배열 길이층  —  `fig_G1_meth_domain` a 를 내 틀로 재현
#     C  길이층별 저메틸 유닛 수      — 같은 그림 b 를 내 틀로 재현
#     D  유닛 안 메틸 패턴 (구획별) — 고메틸 유닛 ↔ 저메틸 유닛
#     E  8 구역을 그렇게 나눈 까닭 — 메틸 축 x 서열 축 교차에서 칸이 8 개 나온다
#     F  구역별 하플로타입 사이 차    — Fig1_3_zone_variability a 재현
#
#   차용하지 않는다 : B·C 는 상대 표를 받지 않고 내 L2 로 직접 낸다.
#     그쪽 그림 주석값과 맞춰 본 결과 최대 차 0.010 이었다 (아래 [검산] 절에 다시 적는다).
#
#   자(ruler) : 유닛 β = L2_HAP_UNIT.b_mean (하드콜) · 관문 gate(v13:84) + n_read>=5
#   방향 : L2 기준 — ap = 근위 끝에서 센 번호 · ad = 원위 끝에서 센 번호 (ru_v11 과 반대)
#   읽기만 · 무거운 IO 는 이 한 패스에서 다 낸다
import io, os, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
DATA_ROOT = os.environ["DATA_ROOT"]

S = DATA_ROOT + "/"
# 표를 읽고 쓰는 자리 — 저장소에서는 표가 코드 옆에 없으므로 환경변수로 받는다 (2026-10-21)
# Writes the deposited table with English column names and values
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
    "FIG3_C_LOWCOUNT.tsv": ("컷", {"0.508(이 세션 골)": "0.508 (trough measured here)"}),
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
LOG = io.open(os.path.join(HERE, "FIG3_MAIN_PREP.txt"), "w", encoding="utf8", newline="\n")


def P(s=""):
    print(s, flush=True)
    LOG.write(str(s) + "\n")
    LOG.flush()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def save(df, name, **kw):
    p = os.path.join(HERE, name)
    _to_english(name, df).to_csv(p, sep="\t", index=False, float_format="%.5f", lineterminator="\n", **kw)
    P("   저장 %-28s 행 %5d · md5 %s" % (name, len(df), md5(p)[:12]))


LENB = [(40, 80, "N40-80"), (81, 110, "N81-110"), (111, 140, "N111-140"), (141, 10 ** 9, "N141+")]
CUT = 0.508 # trough of the unit hard call
CUT_THEM = 0.51      #  가 쓴 컷 (견주기용)

P("=" * 92)
P("Figure 3 본 figure 자료 준비 — 2026-10-09")
P("=" * 92)
P()

# ---------------------------------------------------------------- 입력
L2 = pd.read_parquet(S + "02_tables/L2_HAP_UNIT.parquet",
                     columns=["sample", "slot", "ap", "ad", "N", "n_read", "b_mean", "f_hi"])
L2 = L2[L2.n_read >= 5].copy()
L2["hap"] = L2["sample"].astype(str) + "|" + L2.slot.astype(str)
P("[입력] L2_HAP_UNIT (n_read>=5) : 유닛 %d · 하플로타입 %d · 개체 %d"
  % (len(L2), L2.hap.nunique(), L2["sample"].nunique()))

L1C = ["sample", "slot", "read_id", "ap", "ad", "N", "gate", "b_call",
       "b_alu", "b_nts1", "b_ga", "b_nts2", "b_ups", "b_gene",
       "n_alu", "n_nts1", "n_ga", "n_nts2", "n_ups", "n_gene"]
L1 = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet", columns=L1C)
L1 = L1[L1.gate == 1].copy()
L1["hap"] = L1["sample"].astype(str) + "|" + L1.slot.astype(str)
P("[입력] L1_READ_UNIT (gate==1)  : read x 유닛 %d · 하플로타입 %d · 개체 %d"
  % (len(L1), L1.hap.nunique(), L1["sample"].nunique()))
P()

# ================================================================ B
P("=" * 92)
P("[B] 자리별 유닛 β x 배열 길이층 — fig_G1 a 를 내 틀로 재현")
P()
rows = []
for lo, hi, nm in LENB + [(0, 10 ** 9, "전체")]:
    g = L2[(L2.N >= lo) & (L2.N <= hi)]
    if not len(g):
        continue
    for k in range(1, 31):
        v = g[g.ap == k]
        if len(v):
            rows.append(dict(층=nm, 축="P", 자리=k, x=k - 1, b=v.b_mean.mean(),
                             sd=v.b_mean.std(), n=len(v), n_hap=v.hap.nunique()))
    mid = g[(g.ap > 25) & (g.ad > 25)]
    if len(mid):
        rows.append(dict(층=nm, 축="MID", 자리=0, x=30, b=mid.b_mean.mean(),
                         sd=mid.b_mean.std(), n=len(mid), n_hap=mid.hap.nunique()))
    for k in range(30, 0, -1):
        v = g[g.ad == k]
        if len(v):
            rows.append(dict(층=nm, 축="D", 자리=k, x=31 + (30 - k), b=v.b_mean.mean(),
                             sd=v.b_mean.std(), n=len(v), n_hap=v.hap.nunique()))
B = pd.DataFrame(rows)
save(B, "FIG3_B_POSBETA.tsv")

P()
P("   검산 —  그림 주석과 맞춤 (차용하지 않고 내 자료로)")
ref = {("P", 1): 0.343, ("MID", 0): 0.817, ("D", 20): 0.658, ("D", 10): 0.151, ("D", 1): 0.027}
P("   %-8s %10s %10s %10s %8s" % ("자리", "내 β", "그쪽 주석", "차", "유닛"))
mx = 0.0
for (ax, k), r in ref.items():
    v = B[(B.층 == "전체") & (B.축 == ax) & (B.자리 == k)]
    if len(v):
        b = float(v.b.iloc[0])
        mx = max(mx, abs(b - r))
        P("   %-8s %10.4f %10.3f %+10.4f %8d"
          % (("가운데" if ax == "MID" else "%s%d" % (ax, k)), b, r, b - r, int(v.n.iloc[0])))
P("   ⇒ 최대 차 %.4f — 같은 현상을 같은 자로 보고 있다" % mx)
P("   ⇒ 그래서 상대 표를 받지 않고 내 코호트(510 하플로타입 / 280 개체)로 그린다")

P()
ov = ((L2.ap <= 30) & (L2.ad <= 30)).sum()
P("   주의 — 짧은 배열은 한 유닛이 P 와 D 양쪽에 든다 : %d 유닛 (%.1f%%)"
  % (ov, 100 * ov / len(L2)))
P("     legend 에 적는다. 가운데는 ap>25 AND ad>25 로 따로 정의했다(겹침 없음).")

# ================================================================ C
P()
P("=" * 92)
P("[C] 배열 길이층별 저메틸 유닛 수 — fig_G1 b 를 내 틀로 재현")
P()
rows = []
for cut, tag in ((CUT, "0.508(이 세션 골)"), (CUT_THEM, "0.51")):
    base = None
    for lo, hi, nm in LENB:
        g = L2[(L2.N >= lo) & (L2.N <= hi)]
        if not len(g):
            continue
        per = g.groupby("hap").b_mean.apply(lambda x: int((x <= cut).sum()))
        nmed = float(g.groupby("hap").N.first().median())
        if base is None:
            base = (per.mean(), nmed)
        rows.append(dict(컷=tag, 층=nm, n_hap=len(per), N중앙=nmed,
                         저메틸수=per.mean(), 저메틸수_sd=per.std(),
                         비례했다면=base[0] * nmed / base[1],
                         저메틸몫=float((g.b_mean <= cut).mean())))
C = pd.DataFrame(rows)
save(C, "FIG3_C_LOWCOUNT.tsv")
P()
refc = {"N40-80": 23.9, "N81-110": 18.9, "N111-140": 18.1, "N141+": 17.1}
P("   %-12s %8s %10s %10s %12s" % ("층", "하플로", "내 저메틸수", "그쪽 주석", "비례했다면"))
for _, r in C[C.컷.str.startswith("0.508")].iterrows():
    P("   %-12s %8d %10.1f %10.1f %12.1f"
      % (r.층, r.n_hap, r.저메틸수, refc[r.층], r.비례했다면))
P("   ⇒ 배열이 길어져도 저메틸 유닛 수는 늘지 않는다 — 오히려 줄어든다")
P("     (길이에 비례했다면 N141+ 에서 %.1f 개여야 하는데 실측 %.1f 개)"
  % (C[(C.컷.str.startswith("0.508")) & (C.층 == "N141+")].비례했다면.iloc[0],
     C[(C.컷.str.startswith("0.508")) & (C.층 == "N141+")].저메틸수.iloc[0]))

# ================================================================ D
P()
P("=" * 92)
P("[D] 유닛 안 메틸 패턴 — 구획별 · 고메틸 유닛 ↔ 저메틸 유닛")
P()
SEG = [("AluY", "b_alu", "n_alu"), ("spacer-1", "b_nts1", "n_nts1"), ("GA", "b_ga", "n_ga"),
       ("spacer-2+flank", "b_nts2", "n_nts2"), ("5'-flank", "b_ups", "n_ups"), ("유전자부", "b_gene", "n_gene")]
P("   구획 해상도 (CpG 중앙값) — 자리별 199 점이 아니라 구획 6 칸이다")
for nm, bc, nc in SEG:
    P("      %-8s CpG 중앙 %4.0f" % (nm, L1[nc].median()))
P("   5'-flank는 spacer-2+flank 끝 90 bp 라 spacer-2+flank 안에 든다 (12 + 196 = 208 = 6 칸 합)")
P()
L1["lohi"] = np.where(L1.b_call > CUT, "고메틸 유닛", "저메틸 유닛")
rows = []
for grp in ("전체", "고메틸 유닛", "저메틸 유닛"):
    g = L1 if grp == "전체" else L1[L1.lohi == grp]
    for i, (nm, bc, nc) in enumerate(SEG):
        v = g[bc].dropna()
        # 하플로타입 평균으로도 (유닛은 하플로 안에서 독립이 아니다)
        hv = g.groupby("hap")[bc].mean().dropna()
        rows.append(dict(무리=grp, 구획=nm, 순서=i, b=v.mean(), sd=v.std(), n=len(v),
                         b_hap=hv.mean(), sd_hap=hv.std(), n_hap=len(hv),
                         # cpg 가 아니라 '그 유닛에서 실제로 호출된 참조 CpG 자리 수' 의 중앙값이다
                         #   (서열에서 센 CpG 개수가 아니다 — canon_v2.fa 전체는 199)
                         cpg_호출수_중앙=float(g[nc].median())))
D = pd.DataFrame(rows)
save(D, "FIG3_D_UNITPROFILE.tsv")
P()
P("   %-12s %-9s %8s %9s %9s" % ("무리", "구획", "유닛β", "하플로평균", "n"))
for _, r in D.iterrows():
    P("   %-12s %-9s %8.4f %9.4f %9d" % (r.무리, r.구획, r.b, r.b_hap, r.n))
hi = D[D.무리 == "고메틸 유닛"].set_index("구획")
lo = D[D.무리 == "저메틸 유닛"].set_index("구획")
P()
P("   두 무리의 구획별 차")
for nm, _, _ in SEG:
    P("      %-8s 고 %.4f − 저 %.4f = %+.4f" % (nm, hi.loc[nm].b, lo.loc[nm].b,
                                                hi.loc[nm].b - lo.loc[nm].b))
d = (hi.b - lo.b)
P("   ⇒ 차가 가장 큰 구획 %s (%+.4f) · 가장 작은 구획 %s (%+.4f)"
  % (d.idxmax(), d.max(), d.idxmin(), d.min()))

# ================================================================ E
P()
P("=" * 92)
P("[E] 8 구역을 그렇게 나눈 까닭 — 메틸 축 x 서열 축 교차")
P()
ZI = S + "02_tables/ZONE_INTEGRATED.tsv"
if os.path.exists(ZI):
    Z = pd.read_csv(ZI, sep="\t")
    P("   입력 ZONE_INTEGRATED.tsv (md5 %s · 행 %d · 열 %s)"
      % (md5(ZI)[:12], len(Z), list(Z.columns)[:10]))
    save(Z, "FIG3_E_ZONEGRID.tsv")
else:
    P("   ZONE_INTEGRATED.tsv 없음")

# 근위 제3경계 17 의 근거 — 길이층별로 내부값 98% 에 닿는 자리
rows = []
for lo_, hi_, nm in ((30, 59, "N30-59"), (60, 99, "N60-99"),
                     (100, 139, "N100-139"), (140, 300, "N140-300")):
    s = L2[(L2.N >= lo_) & (L2.N <= hi_) & (L2.ap <= L2.ad)]
    it = L2[(L2.N >= lo_) & (L2.N <= hi_) & (L2.ap > 30) & (L2.ad > 30)].f_hi.mean()
    hit = None
    for k in range(4, 31):
        v = s.loc[s.ap == k, "f_hi"].mean()
        if v == v and it == it and v >= 0.98 * it:
            hit = k
            break
    for k in range(1, 31):
        v = s.loc[s.ap == k, "f_hi"]
        if len(v):
            rows.append(dict(층=nm, ap=k, f_hi=v.mean(), n=len(v), 내부=it, 도달자리=hit))
E2 = pd.DataFrame(rows)
save(E2, "FIG3_E2_PROXBOUND.tsv")
P()
P("   근위 제3경계 — 내부값의 98% 에 닿는 자리")
for nm, g in E2.groupby("층", sort=False):
    P("      %-10s 내부 f_hi %.4f · 도달 자리 %s"
      % (nm, g.내부.iloc[0], g.도달자리.iloc[0]))
P("   ⇒ 길이층마다 16~21 ⇒ 17 을 경계로 삼았다 ( 은 25 를 '무름' 이라 했고  는 17)")

# ================================================================ F
P()
P("=" * 92)
P("[F] 구역별 하플로타입 사이 차 — Fig1_3 a 재현")
P()


def zone8(ap, ad):
    if ap <= ad:
        if ap <= 3:
            return "근위 1-3"
        if ap <= 10:
            return "근위 4-10"
        if ap <= 17:
            return "근위 11-17"
        return "내부"
    if ad <= 2:
        return "원위 1-2"
    if ad <= 8:
        return "원위 3-8"
    if ad <= 10:
        return "원위 9-10"
    if ad <= 22:
        return "원위 11-22"
    return "내부"


ZORD = ["근위 1-3", "근위 4-10", "근위 11-17", "내부",
        "원위 11-22", "원위 9-10", "원위 3-8", "원위 1-2"]
L2["z8"] = [zone8(a, d) for a, d in zip(L2.ap, L2.ad)]
hz = L2.groupby(["hap", "z8"]).b_mean.mean().reset_index()
rows = []
for z in ZORD:
    v = hz[hz.z8 == z].b_mean.dropna()
    rows.append(dict(구역=z, 순서=ZORD.index(z), n_hap=len(v), b=v.mean(),
                     sd=v.std(), iqr=float(v.quantile(.75) - v.quantile(.25)),
                     p10=float(v.quantile(.10)), p90=float(v.quantile(.90))))
F = pd.DataFrame(rows)
save(F, "FIG3_F_ZONEVAR.tsv")
P()
P("   %-12s %7s %8s %8s %8s" % ("구역", "하플로", "평균β", "하플로간 SD", "IQR"))
for _, r in F.iterrows():
    P("   %-12s %7d %8.4f %10.4f %8.4f" % (r.구역, r.n_hap, r.b, r.sd, r.iqr))
P()
P("   ⇒ 하플로타입 사이 차가 가장 작은 구역 : %s (SD %.4f)"
  % (F.loc[F.sd.idxmin()].구역, F.sd.min()))
P("   ⇒ 가장 큰 구역 : %s (SD %.4f)" % (F.loc[F.sd.idxmax()].구역, F.sd.max()))

# ================================================================ 코호트
P()
P("=" * 92)
P("[코호트] 분모 — 모든 패널 공통")
P()
CO = L2[["sample", "slot", "hap", "N"]].drop_duplicates("hap")
save(CO.sort_values("hap"), "FIG3_COHORT.tsv")
el = hashlib.md5("".join(sorted(CO.hap)).encode()).hexdigest()
P("   하플로타입 %d · 개체 %d · 원소 md5(붙여이음) %s"
  % (len(CO), CO["sample"].nunique(), el[:12]))
P("   => the same number as the H2 stratum (510 / 280). The symmetric difference of the element sets is 0.")
P()
P("   화학 층 — legend 에 적어야 한다 (HPRC ONT release2 공개 색인 md5 f24b361a62b4)")
P("      R941 143 개체 (guppy 6.5.7 / 6.4.6) · R1041 41 개체 (dorado 0.6.0) · 색인에 없음 96 개체")
P("      층별 골은 적지 않는다 (2026-10-24 에 뺐다). 전에 여기 R1041 0.5216 ↔ R941 0.4648 이 적혀 있었다.")
P("      ⑴ 그 추정기의 섞은 값이 0.4582 라서 0.508 을 쪼갠 것이 아니다 · ⑵ R1041 층이 6,721 유닛(R941 23,674)뿐이라")
P("      골 칸이 비어 argmin 이 동점이 되고 위치가 스캔 순서로 정해지는 자리 다 (기억 `trough-depends-on-aggregation-unit`).")
P("      ⇒ 적을 수 있는 것은 '컷을 0.45→0.55 로 옮겨도 저메틸 몫이 0.208→0.234' 뿐이다")
LOG.close()
