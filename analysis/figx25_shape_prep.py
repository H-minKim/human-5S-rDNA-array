# -*- coding: utf-8 -*-
# FigX25 자료 준비 — 자료원별 ⑴ 유닛 β 가 이봉인가 ⑵ 원위→근위 자리별 패턴이 닮았나
# 같은 이름의 사본이 저장소에 하나 더 있다 — `analysis/producers_v11_2026-10-11/figx25_shape_prep.py`
#   이름이 같고 논리도 같다 (2026-10-25 바로잡음). 다른 것은 주석 과 경로를 받는 법 뿐이다
#   (`DATA_ROOT`/`TABLE_DIR` ↔ `_need` · 문자열 잇기 ↔ `os.path.join`).
#   증거 : 두 사본의 산출을 `#` 머리글을 뺀 본문 md5 로 재면 같다 —
#     `FIGX25_BIMODAL` `38ff2f78943a` · `FIGX25_DENSITY` `4012e473da5a` ·
#     `FIGX25_PROFILE` `2151a2cc4d17` (뒤의 것은 동결 기록의 값과도 같다).
#   앞판 주석(2026-10-22)은 *'내는 표가 다르다'* 고 적었다 — 그것이 틀렸다. 까닭 둘 :
#     ⒜ 파일 md5 로 비교했다 — `.tsv` 머리글에 날짜가 들어 파일 md5 는 늘 다르다.
#     ⒝ *'자료원 세 줄'* ↔ 남의 대장 *'4 rows'* 를 다른 표의 증거로 읽었다 —
#        앞은 자료 행, 뒤는 머리줄 포함 이었다.
#   그래서 이 글을 돌릴 때 주의할 것은 '다른 표' 가 아니라 덮어쓰기 다 :
#     두 사본을 같은 `TABLE_DIR` 로 돌리면 서로의 표를 덮는다 (이름이 같으므로).
#     ⇒ 검정은 반드시 제 샌드박스(`TABLE_DIR` 을 따로 세워)에서 돌린다.
#   ⇒ `Fig 3`·`S3` 가 읽는 것은 이 글이 내는 쪽 이다 (한 줄기이므로 값은 어느 쪽이든 같다).
#   
#
#   왜 이 두 가지만 보나 : 자료원 사이 절대 수준 은 플랫폼·basecaller 교란과 구분되지 않고,
#   온전히 조립된 array 는 1kGP 에 1 개뿐이라 평균 비교가 성립하지 않는다.
#   ⇒ 모양 (이봉성 · 자리별 윤곽) 만 본다. 모양은 수준 차이에 덜 휘둘린다.
#
#   방향은 이름이 아니라 실측으로 정했다 :
#     `ap = 1` 쪽 유닛 β 가 0.030 (끝에서 10 유닛까지 0.080) · 반대 끝은 0.354 · 가운데 0.815.
#     대장의 "저메틸 끝 길이 원위 14.58 · 근위 2.56" 과 맞는다 ⇒ `ap` = 원위 끝에서의 거리,
#     `ad = N - ap + 1` = 근위 끝에서의 거리. (기억 : "G8 지도 array_pos 1 = 원위")
#   자(ruler) : `ru_v11.b` = 하드콜. 유닛 값은 그 유닛을 덮은 read 들의 평균.
import io, os, sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
DATA_ROOT = os.environ["DATA_ROOT"]

# 표를 읽고 쓰는 자리 — 저장소에서는 표가 코드 옆에 없으므로 환경변수로 받는다 (2026-10-21)
HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
S = DATA_ROOT + "/"
FAMS = ["1kGP", "GIAB", "HPRC"]
EDGE = 20            # 끝에서 몇 유닛까지 그릴까


def fam(s):
    s = str(s)
    # 종양은 GIAB 곡선에 넣지 않는다 (2026-10-25).
    #   `GIAB_SOMATIC_ONT` 는 HG008-T 한 사람이고, GIAB 465 유닛의 237(51%) 이었다.
    #   `S3 F`(정상 vs 종양)에는 남는다 — 거기서는 그 역할이 맞다.
    if s.startswith("GIAB_SOMATIC"):
        return "그 밖"
    return ("HPRC" if s.startswith("HPRC") else "GIAB" if s.startswith("GIAB")
            else "1kGP" if s.startswith("1kGP") else "그 밖")


R = pd.read_parquet(S + "02_tables/ru_v11.parquet",
                    columns=["sample", "slot", "data_source", "ap", "N", "b"])
R["fam"] = R.data_source.map(fam)
R = R[R.fam.isin(FAMS)].copy()
R["hap"] = R["sample"].astype(str) + "|" + R.slot.astype(str)

# 정정 ① : HPRC 를 대장 H2(qc_eligible = 1) 로 자른다.
#   전에는 tier·qc_eligible 관문을 하나도 걸지 않아 HPRC 가 ru_v11 전수(590 슬롯)였다.
#   대장은 전부 HPRC_ONT_submissions 이라 1kGP·GIAB 에는 걸 층이 없다 ⇒ 그 둘은 그대로 둔다.
#   이름으로 거르면 안 된다 — 1kGP 18 중 14 가 H2 키와 이름이 겹친다(같은 사람이 두 자료원에).
_TI = pd.read_csv(S + "02_tables/TIERS_H0H3_SLOTS_2026-10-05.tsv", sep="	", dtype=str)
_H2 = set(_TI[_TI.tier == "H2"]["sample"] + "|" + _TI[_TI.tier == "H2"]["slot"])
_before = R[R.fam == "HPRC"].hap.nunique()
R = R[(R.fam != "HPRC") | (R.hap.isin(_H2))].copy()
print("[H2 관문] HPRC 하플로 %d -> %d (대장 H2 %d) · 1kGP·GIAB 는 층 없어 그대로"
      % (_before, R[R.fam == "HPRC"].hap.nunique(), len(_H2)))
U = R.groupby(["fam", "sample", "hap", "ap", "N"], as_index=False).b.mean()
U["ad"] = U.N - U.ap + 1
U["rel"] = (U.ap - 1) / (U.N - 1).replace(0, np.nan)      # 0 = 원위 끝, 1 = 근위 끝
U.to_csv(os.path.join(HERE, "FIGX25_UNIT.tsv"), sep="\t", index=False, float_format="%.5f")
print("유닛 %d · " % len(U) +
      " · ".join("%s 하플로 %d/개체 %d"
                 % (f, U[U.fam == f].hap.nunique(), U[U.fam == f]["sample"].nunique())
                 for f in FAMS))
print()

# ---------- ⑴ 이봉인가 ----------
G = np.linspace(0.01, 0.99, 4901)


def peaks(y):
    """국소 최대 (양 옆보다 큰 자리) · 높이 순"""
    i = np.where((y[1:-1] > y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    return sorted(i, key=lambda k: -y[k])


def shape(v, bw=0.03):
    y = gaussian_kde(v, bw_method=bw)(G)
    pk = peaks(y)
    big = [k for k in pk if y[k] >= 0.05 * y[pk[0]]]      # 가장 높은 봉의 5% 이상만 센다
    lo = G[int(np.argmax(np.where(G < 0.4, y, -1)))]
    hi = G[int(np.argmax(np.where(G > 0.6, y, -1)))]
    m = (G > lo) & (G < hi)
    j = int(np.argmin(np.where(m, y, 1e18)))
    return dict(봉우리수=len(big), 저봉=round(float(lo), 3), 고봉=round(float(hi), 3),
                골=round(float(G[j]), 4),
                골봉비=round(float(y[j] / max(y[int(np.argmax(np.where(G < 0.4, y, -1)))],
                                            y[int(np.argmax(np.where(G > 0.6, y, -1)))])), 4),
                y=y)


rows, curves = [], []
rng = np.random.default_rng(20261008)
for f in FAMS:
    d = U[U.fam == f]
    v = d.b.dropna().values
    s = shape(v)
    # split-half : 하플로타입을 반으로 갈라 골이 재현되나 (50 회)
    haps = d.hap.unique()
    by = {h: x.b.dropna().values for h, x in d.groupby("hap")}
    sh = []
    for _ in range(50):
        p = rng.permutation(haps)
        a = np.concatenate([by[h] for h in p[:len(p) // 2]])
        b2 = np.concatenate([by[h] for h in p[len(p) // 2:]])
        try:
            sh.append((shape(a)["골"], shape(b2)["골"]))
        except Exception:
            pass
    sh = np.array(sh)
    diff = np.abs(sh[:, 0] - sh[:, 1]) if len(sh) else np.array([np.nan])
    rows.append(dict(자료원=f, 하플로타입=len(haps), 개체=d["sample"].nunique(), 유닛=len(v),
                     전체평균=round(float(v.mean()), 4),
                     봉우리수=s["봉우리수"], 저봉=s["저봉"], 고봉=s["고봉"], 골=s["골"],
                     골봉비=s["골봉비"],
                     양끝몫=round(float(((v <= 0.2) | (v >= 0.8)).mean()), 4),
                     interior_share=round(float(((v > 0.2) & (v < 0.8)).mean()), 4),
                     splithalf_골차_중앙=round(float(np.median(diff)), 4),
                     splithalf_골차_90p=round(float(np.percentile(diff, 90)), 4)))
    for x, yy in zip(G[::10], s["y"][::10]):
        curves.append(dict(fam=f, beta=round(float(x), 4), density=round(float(yy), 5)))
B = pd.DataFrame(rows)
B.to_csv(os.path.join(HERE, "FIGX25_BIMODAL.tsv"), sep="\t", index=False)
pd.DataFrame(curves).to_csv(os.path.join(HERE, "FIGX25_DENSITY.tsv"), sep="\t", index=False)
pd.set_option("display.width", 240)
print("## ⑴ 유닛 β 가 이봉인가 (자료원마다)")
print(B.to_string(index=False))
print()

# ---------- ⑵ 원위 → 근위 자리별 패턴 ----------
#   원위 끝은 `ap` 로, 근위 끝은 `ad` 로 묶는다. 근위 끝은 array 끝을 알아야 하므로
#   관측 최대 자리가 N 과 같은 하플로타입만 쓴다 (그렇지 않으면 끝이 어디인지 모른다).
# 정정 ② : obs_max 를 `hap`(= sample|slot) 만으로 묶으면 자료원이 안 들어가서
#   같은 사람이 1kGP·HPRC 에 다 있을 때 다른 자료원 줄이 그 사람 array 끝을 채워 준다.
#   그 결과 근위 끝을 쓸 수 있는 1kGP 하플로가 17/18 로 부풀어 있었다 (실제 6/18).
#   적발 단서 : 1kGP 근위 n_hap 이 거리에 따라 6·6·7·9·10·11 로 늘어났다 —
#   끝까지 관측된 하플로라면 늘 수 없다 (근위 끝 유닛이 없는 것이 먼 거리에만 나타난 것).
obs = U.groupby(["fam", "hap"]).ap.max().rename("obs_max").reset_index()
U2 = U.merge(obs, on=["fam", "hap"])
U2["_ok"] = U2.obs_max == U2.N

print("   근위 끝을 쓸 수 있는 하플로타입 : " +
      " · ".join("%s %d / %d" % (f, U2[(U2.fam == f) & U2._ok].hap.nunique(),
                                 U2[U2.fam == f].hap.nunique()) for f in FAMS))

prof = []
for f in FAMS:
    d = U2[U2.fam == f]
    for k in range(1, EDGE + 1):
        v = d[d.ap == k].b
        prof.append(dict(fam=f, 쪽="원위", 거리=k, n_unit=len(v), n_hap=d[d.ap == k].hap.nunique(),
                         mean=float(v.mean()) if len(v) else np.nan,
                         median=float(v.median()) if len(v) else np.nan))
    dp = d[d._ok]
    for k in range(1, EDGE + 1):
        v = dp[dp.ad == k].b
        prof.append(dict(fam=f, 쪽="근위", 거리=k, n_unit=len(v), n_hap=dp[dp.ad == k].hap.nunique(),
                         mean=float(v.mean()) if len(v) else np.nan,
                         median=float(v.median()) if len(v) else np.nan))
    mid = d[(d.ap > EDGE) & (d.ad > EDGE)].b
    prof.append(dict(fam=f, 쪽="가운데", 거리=0, n_unit=len(mid),
                     n_hap=d[(d.ap > EDGE) & (d.ad > EDGE)].hap.nunique(),
                     mean=float(mid.mean()), median=float(mid.median())))
PR = pd.DataFrame(prof)
PR.to_csv(os.path.join(HERE, "FIGX25_PROFILE.tsv"), sep="\t", index=False, float_format="%.5f")
print()
print("## ⑵ 원위 끝에서 센 자리별 평균 β (끝 1-10)")
print(PR[(PR.쪽 == "원위") & (PR.거리 <= 10)]
      .pivot_table(index="거리", columns="fam", values="mean").round(3).to_string())
print()
print("## 근위 끝에서 센 자리별 평균 β (끝 1-10)")
print(PR[(PR.쪽 == "근위") & (PR.거리 <= 10)]
      .pivot_table(index="거리", columns="fam", values="mean").round(3).to_string())
print()
print("## 가운데")
print(PR[PR.쪽 == "가운데"][["fam", "n_unit", "n_hap", "mean"]].round(3).to_string(index=False))

# 윤곽이 닮았나 — 자료원 쌍마다 상관
print()
print("## ⑶ 윤곽이 닮았나 (끝 1-20 의 평균 β 를 자료원끼리 상관)")
for side in ("원위", "근위"):
    w = PR[(PR.쪽 == side) & (PR.거리 <= EDGE)].pivot_table(index="거리", columns="fam",
                                                          values="mean")
    out = []
    for a in range(len(FAMS)):
        for b2 in range(a + 1, len(FAMS)):
            x, y = FAMS[a], FAMS[b2]
            m = w[[x, y]].dropna()
            out.append("%s-%s rho %+.3f (n %d)" % (x, y, m[x].corr(m[y], method="spearman"), len(m)))
    print("   %-4s : %s" % (side, " · ".join(out)))
print()
print("[끝] FIGX25_UNIT · FIGX25_BIMODAL · FIGX25_DENSITY · FIGX25_PROFILE")
