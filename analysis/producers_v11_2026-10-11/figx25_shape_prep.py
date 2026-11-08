# -*- coding: utf-8 -*-
# FigX25 자료 준비 — 자료원별 ⑴ 유닛 β 가 이봉인가 ⑵ 원위→근위 자리별 패턴이 닮았나
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


HERE = os.path.dirname(os.path.abspath(__file__))
S = _need("RDNA_V11_ROOT")
FAMS = ["1kGP", "GIAB", "HPRC"]
EDGE = 20            # 끝에서 몇 유닛까지 그릴까


def fam(s):
    # `GIAB_SOMATIC_ONT` is excluded.
    #    까닭 — 앞머리가 "GIAB" 라는 이유만으로 종양 시료 `HG008-T` 가 GIAB 묶음에
    #    들어와 있었다. 실측 : 이 묶음 465 유닛 가운데 237(51%)이 그 한 사람이었다.
    # The other outputs exclude that sample as well (`arrayprof_noHG008T_…py`).
    #    뺀 뒤 : hg002·hg003·hg004 3 개체 / 4 슬롯 / 228 유닛.
    s = str(s)
    if s.startswith("GIAB_SOMATIC"):
        return "그 밖"
    return ("HPRC" if s.startswith("HPRC") else "GIAB" if s.startswith("GIAB")
            else "1kGP" if s.startswith("1kGP") else "그 밖")


#  원본은 `S = "$DATA_ROOT/"` 로 끝에 `/` 가 붙어 있었고
#    여기서 문자열을 이어붙였다. 경로치환이 그 `/` 를 없애서 `…02_tables/…` 가 되었다.
#    ⇒ 이어붙이기를 `os.path.join` 으로 바꾼다. 가리키는 파일은 같다.
R = pd.read_parquet(os.path.join(S, "02_tables", "ru_v11.parquet"),
                    columns=["sample", "slot", "data_source", "ap", "N", "b"])
R["fam"] = R.data_source.map(fam)
R = R[R.fam.isin(FAMS)].copy()
R["hap"] = R["sample"].astype(str) + "|" + R.slot.astype(str)

# 정정 ① : HPRC 를 대장 H2(qc_eligible = 1) 로 자른다.
#   전에는 tier·qc_eligible 관문을 하나도 걸지 않아 HPRC 가 ru_v11 전수(590 슬롯)였다.
#   대장은 전부 HPRC_ONT_submissions 이라 1kGP·GIAB 에는 걸 층이 없다 ⇒ 그 둘은 그대로 둔다.
#   이름으로 거르면 안 된다 — 1kGP 18 중 14 가 H2 키와 이름이 겹친다(같은 사람이 두 자료원에).
_TI = pd.read_csv(os.path.join(S, "02_tables", "TIERS_H0H3_SLOTS_2026-10-05.tsv"),
                  sep="	", dtype=str)
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
                     가운데몫=round(float(((v > 0.2) & (v < 0.8)).mean()), 4),
                     splithalf_골차_중앙=round(float(np.median(diff)), 4),
                     splithalf_골차_90p=round(float(np.percentile(diff, 90)), 4)))
    for x, yy in zip(G[::10], s["y"][::10]):
        curves.append(dict(fam=f, beta=round(float(x), 4), density=round(float(yy), 5)))
B = pd.DataFrame(rows)
#  FIGX25_BIMODAL.tsv is written with English column names (the source values were already
#    English). B keeps its names for the printout below. Reader : figS4_A_unit_beta_three_sources.py
B.rename(columns={"자료원": "source", "하플로타입": "haplotypes", "개체": "individuals", "유닛": "units",
                  "전체평균": "mean_beta", "봉우리수": "n_peaks", "저봉": "low_peak", "고봉": "high_peak",
                  "골": "valley", "골봉비": "valley_to_peak_ratio", "양끝몫": "fraction_at_ends",
                  "가운데몫": "fraction_in_middle", "splithalf_골차_중앙": "splithalf_valley_difference_median",
                  "splithalf_골차_90p": "splithalf_valley_difference_p90"}
         ).to_csv(os.path.join(HERE, "FIGX25_BIMODAL.tsv"), sep="\t", index=False)
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
