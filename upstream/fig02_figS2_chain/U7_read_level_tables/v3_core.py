# -*- coding: utf-8 -*-
# v11 core remeasurement, corresponding to the old sections 5-1 to 5-7
#   [0] do array_only and junction give the same value at the same position (a new comparison)
#   [1] measurement noise floor (odd/even split of canonical CpGs)
#   [2] raw beta distribution, is it bimodal
#   [3] per-position profile (with the confirmed proximal/distal names)
#   [4] two-state model mean beta = L + (H-L)f
#   [5] length law, boundaries Bp/Bd
#   [6] co-variation (autocovariance)
#
# Polarity confirmed (v2 measured): old 'distal' = array_pos 1 side (genomic distal, hypomethylated run 14.6)
#                        old 'proximal' = N-array_pos+1 side (genomic proximal, run 2.6)
#   -> ad = array_pos (distal index), ap = N-array_pos+1 (proximal index). Same names as the old version.
import io, os
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

OUT = __import__("os").environ["RDNA_DATA_ROOT"] + r"\v11_readvar"
LOG = io.open(os.path.join(OUT, r"04_logs\V3_CORE.txt"), "w", encoding="utf8", newline="\n")
rng = np.random.RandomState(20260913)


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


A = pd.read_parquet(os.path.join(OUT, r"02_tables\ru_v11.parquet"))
A = A[(A.frac_canon >= 0.9) & (A.n_canon >= 150)]
A = A[A.data_source.astype(str) == "HPRC_ONT_submissions"].copy()
A["ad"] = A.ap_dj.astype("int32")          # distal index (1 = genomic distal end)
A["ap"] = A.ap_pj.astype("int32")          # proximal index (1 = genomic proximal end)
A = A[(A.ap >= 1) & (A.ad >= 1)]
A["hi"] = (A.b > 0.5).astype(float)
A["kindS"] = A.kind.astype(str)

P("=" * 112)
P("v11 핵심 재측정 (HPRC_ONT_submissions · 관문 통과)")
P("=" * 112)
P()
P("  read-유닛 %s · 개체 %d · 슬롯 %d · read %s"
  % (format(len(A), ","), A["sample"].nunique(),
     A.groupby(["sample", "slot"], observed=True).ngroups, format(A.read_id.nunique(), ",")))
P("  극성: ad = array_pos (게놈 원위 순번) · ap = N-array_pos+1 (게놈 근위 순번)")
P("    옛 판 이름과 같다 — v2 에서 저메틸 달림 지문으로 확인 (원위 14.6 vs 근위 2.6)")
P()

P("=" * 112)
P("[0] 새 대조 — array_only read 와 junction read 가 같은 자리에서 같은 값을 주나")
P("=" * 112)
P()
P("  옛 판에는 flank(junction) read 밖에 없었다. v11 은 array 안쪽 read 가 따로 있다.")
P("  같은 (개체,슬롯,자리) 를 둘 다 덮은 곳에서 자리 평균 β 를 견준다.")
P()
g = A.groupby(["sample", "slot", "ap", "kindS"], observed=True).agg(
    n=("b", "size"), b=("b", "mean")).reset_index()
g = g[g.n >= 3]
pv = g.pivot_table(index=["sample", "slot", "ap"], columns="kindS", values="b")
pv = pv.dropna()
P("  둘 다 있는 자리 %s" % format(len(pv), ","))
if len(pv) > 50:
    d = pv["array_only"] - pv["junction"]
    P("  차(array_only − junction) : 평균 %+.4f · sd %.4f · r = %.4f"
      % (d.mean(), d.std(), float(np.corrcoef(pv["array_only"], pv["junction"])[0, 1])))
    P("  ⇒ 차가 0 근처이고 상관이 높으면 두 경로를 합쳐 써도 된다.")
P()

P("=" * 112)
P("[1] 측정잡음 바닥 — 정본 CpG 홀짝 반분할")
P("=" * 112)
P()
ok = A.nA.gt(0) & A.nB.gt(0) & A.bA.notna() & A.bB.notna()
vn = ((A.loc[ok, "bA"] - A.loc[ok, "bB"]) ** 2) / 4.0
P("  쓸 수 있는 read-유닛 %s" % format(int(ok.sum()), ","))
P("  var_noise = E[(bA-bB)^2]/4 = %.6f  ->  잡음 sd = %.4f" % (vn.mean(), np.sqrt(vn.mean())))
P("  관측 β 의 sd = %.4f · 분산 중 잡음 몫 = %.4f"
  % (A.b.std(), vn.mean() / A.b.var()))
P("  (옛 판: 잡음 sd 0.018 · 분산의 2.8%%)")
P()

P("=" * 112)
P("[2] 원 β 분포 — 이봉인가")
P("=" * 112)
P()
xs = np.linspace(0.01, 0.99, 300)
k = gaussian_kde(A.b.values[rng.choice(len(A), min(len(A), 200000), replace=False)],
                 bw_method=0.10)(xs)
pk = [xs[i] for i in range(1, len(xs) - 1)
      if k[i] > k[i - 1] and k[i] >= k[i + 1] and k[i] > 0.10 * k.max()]
P("  전체 KDE 봉우리: %s" % "  ".join("%.2f" % x for x in pk))
h, _ = np.histogram(A.b.values, bins=25, range=(0, 1))
h = h / h.max()
P("  밀도 " + "".join(" .:-=+*#%@"[min(int(x * 9), 9)] for x in h) + "   (0 ... 1)")
P("  β<0.35 %.3f · 0.35-0.65 %.3f · >0.65 %.3f"
  % ((A.b < 0.35).mean(), ((A.b >= 0.35) & (A.b <= 0.65)).mean(), (A.b > 0.65).mean()))
P()

P("=" * 112)
P("[3] 자리별 옆모습")
P("=" * 112)
P()
P("   자리 |            근위 쪽 (ap)             |            원위 쪽 (ad)")
P("        |  슬롯   read-유닛   평균β    f      |  슬롯   read-유닛   평균β    f")
PR = []
for kk in list(range(1, 31)) + [35, 40, 50]:
    a = A[(A.ap == kk) & (A.ad >= 28)]
    b = A[(A.ad == kk) & (A.ap >= 28)]
    if len(a) < 100 and len(b) < 100:
        continue
    P("   %4d |  %5d  %9s   %.3f  %.3f   |  %5d  %9s   %.3f  %.3f"
      % (kk, a.groupby(["sample", "slot"], observed=True).ngroups, format(len(a), ","),
         a.b.mean(), a.hi.mean(),
         b.groupby(["sample", "slot"], observed=True).ngroups, format(len(b), ","),
         b.b.mean(), b.hi.mean()))
    PR.append((kk, a.b.mean(), a.hi.mean(), len(a), b.b.mean(), b.hi.mean(), len(b)))
mid = A[(A.ap >= 30) & (A.ad >= 30)]
P("   가운데 (양쪽 >=30) : 슬롯 %d · read-유닛 %s · 평균β %.3f · f %.3f"
  % (mid.groupby(["sample", "slot"], observed=True).ngroups, format(len(mid), ","),
     mid.b.mean(), mid.hi.mean()))
pd.DataFrame(PR, columns=["pos", "b_prox", "f_prox", "n_prox", "b_dist", "f_dist",
                          "n_dist"]).to_csv(
    os.path.join(OUT, r"02_tables\PROFILE.tsv"), sep="\t", index=False, float_format="%.4f")
P()

P("=" * 112)
P("[4] 두 상태 모형 — 평균β = L + (H−L)·f")
P("=" * 112)
P()
S = A.copy()
S["zone"] = np.where(S.ap <= 40, "P" + S.ap.astype(str),
                     np.where(S.ad <= 40, "D" + S.ad.astype(str), "MID"))
Z = S.groupby("zone", observed=True).agg(n=("b", "size"), b=("b", "mean"),
                                         f=("hi", "mean")).reset_index()
Z = Z[Z.n >= 500]
x = Z.f.values; y = Z.b.values
Xd = np.column_stack([np.ones(len(x)), x])
bb, *_ = np.linalg.lstsq(Xd, y, rcond=None)
pred = Xd @ bb
r2 = 1 - np.var(y - pred, ddof=1) / np.var(y, ddof=1)
P("  자리구역 %d 개로 회귀 : L = %.3f · H = %.3f · R² = %.4f" % (len(Z), bb[0], bb[0] + bb[1], r2))
P("  (옛 판: 자리 81 개 · L 0.196 · H 0.820 · R² 0.9956)")
P()
lowz = A[(A.ad <= 8) & (A.ap >= 28)]
higz = mid
P("  따로 잰 상태값 : 원위 1-8 의 최빈 β %.3f · 가운데 최빈 β %.3f"
  % (np.median(lowz.b), np.median(higz.b)))
P()

P("=" * 112)
P("[5] 길이 법칙 · 경계")
P("=" * 112)
P()
pp = A.groupby(["sample", "slot", "ap"], observed=True).agg(
    n=("hi", "size"), f=("hi", "mean"), N=("N", "first")).reset_index()
pp = pp[pp.n >= 3]
rows = []
for (s, sl), gg in pp.groupby(["sample", "slot"], observed=True):
    N = int(gg.N.iloc[0])
    if N < 40 or len(gg) < 20:
        continue
    fw = gg.f.mean()
    gg = gg.sort_values("ap")
    v = gg.f.values
    kp = 0
    while kp < len(v) and v[kp] < 0.5:
        kp += 1
    gg2 = gg.sort_values("ap", ascending=False)
    v2 = gg2.f.values
    kd = 0
    while kd < len(v2) and v2[kd] < 0.5:
        kd += 1
    rows.append((s, sl, N, fw, kp, kd, (1 - fw) * N))
L = pd.DataFrame(rows, columns=["sample", "slot", "N", "f", "Bp", "Bd", "n_low"])
P("  하플로타입 %d" % len(L))
P("  f 와 N : r(f, log N) = %+.3f" % float(np.corrcoef(L.f, np.log(L.N))[0, 1]))
P("  저메틸 유닛 개수 (1-f)·N : 평균 %.2f · sd %.2f · r(개수, N) = %+.3f"
  % (L.n_low.mean(), L.n_low.std(), float(np.corrcoef(L.n_low, L.N)[0, 1])))
P("  (옛 판: 21.0 ± 5.3 · r −0.368 · f = 1 − 23.0/N R² 0.726)")
kbest = None
for kk in np.arange(5, 60, 0.5):
    pr = 1 - kk / L.N.values
    r2_ = 1 - np.var(L.f.values - pr, ddof=1) / np.var(L.f.values, ddof=1)
    if kbest is None or r2_ > kbest[1]:
        kbest = (kk, r2_)
P("  f = 1 − k/N 최적 k = %.1f · R² = %.3f" % kbest)
P()
P("  경계 : Bp(근위) 평균 %.2f ± %.2f · Bd(원위) 평균 %.2f ± %.2f"
  % (L.Bp.mean(), L.Bp.std(), L.Bd.mean(), L.Bd.std()))
P("  (옛 판 N>=81: Bp 3.47±2.14 · Bd 15.39±3.18)")
L2 = L[L.N >= 81]
P("  N>=81 만 (%d): Bp %.2f±%.2f · Bd %.2f±%.2f · r(Bp,Bd) = %+.3f"
  % (len(L2), L2.Bp.mean(), L2.Bp.std(), L2.Bd.mean(), L2.Bd.std(),
     float(np.corrcoef(L2.Bp, L2.Bd)[0, 1])))
P("  (옛 판: r(Bp,Bd) = +0.003)")
L.to_csv(os.path.join(OUT, r"02_tables\HAP_LENGTH_BOUNDARY.tsv"), sep="\t", index=False,
         float_format="%.4f")
P()

P("=" * 112)
P("[6] 공동변동 — 같은 read 안 자리 간 자기공분산")
P("=" * 112)
P()
mu = A.groupby(["sample", "slot", "ap"], observed=True)["b"].transform("mean")
A["dev"] = A.b - mu
acc = {}
for rid, gg in A.groupby("read_id", observed=True):
    if len(gg) < 4:
        continue
    gg = gg.sort_values("ap")
    d = gg.dev.values.astype(float)
    ap = gg.ap.values.astype(int)
    for i in range(len(d)):
        for j in range(i, min(i + 13, len(d))):
            s_ = ap[j] - ap[i]
            if s_ > 12:
                break
            a_ = acc.setdefault(s_, [0.0, 0])
            a_[0] += d[i] * d[j]; a_[1] += 1
P("   간격 s   쌍        C(s)")
C0 = None
for s_ in sorted(acc):
    if acc[s_][1] < 500:
        continue
    c = acc[s_][0] / acc[s_][1]
    if s_ == 0:
        C0 = c
    P("   %4d   %9s   %+.5f" % (s_, format(acc[s_][1], ","), c))
P()
P("  (옛 판: C(s)=0.00620·exp(−s/5.00)+0.00187 · 반감 3.5 유닛)")
P()
LOG.close()
print("done")
