# -*- coding: utf-8 -*-
# v11 core measurements
#   [0] do array_only and junction give the same value at the same position
#   [1] measurement noise floor (odd/even split of canonical CpGs)
#   [2] raw beta distribution, is it bimodal
#   [3] per-position profile
#   [4] two-state model mean beta = L + (H-L)f
#   [5] length law, boundaries Bp/Bd
#   [6] co-variation (autocovariance)
#
# Polarity (measured): 'distal' = array_pos 1 side (genomic distal, hypomethylated run 14.6)
#                        'proximal' = N-array_pos+1 side (genomic proximal, run 2.6)
#   -> ad = array_pos (distal index), ap = N-array_pos+1 (proximal index).
import io, os
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

OUT = __import__("os").environ["RDNA_DATA_ROOT"] + r"\v11_readvar"
LOG = io.open(os.path.join(OUT, r"04_logs/V3_CORE.txt"), "w", encoding="utf8", newline="\n")
rng = np.random.RandomState(20260913)


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


A = pd.read_parquet(os.path.join(OUT, r"02_tables/ru_v11.parquet"))
A = A[(A.frac_canon >= 0.9) & (A.n_canon >= 150)]
A = A[A.data_source.astype(str) == "HPRC_ONT_submissions"].copy()
A["ad"] = A.ap_dj.astype("int32")          # distal index (1 = genomic distal end)
A["ap"] = A.ap_pj.astype("int32")          # proximal index (1 = genomic proximal end)
A = A[(A.ap >= 1) & (A.ad >= 1)]
A["hi"] = (A.b > 0.5).astype(float)
A["kindS"] = A.kind.astype(str)

P("=" * 112)
P("v11 core measurements (HPRC_ONT_submissions, gate passed)")
P("=" * 112)
P()
P("  read-units %s, individuals %d, slots %d, reads %s"
  % (format(len(A), ","), A["sample"].nunique(),
     A.groupby(["sample", "slot"], observed=True).ngroups, format(A.read_id.nunique(), ",")))
P("  polarity: ad = array_pos (genomic distal index), ap = N-array_pos+1 (genomic proximal index)")
P("    names confirmed in v2 by the low-methylation run fingerprint (distal 14.6 vs proximal 2.6)")
P()

P("=" * 112)
P("[0] control: do array_only reads and junction reads give the same value at the same position")
P("=" * 112)
P()
P("  v11 has array interior reads in addition to flank (junction) reads.")
P("  the position mean β is compared where both cover the same (individual,slot,position).")
P()
g = A.groupby(["sample", "slot", "ap", "kindS"], observed=True).agg(
    n=("b", "size"), b=("b", "mean")).reset_index()
g = g[g.n >= 3]
pv = g.pivot_table(index=["sample", "slot", "ap"], columns="kindS", values="b")
pv = pv.dropna()
P("  positions with both %s" % format(len(pv), ","))
if len(pv) > 50:
    d = pv["array_only"] - pv["junction"]
    P("  difference (array_only − junction) : mean %+.4f, sd %.4f, r = %.4f"
      % (d.mean(), d.std(), float(np.corrcoef(pv["array_only"], pv["junction"])[0, 1])))
    P("  -> if the difference is near 0 and the correlation is high, the two routes can be pooled.")
P()

P("=" * 112)
P("[1] measurement noise floor: odd/even split-half of canonical CpGs")
P("=" * 112)
P()
ok = A.nA.gt(0) & A.nB.gt(0) & A.bA.notna() & A.bB.notna()
vn = ((A.loc[ok, "bA"] - A.loc[ok, "bB"]) ** 2) / 4.0
P("  usable read-units %s" % format(int(ok.sum()), ","))
P("  var_noise = E[(bA-bB)^2]/4 = %.6f  ->  noise sd = %.4f" % (vn.mean(), np.sqrt(vn.mean())))
P("  sd of observed β = %.4f, noise share of variance = %.4f"
  % (A.b.std(), vn.mean() / A.b.var()))
P("  (reference: noise sd 0.018, 2.8%% of variance)")
P()

P("=" * 112)
P("[2] raw β distribution: is it bimodal")
P("=" * 112)
P()
xs = np.linspace(0.01, 0.99, 300)
k = gaussian_kde(A.b.values[rng.choice(len(A), min(len(A), 200000), replace=False)],
                 bw_method=0.10)(xs)
pk = [xs[i] for i in range(1, len(xs) - 1)
      if k[i] > k[i - 1] and k[i] >= k[i + 1] and k[i] > 0.10 * k.max()]
P("  overall KDE peaks: %s" % "  ".join("%.2f" % x for x in pk))
h, _ = np.histogram(A.b.values, bins=25, range=(0, 1))
h = h / h.max()
P("  density " + "".join(" .:-=+*#%@"[min(int(x * 9), 9)] for x in h) + "   (0 ... 1)")
P("  β<0.35 %.3f · 0.35-0.65 %.3f · >0.65 %.3f"
  % ((A.b < 0.35).mean(), ((A.b >= 0.35) & (A.b <= 0.65)).mean(), (A.b > 0.65).mean()))
P()

P("=" * 112)
P("[3] profile by position")
P("=" * 112)
P()
P("   pos  |         proximal side (ap)          |          distal side (ad)")
P("        |  slot   read-unit   meanβ    f      |  slot   read-unit   meanβ    f")
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
P("   interior (both sides >=30) : slots %d, read-units %s, mean β %.3f, f %.3f"
  % (mid.groupby(["sample", "slot"], observed=True).ngroups, format(len(mid), ","),
     mid.b.mean(), mid.hi.mean()))
pd.DataFrame(PR, columns=["pos", "b_prox", "f_prox", "n_prox", "b_dist", "f_dist",
                          "n_dist"]).to_csv(
    os.path.join(OUT, r"02_tables/PROFILE.tsv"), sep="\t", index=False, float_format="%.4f")
P()

P("=" * 112)
P("[4] two-state model: mean β = L + (H−L)×f")
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
P("  regression over %d position zones : L = %.3f, H = %.3f, R² = %.4f" % (len(Z), bb[0], bb[0] + bb[1], r2))
P("  (reference: 81 positions, L 0.196, H 0.820, R² 0.9956)")
P()
lowz = A[(A.ad <= 8) & (A.ap >= 28)]
higz = mid
P("  state values measured separately : modal β of distal 1-8 %.3f, interior modal β %.3f"
  % (np.median(lowz.b), np.median(higz.b)))
P()

P("=" * 112)
P("[5] length law, boundaries")
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
P("  f vs N : r(f, log N) = %+.3f" % float(np.corrcoef(L.f, np.log(L.N))[0, 1]))
P("  low-methylation unit count (1-f)×N : mean %.2f, sd %.2f, r(count, N) = %+.3f"
  % (L.n_low.mean(), L.n_low.std(), float(np.corrcoef(L.n_low, L.N)[0, 1])))
P("  (reference: 21.0 ± 5.3, r −0.368, f = 1 − 23.0/N R² 0.726)")
kbest = None
for kk in np.arange(5, 60, 0.5):
    pr = 1 - kk / L.N.values
    r2_ = 1 - np.var(L.f.values - pr, ddof=1) / np.var(L.f.values, ddof=1)
    if kbest is None or r2_ > kbest[1]:
        kbest = (kk, r2_)
P("  f = 1 − k/N best-fit k = %.1f, R² = %.3f" % kbest)
P()
P("  boundaries : Bp (proximal) mean %.2f ± %.2f, Bd (distal) mean %.2f ± %.2f"
  % (L.Bp.mean(), L.Bp.std(), L.Bd.mean(), L.Bd.std()))
P("  (reference N>=81: Bp 3.47±2.14, Bd 15.39±3.18)")
L2 = L[L.N >= 81]
P("  N>=81 only (%d): Bp %.2f±%.2f, Bd %.2f±%.2f, r(Bp,Bd) = %+.3f"
  % (len(L2), L2.Bp.mean(), L2.Bp.std(), L2.Bd.mean(), L2.Bd.std(),
     float(np.corrcoef(L2.Bp, L2.Bd)[0, 1])))
P("  (reference: r(Bp,Bd) = +0.003)")
L.to_csv(os.path.join(OUT, r"02_tables/HAP_LENGTH_BOUNDARY.tsv"), sep="\t", index=False,
         float_format="%.4f")
P()

P("=" * 112)
P("[6] co-variation: autocovariance between positions within the same read")
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
P("   lag s    pairs     C(s)")
C0 = None
for s_ in sorted(acc):
    if acc[s_][1] < 500:
        continue
    c = acc[s_][0] / acc[s_][1]
    if s_ == 0:
        C0 = c
    P("   %4d   %9s   %+.5f" % (s_, format(acc[s_][1], ","), c))
P()
P("  (reference: C(s)=0.00620×exp(−s/5.00)+0.00187, half-life 3.5 units)")
P()
LOG.close()
print("done")
