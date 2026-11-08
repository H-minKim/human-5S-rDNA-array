# -*- coding: utf-8 -*-
"""
S24: split the array into zones **as the data decide** (optimal segmentation).

Why it is needed
  The zoning used so far (terminal 1-8, middle >30, 9 zones) was a **threshold set by convention**.
  Change-point analysis showed that the boundary depth differs between proximal and distal,
  so the segmentation is applied to the per-depth profile itself to find the optimal boundaries.

Method
  1) For each depth d = 1..DMAX build a **feature vector**.
     features = 20 variant characters + methylation beta + total variants per unit.
     value at each depth = mean over the units at that depth (each haplotype has 1 unit at that depth, so it equals the haplotype mean).
     A depth is discarded if fewer than MINHAP haplotypes contribute (prevents bias from only long arrays remaining).
  2) z-standardise each feature (along the depth axis) so that features with a large scale do not dominate the split.
  3) Binary segmentation of proximal and distal **separately**.
     cost = sum of squared deviations within segments, penalty = PEN x number of features x log(DMAX), minimum segment 2 depths.
     The number of zones is chosen by BIC while varying PEN.
  4) Null: shuffle the depth order and run the same procedure, without order structure no boundary should appear.
  5) Also report what distinguishes each zone (the most prominent feature).
Output: S24_optimal_zoning.tsv, fig_v11_optimal_zoning_<preset>.png
"""
import sys, io, os, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as _fm
for _p in ((__import__("os").environ.get("RDNA_FONT_DIR", __import__("os").path.join(__import__("os").environ.get("SystemRoot", ""), "Fonts")) + "/malgun.ttf"), (__import__("os").environ.get("RDNA_FONT_DIR", __import__("os").path.join(__import__("os").environ.get("SystemRoot", ""), "Fonts")) + "/gulim.ttc")):
    if os.path.exists(_p):
        _fm.fontManager.addfont(_p)
        plt.rcParams["font.family"] = _fm.FontProperties(fname=_p).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
import unitdef_v11 as FC
import v11_io as G

PRESET = sys.argv[1] if len(sys.argv) > 1 else "primary"
OUT = os.path.join(G.OUT, PRESET); os.makedirs(OUT, exist_ok=True)
DMAX = 40; MINHAP = 100; MINSEG = 2
U, A, M = G.load(PRESET)
MT = G.load_meth(min_reads=1, verbose=False)
rng = np.random.default_rng(20260915)

FB = ["deamination_CpG", "C>T_nonCpG", "oxidation_C>A_G>T", "transition_T>C_A>G", "tv_C>G", "tv_T>A",
      "tv_T>G", "del_repeat", "del_unique", "del_ge10", "ins_repeat", "ins_unique",
      "coding_snv_xU84C", "U84C", "coding_indel", "AluY", "NTS1", "GArich", "NTS2", "n_var"]
FI = {f: i for i, f in enumerate(FB)}
NF = len(FB) + 1                      # + beta
def feat(vs):
    v = np.zeros(len(FB))
    for t in vs:
        cl, pos, ref, alt, ln, blk, rep = FC.parse(t)
        inrep = bool(rep or blk)
        v[FI["n_var"]] += 1
        if cl == "snv":
            if (ref, alt) in (("C", "T"), ("G", "A")):
                v[FI["deamination_CpG" if FC.iscpg(pos) else "C>T_nonCpG"]] += 1
            elif (ref, alt) in (("C", "A"), ("G", "T")): v[FI["oxidation_C>A_G>T"]] += 1
            elif (ref, alt) in (("T", "C"), ("A", "G")): v[FI["transition_T>C_A>G"]] += 1
            elif (ref, alt) in (("C", "G"), ("G", "C")): v[FI["tv_C>G"]] += 1
            elif (ref, alt) in (("T", "A"), ("A", "T")): v[FI["tv_T>A"]] += 1
            elif (ref, alt) in (("T", "G"), ("A", "C")): v[FI["tv_T>G"]] += 1
            if pos is not None and pos >= 2122: v[FI["U84C" if pos == 2205 else "coding_snv_xU84C"]] += 1
        else:
            v[FI[("del_repeat" if inrep else "del_unique") if cl == "del" else ("ins_repeat" if inrep else "ins_unique")]] += 1
            if cl == "del" and ln >= 10: v[FI["del_ge10"]] += 1
            if pos is not None and pos >= 2122: v[FI["coding_indel"]] += 1
        e = FC.element(pos)
        if e and e != "coding": v[FI[e]] += 1
    return v

# ---------- per-depth profile ----------
acc = {("prox", d): [] for d in range(1, DMAX + 1)}
acc.update({("dist", d): [] for d in range(1, DMAX + 1)})
accb = {k: [] for k in acc}
mid = []; midb = []
for k, rs in U.items():
    for r in rs:
        v = feat(r[4])
        m = MT.get((k[0], k[1], r[6]["array_pos"]))
        b = m["beta_call"] if (m and m["n_reads"] >= 3) else np.nan
        if r[0] <= DMAX: acc[("prox", r[0])].append(v); accb[("prox", r[0])].append(b)
        if r[1] <= DMAX: acc[("dist", r[1])].append(v); accb[("dist", r[1])].append(b)
        if r[0] > 30 and r[1] > 30: mid.append(v); midb.append(b)
MIDV = np.vstack(mid).mean(0); MIDB = np.nanmean(midb)
print("가운데 기준 : 유닛 %d · beta %.3f" % (len(mid), MIDB))

def profile(side):
    ds = []; X = []
    for d in range(1, DMAX + 1):
        a = acc[(side, d)]
        if len(a) < MINHAP: continue
        v = np.vstack(a).mean(0)
        b = np.nanmean(accb[(side, d)])
        ds.append(d); X.append(np.concatenate([v, [b]]))
    return np.array(ds), np.array(X)

def zscore(X):
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1
    return (X - mu) / sd

def best_split(S, ms=MINSEG):
    m = len(S)
    if m < 2 * ms: return None
    cs1 = np.cumsum(S, 0); cs2 = np.cumsum(S * S, 0); t1 = cs1[-1]; t2 = cs2[-1]
    mu = t1 / m; c0 = float((t2 - 2 * mu * t1 + m * mu * mu).sum())
    best = None
    for t in range(ms, m - ms + 1):
        l1 = cs1[t - 1]; l2 = cs2[t - 1]
        cl = float((l2 - l1 * l1 / t).sum())
        r1 = t1 - l1; r2 = t2 - l2; rr = m - t
        cr = float((r2 - r1 * r1 / rr).sum())
        g = c0 - (cl + cr)
        if best is None or g > best[0]: best = (g, t)
    return best

def binseg(X, pen=None, K=None):
    n = len(X); segs = [(0, n)]; cps = []
    while True:
        if K is not None and len(cps) >= K: break
        cand = None
        for si, (a, b) in enumerate(segs):
            r = best_split(X[a:b])
            if r is None: continue
            if cand is None or r[0] > cand[0]: cand = (r[0], si, r[1])
        if cand is None: break
        if pen is not None and cand[0] <= pen: break
        g, si, t = cand
        a, b = segs.pop(si); segs.append((a, a + t)); segs.append((a + t, b)); segs.sort()
        cps.append(a + t)
    return sorted(cps)

g = io.open(os.path.join(OUT, "S24_optimal_zoning.tsv"), "w", encoding="utf8", newline="\n")
g.write("side\tpen\tn_zone\tzone\tdepth_from\tdepth_to\ttop_features\n")
RES = {}
for side, nm in (("prox", "근위"), ("dist", "원위")):
    ds, X = profile(side)
    Z = zscore(X)
    print("")
    print("=" * 84)
    print("=== %s 쪽 : 깊이 %d-%d (하플로타입 %d 개 이상 기여) · 특징 %d 개 ===" % (nm, ds.min(), ds.max(), MINHAP, NF))
    for PEN in (0.5, 1.0, 2.0):
        pen = PEN * NF * np.log(len(ds))
        cps = binseg(Z, pen=pen)
        nul = [len(binseg(Z[rng.permutation(len(ds))], pen=pen)) for _ in range(200)]
        bounds = [0] + cps + [len(ds)]
        zones = [(int(ds[bounds[i]]), int(ds[bounds[i + 1] - 1])) for i in range(len(bounds) - 1)]
        print("  [PEN=%.1f] 경계 %d 개 · 귀무 평균 %.2f  ->  구역 %s" % (
            PEN, len(cps), float(np.mean(nul)), " | ".join("%d-%d" % z for z in zones)))
        for i, (a, b) in enumerate(zones):
            m = (ds >= a) & (ds <= b)
            mu = Z[m].mean(0)
            o = np.argsort(-np.abs(mu))[:4]
            names = FB + ["beta_meth"]
            ft = ", ".join("%s%+.1f" % (names[j], mu[j]) for j in o)
            g.write("%s\t%.1f\t%d\t%d\t%d\t%d\t%s\n" % (side, PEN, len(zones), i + 1, a, b, ft))
            if PEN == 1.0: print("      구역%d  깊이 %2d-%2d : %s" % (i + 1, a, b, ft))
        if PEN == 1.0: RES[side] = (ds, X, Z, cps, zones)
g.close()

# ---------- figure ----------
fig = plt.figure(figsize=(16, 9)); gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.2)
names = FB + ["beta_meth"]
KEY = ["beta_meth", "n_var", "GArich", "NTS2", "U84C", "del_repeat", "del_ge10", "coding_snv_xU84C"]
for gi, (side, nm) in enumerate((("prox", "근위 (PJ 에서)"), ("dist", "원위 (DJ 에서)"))):
    ds, X, Z, cps, zones = RES[side]
    ax = fig.add_subplot(gs[0, gi])
    for f, c in zip(KEY, ["#B03A2E", "#2E5E8E", "#4F7942", "#7D6608", "#2A9D8F", "#C06C3E", "#6A4C93", "#8D5524"]):
        j = names.index(f)
        ax.plot(ds, Z[:, j], "-", lw=1.4, color=c, label=f)
    for x in cps:
        ax.axvline(ds[x] - 0.5, color="k", lw=2, ls="--")
    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_xlabel("접합에서의 깊이 (유닛)", fontsize=9)
    ax.set_ylabel("깊이축 z 값", fontsize=9)
    ax.set_title("%s  %s — 최적 경계 %d 개\n구역 : %s" % ("AB"[gi], nm, len(cps), " | ".join("%d-%d" % z for z in zones)),
                 loc="left", fontsize=10)
    if gi == 0: ax.legend(fontsize=7, ncol=2)
    ax = fig.add_subplot(gs[1, gi])
    j = names.index("beta_meth")
    ax.plot(ds, X[:, j], "-o", color="#B03A2E", ms=3.5, lw=1.5)
    ax.axhline(MIDB, color="k", ls=":", lw=1.2)
    ax.text(ds.max() * 0.7, MIDB + 0.02, "가운데 %.2f" % MIDB, fontsize=8)
    for x in cps: ax.axvline(ds[x] - 0.5, color="k", lw=2, ls="--")
    ax.set_ylim(0, 0.95)
    ax.set_xlabel("접합에서의 깊이 (유닛)", fontsize=9); ax.set_ylabel("메틸 beta", fontsize=9)
    ax.set_title("%s  같은 경계를 메틸 위에 얹으면" % "CD"[gi], loc="left", fontsize=10)
fig.suptitle("5S rDNA — 자료가 정한 최적 구역 (깊이별 프로파일에 변화점 탐지 · 근위/원위 따로)\n"
             "v11 · %s · 특징 = 변이 성격 20 + 메틸 + 총변이 · 벌점 PEN=1.0" % PRESET, fontsize=12, y=0.98)
p = os.path.join(G.OUT, "fig_v11_optimal_zoning_%s.png" % PRESET)
fig.savefig(p, dpi=150, bbox_inches="tight"); print("\nsaved:", p)
