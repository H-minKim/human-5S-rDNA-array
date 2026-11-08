# -*- coding: utf-8 -*-
r"""
W2 - why distal terminal-type units sit 1-2 positions apart, as in '2, 5, 7' (v11-fix2, primary)

What is measured
  0. Marker check: whether the terminal-type markers D1096:1 and 635CT follow the same carrier-fraction curve by depth, and whether 1096TG is its complement.
  1. Spacing distribution vs null: per haplotype, keep the number of terminal-type units and shuffle only their positions within depths 2..15
     (2,000 sets) to build the expected spacing 1, 2, 3 and so on. If observations concentrate at 2 and 3, the arrangement is not by chance.
  2. Similarity vs spacing (lag): Jaccard of the variant sets of two units at distal depth <=15 by lag, terminal with terminal, normal with normal, and mixed.
     If a block (cassette) was repeated, similarity should recover at spacing 2-3 (periodicity).
  3. Cassette test: when terminal-type units are at depths d and d+g, how similar are the normal unit in between (d+1) and the next position (d+g+1),
     compared with random normal pairs in the same depth range.
  4. Control: the same calculation on middle units (>30 on both sides).
Output: $OUTDIR\W2_*.tsv and log
"""
import io, os, sys, collections
import numpy as np
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of 12 diagnostic variants)
OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
rng = np.random.default_rng(20260923)
NPERM = 2000
MARK = "D1096:1"
U, A_, M_ = G.load("primary", verbose=False)
HK = sorted(U)
# per haplotype: distal depth d -> variant set
DIST = {}
MIDU = {}
for k in HK:
    d = {}; mid = []
    for r in U[k]:
        ap, ad, N = r[0], r[1], r[2]
        if ad <= 20: d[ad] = r[4]
        if ap > 30 and ad > 30: mid.append((ap, r[4]))
    DIST[k] = d; MIDU[k] = sorted(mid)
# 0. marker curves
g = io.open(os.path.join(OUT, "W2_marker_by_depth.tsv"), "w", encoding="utf8", newline="\n")
g.write("depth\tn_units\tfrac_D1096_1\tfrac_635CT\tfrac_1096TG\tfrac_1164GA\tfrac_term\n")   # frac_term = share of terminal-type units (TD.is_term), frac_D1096_1 is the variant D1096:1 itself
print("=== 0. 깊이별 표지 보유율 (전 하플로타입) ===")
print("  깊이  유닛수  D1096:1  635CT  1096TG  1164GA")
prof = {}
for depth in range(1, 16):
    vs = [DIST[k][depth] for k in HK if depth in DIST[k]]
    n = len(vs)
    f = [np.mean([MARK in v for v in vs]), np.mean(["635CT" in v for v in vs]),
         np.mean(["1096TG" in v for v in vs]), np.mean(["1164GA" in v for v in vs])]
    prof[depth] = f
    g.write("%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (depth, n, f[0], f[1], f[2], f[3], np.mean([TD.is_term(v) for v in vs])))
    print("  %3d  %5d  %6.3f %6.3f %6.3f %6.3f" % (depth, n, f[0], f[1], f[2], f[3]))
g.close()
a = np.array([prof[d][0] for d in range(1, 16)]); b = np.array([prof[d][1] for d in range(1, 16)])
c = np.array([prof[d][2] for d in range(1, 16)])
print("  D1096:1 과 635CT 상관 %.4f · D1096:1 과 1096TG 상관 %.4f · (D1096:1 + 1096TG) 평균 %.3f (1 이면 완전 여집합)"
      % (np.corrcoef(a, b)[0, 1], np.corrcoef(a, c)[0, 1], float(np.mean(a + c))))
# 1. spacing distribution vs null
obs = collections.Counter(); null = np.zeros((NPERM, 16))
pos_by_hap = []
for k in HK:
    d = DIST[k]
    P = sorted([x for x in range(2, 16) if x in d and TD.is_term(d[x])])
    if len(P) >= 2: pos_by_hap.append((len(P), sorted(d.keys() & set(range(2, 16)))))
    for i in range(len(P) - 1): obs[P[i + 1] - P[i]] += 1
for j in range(NPERM):
    for npos, avail in pos_by_hap:
        pick = np.sort(rng.choice(avail, npos, replace=False))
        for i in range(npos - 1):
            gap = pick[i + 1] - pick[i]
            if gap <= 15: null[j, gap] += 1
tot_obs = sum(obs.values())
g = io.open(os.path.join(OUT, "W2_gap_vs_null.tsv"), "w", encoding="utf8", newline="\n")
g.write("gap\tn_obs\tfrac_obs\tfrac_null_mean\tnull_lo\tnull_hi\tenrichment\tp_two_sided\n")
print("\n=== 1. 말단형 유닛 사이 간격 : 관측 대 무작위 배치 귀무 (짝 %d · 하플로타입 %d) ===" % (tot_obs, len(pos_by_hap)))
print("  간격  관측   관측비   귀무비 [95%]           배수   p")
for gap in range(1, 9):
    o = obs.get(gap, 0); fo = o / max(tot_obs, 1)
    nf = null[:, gap] / np.maximum(null.sum(1), 1)
    p = (np.sum(np.abs(nf - nf.mean()) >= abs(fo - nf.mean())) + 1) / (NPERM + 1.0)
    g.write("%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.3f\t%.4g\n" % (gap, o, fo, nf.mean(), np.percentile(nf, 2.5), np.percentile(nf, 97.5), fo / max(nf.mean(), 1e-9), p))
    print("  %3d  %5d  %.3f   %.3f [%.3f–%.3f]  %5.2f  %.3g" % (gap, o, fo, nf.mean(), np.percentile(nf, 2.5), np.percentile(nf, 97.5), fo / max(nf.mean(), 1e-9), p))
g.close()
# 2. similarity vs spacing
def jac(a, b):
    if not a and not b: return 1.0
    return len(a & b) / max(len(a | b), 1)
pairs = collections.defaultdict(list)
for k in HK:
    d = DIST[k]
    for d1 in range(1, 16):
        for lag in range(1, 7):
            d2 = d1 + lag
            if d1 in d and d2 in d and d2 <= 15:
                t1, t2 = TD.is_term(d[d1]), TD.is_term(d[d2])
                kind = "말단형–말단형" if (t1 and t2) else ("보통–보통" if (not t1 and not t2) else "섞임")
                pairs[(lag, kind)].append(jac(d[d1], d[d2]))
    M = MIDU[k]
    for i in range(len(M)):
        for lag in range(1, 7):
            j = i + lag
            if j < len(M) and M[j][0] - M[i][0] == lag: pairs[(lag, "가운데 대조")].append(jac(M[i][1], M[j][1]))
g = io.open(os.path.join(OUT, "W2_similarity_by_lag.tsv"), "w", encoding="utf8", newline="\n")
g.write("lag\tkind\tn_pairs\tmean_jaccard\tmedian\n")
print("\n=== 2. 간격(lag)별 유닛 닮음 (자카드 · 원위 깊이 ≤15) ===")
print("  간격  말단형–말단형        보통–보통          섞임              가운데 대조")
for lag in range(1, 7):
    row = []
    for kind in ("말단형–말단형", "보통–보통", "섞임", "가운데 대조"):
        v = pairs.get((lag, kind), [])
        g.write("%d\t%s\t%d\t%.4f\t%.4f\n" % (lag, kind, len(v), np.mean(v) if v else float("nan"), np.median(v) if v else float("nan")))
        row.append("%.3f (n=%d)" % (np.mean(v), len(v)) if len(v) >= 5 else "—")
    print("  %3d  %-18s %-18s %-18s %s" % (lag, row[0], row[1], row[2], row[3]))
g.close()
# 3. cassette test
cas = collections.defaultdict(list)
for k in HK:
    d = DIST[k]
    P = sorted([x for x in range(2, 16) if x in d and TD.is_term(d[x])])
    for i in range(len(P) - 1):
        d1, d2 = P[i], P[i + 1]; gap = d2 - d1
        if gap not in (2, 3): continue
        if (d1 + 1) in d and (d2 + 1) in d and not TD.is_term(d[d1 + 1]) and not TD.is_term(d[d2 + 1]):
            cas[("카세트 안 보통 짝", gap)].append(jac(d[d1 + 1], d[d2 + 1]))
        ords_ = [x for x in d if 2 <= x <= 15 and not TD.is_term(d[x])]
        if len(ords_) >= 2:
            a1, a2 = rng.choice(ords_, 2, replace=False)
            cas[("무작위 보통 짝", gap)].append(jac(d[a1], d[a2]))
g = io.open(os.path.join(OUT, "W2_cassette_test.tsv"), "w", encoding="utf8", newline="\n")
g.write("kind\tgap\tn\tmean_jaccard\n")
print("\n=== 3. 카세트 검정 : 말단형 사이 보통 유닛이 되풀이된 것인가 ===")
for (kind, gap), v in sorted(cas.items()):
    g.write("%s\t%d\t%d\t%.4f\n" % (kind, gap, len(v), np.mean(v)))
    print("  %-14s 간격 %d · 짝 %4d · 자카드 %.3f" % (kind, gap, len(v), np.mean(v)))
g.close()
print("\n저장 :", OUT)
