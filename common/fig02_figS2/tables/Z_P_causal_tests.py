# -*- coding: utf-8 -*-
r"""
Z_P, the three key tests that support the paper's claims

P1 (most important), dose-response: is 'reduced exchange' a **function of sequence difference**, or a special property of terminal types only
    To avoid circular reasoning each unit is split into two segments.
      x = **sequence difference** measured in the first half (1..1120)
      y = **variant sharing** measured in the second half (1121..2240)
    All pairs of distal 1-15 units within the same haplotype are used, and (ordinary, ordinary) pairs and (terminal, ordinary) pairs are plotted separately.
    If the two curves overlap -> the difference alone explains it (terminal types are not special) = supports the mechanism claim of the model
      If terminal types fall below the curve -> something besides the difference is at work
P2, read-level mixing test: redo the mixing test on **read consensus** (READMAP read_varlist) instead of assemblies
P3, shape of the spacing distribution: is it 'anything 2 or more', or 'a peak at 2-3'
Output : $OUTDIR\ZP_*.tsv
"""
import io, os, sys, collections
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of the 12 diagnostic variants)

OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
rng = np.random.default_rng(20261005)
MARK = "D1096:1"
DIAG = {MARK, "635CT", "1164GA", "740GA", "DB43:2", "DB16:17", "488GC", "DB46:10", "DB58:1", "R1752CG", "R405CG", "R347CG"}
UNIT, HALF = 2240, 1120
U, A_, M_ = G.load("primary", verbose=False)
HK = sorted(U)
BLOCKS = G.BLOCKS


def pz(t):
    if t.startswith(("DB", "IB")):
        b = int(t[2:].split(":")[0])
        return BLOCKS[b - 1][0] if 1 <= b <= len(BLOCKS) else None
    if t[0] in "DI":
        return int(t[1:].split(":")[0])
    if t[0] == "R":
        return int(t[1:-2])
    try:
        return int(t[:-2])
    except ValueError:
        return None


def bp(t):
    if t.startswith(("DB", "IB")):
        b = int(t[2:].split(":")[0])
        a, e = BLOCKS[b - 1]
        return e - a + 1
    if t[0] == "D":
        p = t[1:].split(":")
        return int(p[1]) if len(p) > 1 and p[1].isdigit() else 1
    if t[0] == "I":
        p = t[1:].split(":")
        return len(p[2]) if len(p) > 2 and p[2] not in ("?", "") else (int(p[1]) if len(p) > 1 and p[1].isdigit() else 1)
    return 1


# ---------------- P1 ----------------
print("=== P1. 용량-반응 : 앞쪽 절반의 서열 차이 → 뒤쪽 절반의 변이 공유 ===")
print("  (순환논증 회피 : 차이와 공유를 유닛의 서로 다른 구간에서 잰다)")
rows = []
for k in HK:
    d = {}
    for r in U[k]:
        if r[1] <= 15 and r[1] <= r[0]:
            v = set(r[4]) - DIAG
            d[r[1]] = (set(t for t in v if (pz(t) or 0) <= HALF), set(t for t in v if (pz(t) or 0) > HALF), TD.is_term(r[4]))
    ps = sorted(d)
    for i, a in enumerate(ps):
        for b in ps[i + 1:]:
            A, B = d[a], d[b]
            # x : first-half difference % (symmetric difference bp / 1120)
            x = 100.0 * sum(bp(t) for t in (A[0] ^ B[0])) / HALF
            # y : second-half conditional sharing (mean of both directions)
            if not A[1] and not B[1]:
                continue
            y1 = len(A[1] & B[1]) / len(A[1]) if A[1] else np.nan
            y2 = len(A[1] & B[1]) / len(B[1]) if B[1] else np.nan
            y = np.nanmean([y1, y2])
            kind = "말단형-말단형" if (A[2] and B[2]) else ("말단형-보통" if (A[2] or B[2]) else "보통-보통")
            rows.append((k, x, y, kind, b - a))
print("  쌍 %d (보통-보통 %d · 말단형-보통 %d · 말단형-말단형 %d)" % (
    len(rows), sum(1 for r in rows if r[3] == "보통-보통"), sum(1 for r in rows if r[3] == "말단형-보통"),
    sum(1 for r in rows if r[3] == "말단형-말단형")))
BINS = [0, 0.2, 0.4, 0.6, 0.9, 1.3, 1.8, 2.5, 3.5, 100]
g = io.open(os.path.join(OUT, "ZP_dose_response.tsv"), "w", encoding="utf8", newline="\n")
g.write("bin_lo\tbin_hi\tkind\tn_pairs\tmean_div_front\tmean_share_back\n")
print("\n  앞절반 차이%%     보통-보통 (n · 뒤절반 공유)      말단형-보통 (n · 공유)     차이")
for lo, hi in zip(BINS, BINS[1:]):
    o = [r for r in rows if r[3] == "보통-보통" and lo <= r[1] < hi]
    t = [r for r in rows if r[3] == "말단형-보통" and lo <= r[1] < hi]
    so = np.nanmean([r[2] for r in o]) if o else np.nan
    st = np.nanmean([r[2] for r in t]) if t else np.nan
    for nm, s, arr in (("보통-보통", so, o), ("말단형-보통", st, t)):
        if arr:
            g.write("%.2f\t%.2f\t%s\t%d\t%.4f\t%.4f\n" % (lo, hi, nm, len(arr), np.mean([r[1] for r in arr]), s))
    print("   %4.1f-%-5.1f    %6d  %.3f              %6d  %s        %s" % (
        lo, hi, len(o), so if o else float("nan"), len(t),
        ("%.3f" % st) if t else "  -  ", ("%+.3f" % (st - so)) if (o and t) else "   -"))
g.close()
oo = [r for r in rows if r[3] == "보통-보통"]
print("\n  보통-보통 쌍 안에서 Spearman(앞절반 차이, 뒤절반 공유) = %+.3f (p %.3g · n %d)" % (
    *spearmanr([r[1] for r in oo], [r[2] for r in oo])[:2], len(oo)))
# Difference-matched comparison: are terminal pairs lower within the same difference bin
dif = []
for lo, hi in zip(BINS, BINS[1:]):
    o = [r[2] for r in rows if r[3] == "보통-보통" and lo <= r[1] < hi]
    t = [r[2] for r in rows if r[3] == "말단형-보통" and lo <= r[1] < hi]
    if len(o) >= 30 and len(t) >= 30:
        dif.append((lo, hi, np.nanmean(t) - np.nanmean(o), mannwhitneyu(t, o, nan_policy="omit")[1], len(t), len(o)))
print("\n  차이 맞춤 칸(양쪽 30 쌍 이상) : 말단형 − 보통")
for lo, hi, dd, p, nt, no in dif:
    print("    %.1f-%.1f%%  차 %+.3f  (p %.3g · 말단형 %d · 보통 %d)" % (lo, hi, dd, p, nt, no))
if dif:
    print("  → 칸 %d 개 중 |차| < 0.05 인 칸 %d 개" % (len(dif), sum(1 for x in dif if abs(x[2]) < 0.05)))

# ---------------- P2 ----------------
print("\n=== P2. read 합의(조립본 아님)로 섞임 검정 다시 ===")
RV = collections.defaultdict(dict)
for f, H in G._rows(G.READMAPF):
    if int(f[H["depth"]]) < 10:
        continue
    s = set()
    for t in (f[H["read_varlist"]].split(",") if f[H["read_varlist"]] not in ("", "-") else []):
        try:
            v = G.to_v4(t)
            if v:
                s.add(v)
        except ValueError:
            pass
    RV[(f[H["sample"]], f[H["slot"]])][int(f[H["array_pos"]])] = s
AP = {}
for k in HK:
    for r in U[k]:
        AP[(k[0], k[1], r[6]["array_pos"])] = r[1]
obs = collections.Counter(); perhap = []
for k in HK:
    ap = RV.get(k, {})
    d = {}
    for a, s in ap.items():
        dep = AP.get((k[0], k[1], a))
        if dep is not None and dep <= 15:
            d[dep] = s
    pos = sorted(d)
    if len(pos) < 6:
        continue
    isT = [TD.is_term(d[p]) for p in pos]
    if not any(isT):
        continue
    n = len(pos)
    cnt = collections.Counter()
    for p in pos:
        cnt.update(t for t in d[p] if t not in DIAG)
    toks = []
    for t, c in cnt.items():
        if not (2 <= c < n / 2.0):
            continue
        car = [i for i, p in enumerate(pos) if t in d[p]]
        toks.append(car)
        kt = sum(isT[i] for i in car)
        obs["보통만" if kt == 0 else ("말단형만" if kt == len(car) else "섞임")] += 1
    if toks:
        perhap.append((isT, toks))
tot = sum(obs.values())
print("  하플로타입 %d · 소수 토큰 %d" % (len(perhap), tot))
if tot:
    print("  관측 : 보통만 %.3f · 말단형만 %.3f · **섞임 %.3f**" % (
        obs["보통만"] / tot, obs["말단형만"] / tot, obs["섞임"] / tot))
    nul = []
    for _ in range(1000):
        c = collections.Counter()
        for isT, toks in perhap:
            sh = list(isT); rng.shuffle(sh)
            for car in toks:
                kt = sum(sh[i] for i in car)
                c["보통만" if kt == 0 else ("말단형만" if kt == len(car) else "섞임")] += 1
        nul.append(c["섞임"] / sum(c.values()))
    nul = np.array(nul)
    p = (np.sum(nul <= obs["섞임"] / tot) + 1) / 1001.0
    print("  귀무 %.3f [%.3f-%.3f] · p = %.4f  (조립본 판 : 관측 0.052 · 귀무 0.530)" % (
        nul.mean(), np.percentile(nul, 2.5), np.percentile(nul, 97.5), p))
    g = io.open(os.path.join(OUT, "ZP_readlevel_mixing.tsv"), "w", encoding="utf8", newline="\n")
    g.write("layer\tn_hap\tn_tokens\tobs_mixed\tnull_mixed\tnull_lo\tnull_hi\tp\n")
    g.write("read_consensus\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (
        len(perhap), tot, obs["섞임"] / tot, nul.mean(), np.percentile(nul, 2.5), np.percentile(nul, 97.5), p))
    g.close()

# ---------------- P3 ----------------
print("\n=== P3. 간격 분포의 모양 — '2 이상 아무거나' vs '2-3 봉우리' ===")
obsg = collections.Counter(); n_t = []
for k in HK:
    d = {r[1]: TD.is_term(r[4]) for r in U[k] if r[1] <= 15 and r[1] <= r[0]}
    ts = sorted(p for p in d if d[p])
    if len(ts) < 2:
        continue
    n_t.append((len(d), len(ts)))
    for x, y in zip(ts, ts[1:]):
        obsg[y - x] += 1
tot = sum(obsg.values())
# Null A: fully random, null B: random allowing only spacing 2 or more
nA = collections.Counter(); nB = collections.Counter()
for _ in range(2000):
    for N, kk in n_t:
        s = sorted(rng.choice(N, kk, replace=False))
        for x, y in zip(s, s[1:]):
            nA[y - x] += 1
        for _try in range(60):
            s2 = sorted(rng.choice(N, kk, replace=False))
            if all(y - x >= 2 for x, y in zip(s2, s2[1:])):
                for x, y in zip(s2, s2[1:]):
                    nB[y - x] += 1
                break
tA = sum(nA.values()); tB = sum(nB.values())
g = io.open(os.path.join(OUT, "ZP_gap_shape.tsv"), "w", encoding="utf8", newline="\n")
g.write("gap\tobs\tobs_frac\tnull_free\tnull_min2\tenrich_vs_min2\n")
print("  간격  관측    관측비율   귀무A(자유)   귀무B(간격>=2 강제)   B 대비 배수")
for gp in range(1, 9):
    o = obsg[gp] / max(tot, 1); a = nA[gp] / max(tA, 1); b = nB[gp] / max(tB, 1)
    g.write("%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n" % (gp, obsg[gp], o, a, b, o / max(b, 1e-9)))
    print("   %d   %5d   %.3f      %.3f        %.3f            %.2f" % (gp, obsg[gp], o, a, b, o / max(b, 1e-9)))
g.close()
print("\n저장 :", OUT)
