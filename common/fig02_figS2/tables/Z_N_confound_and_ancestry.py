# -*- coding: utf-8 -*-
r"""
Z_N, (C) is assembly source and quality entangled with population, (D) ancestral/derived state of terminal-type placement sets, (E) placement by population

C1. canon_src (pub public assemblies / hifi in-house / span / cons2) x population x terminal-type carriage
C2. Stratified test: is AFR lower even within the same canon_src (Mantel-Haenszel style + permutation)
C3. Are quality metrics (cov3%, var_ident%, n_unit) entangled with carriage
D1. Are the depth sets **nested**: 2 subset of 2,5 subset of 2,5,7 ...
D2. Population distribution by set size (number of terminal types), are smaller sets older
D3. Phylogenetic proximity between haplotypes carrying a set (similarity measured with ordinary unit variants)
E1. Full table of population x placement pattern
Output : $OUTDIR\ZN_*.tsv
"""
import io, os, sys, collections
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of the 12 diagnostic variants)

OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
POPF = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\v11_readvar\05_reference\1kGP_3202_ped_population_2026-09-23.txt")
STAGEF = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\5S_position_variant_work\2026-09-17_v11_fix2\STAGE_ALL_V11F2.tsv")
rng = np.random.default_rng(20261003)
NPERM = 2000
MARK = "D1096:1"
DIAG = {MARK, "635CT", "1164GA", "740GA", "DB43:2", "DB16:17", "488GC", "DB46:10", "DB58:1", "R1752CG", "R405CG", "R347CG"}
U, A_, M_ = G.load("primary", verbose=False)
HK = sorted(U)

POP = {}
for i, Ln in enumerate(io.open(POPF, encoding="utf8")):
    f = Ln.split()
    if i:
        POP[f[1]] = (f[5], f[6])

ST = {}
H = None
for Ln in io.open(STAGEF, encoding="utf8"):
    f = Ln.rstrip("\n").split("\t")
    if H is None:
        H = {c: i for i, c in enumerate(f)}; continue
    ST[(f[0], f[1])] = f

rec = []
for k in HK:
    d = [r for r in U[k] if r[1] <= 15 and r[1] <= r[0]]
    ds = tuple(sorted(r[1] for r in d if TD.is_term(r[4])))
    s = ST.get(k)
    if s is None:
        continue
    rec.append(dict(key=k, sample=k[0], pop=POP.get(k[0], ("?", "?"))[0], sup=POP.get(k[0], ("?", "?"))[1],
                    src=s[H["canon_src"]], cov3=float(s[H["cov3%"]]), vid=float(s[H["var_ident%"]]),
                    n=int(s[H["n_unit"]]), carrier=int(bool(ds)), dset=ds))
print("하플로타입 %d (STAGE 대응) · 집단표 있는 것 %d" % (len(rec), sum(1 for r in rec if r["sup"] != "?")))

# ---------- C1 ----------
print("\n=== C1. 조립 출처 × 집단 × 말단형 보유 ===")
g = io.open(os.path.join(OUT, "ZN_source_population.tsv"), "w", encoding="utf8", newline="\n")
g.write("source\tsuperpop\tn_hap\tn_carrier\tfrac\n")
srcs = sorted(set(r["src"] for r in rec))
sups = [s for s in sorted(set(r["sup"] for r in rec)) if s != "?"]
print("  출처별 전체 :")
for s in srcs:
    sub = [r for r in rec if r["src"] == s]
    print("    %-7s 하플로 %3d · 보유 %3d (%.3f)" % (s, len(sub), sum(r["carrier"] for r in sub), np.mean([r["carrier"] for r in sub])))
print("\n  출처 × 집단 (하플로타입 수 / 보유율) :")
print("    출처    " + "".join("%16s" % x for x in sups))
for s in srcs:
    line = "    %-7s" % s
    for sp in sups:
        sub = [r for r in rec if r["src"] == s and r["sup"] == sp]
        line += "%16s" % ("%d / %.3f" % (len(sub), np.mean([r["carrier"] for r in sub])) if sub else "-")
        if sub:
            g.write("%s\t%s\t%d\t%d\t%.4f\n" % (s, sp, len(sub), sum(r["carrier"] for r in sub), np.mean([r["carrier"] for r in sub])))
    print(line)
# Is the source entangled with population
tab = np.array([[sum(1 for r in rec if r["src"] == s and r["sup"] == sp) for sp in sups] for s in srcs], float)
from scipy.stats import chi2_contingency
c2, pc, dof, _ = chi2_contingency(tab[tab.sum(1) > 0])
print("\n  출처 × 집단 독립성 : chi2 %.1f · dof %d · p = %.3g  (%s)" % (
    c2, dof, pc, "얽혀 있음" if pc < 0.05 else "독립"))
g.write("#src_x_pop_chi2\t%.4f\tp\t%.6g\tdof\t%d\n" % (c2, pc, dof))

# ---------- C2 stratification ----------
print("\n=== C2. 같은 출처 안에서 AFR vs 비AFR ===")
print("  출처    AFR 하플로 보유율   비AFR 하플로 보유율   Fisher p   오즈비")
tot_or = []
for s in srcs:
    a = [r for r in rec if r["src"] == s and r["sup"] == "AFR"]
    b = [r for r in rec if r["src"] == s and r["sup"] not in ("AFR", "?")]
    if len(a) < 5 or len(b) < 5:
        print("    %-7s 표본 부족 (AFR %d · 비AFR %d)" % (s, len(a), len(b))); continue
    t = [[sum(r["carrier"] for r in a), len(a) - sum(r["carrier"] for r in a)],
         [sum(r["carrier"] for r in b), len(b) - sum(r["carrier"] for r in b)]]
    orr, pf = fisher_exact(t)
    tot_or.append((s, orr, pf, len(a), len(b)))
    print("    %-7s %3d  %.3f          %3d  %.3f          %.3g    %.3f" % (
        s, len(a), np.mean([r["carrier"] for r in a]), len(b), np.mean([r["carrier"] for r in b]), pf, orr))
    g.write("#strat\t%s\tAFR_n\t%d\tAFR_frac\t%.4f\tnonAFR_n\t%d\tnonAFR_frac\t%.4f\tfisher_p\t%.6g\tOR\t%.4f\n" % (
        s, len(a), np.mean([r["carrier"] for r in a]), len(b), np.mean([r["carrier"] for r in b]), pf, orr))

# ---------- C3 quality ----------
print("\n=== C3. 품질 지표가 보유와 얽혀 있나 (보유 vs 비보유) ===")
for nm, key in (("cov3%", "cov3"), ("var_ident%", "vid"), ("n_unit", "n")):
    a = [r[key] for r in rec if r["carrier"]]
    b = [r[key] for r in rec if not r["carrier"]]
    u, p = mannwhitneyu(a, b)
    print("    %-12s 보유 중앙 %8.2f · 비보유 %8.2f · Mann-Whitney p = %.3g" % (nm, np.median(a), np.median(b), p))
    g.write("#quality\t%s\tcarrier_median\t%.4f\tnoncarrier_median\t%.4f\tp\t%.6g\n" % (nm, np.median(a), np.median(b), p))
print("\n    같은 항목을 AFR vs 비AFR 로 :")
for nm, key in (("cov3%", "cov3"), ("var_ident%", "vid"), ("n_unit", "n")):
    a = [r[key] for r in rec if r["sup"] == "AFR"]
    b = [r[key] for r in rec if r["sup"] not in ("AFR", "?")]
    u, p = mannwhitneyu(a, b)
    print("    %-12s AFR 중앙 %8.2f · 비AFR %8.2f · p = %.3g" % (nm, np.median(a), np.median(b), p))
g.close()

# ---------- D1 nesting ----------
print("\n=== D1. 깊이 집합이 포개져 있는가 (nested) ===")
sets = collections.Counter(r["dset"] for r in rec if r["dset"])
uniq = sorted(sets, key=lambda x: (len(x), x))
n_nested = 0; n_tot = 0
parent = {}
for a in uniq:
    if len(a) < 2:
        continue
    n_tot += 1
    cands = [b for b in uniq if len(b) == len(a) - 1 and set(b) < set(a)]
    if cands:
        n_nested += 1
        parent[a] = cands[0]
print("  집합 %d 가지 · 크기 2 이상 %d 개 중 **한 개 뺀 집합이 실제로 존재** %d (%.3f)" % (
    len(uniq), n_tot, n_nested, n_nested / max(n_tot, 1)))
# Null: random sets from 1..15 with the same size distribution
kk = [len(a) for a in uniq if len(a) >= 2]
nul = []
for _ in range(500):
    us = set()
    for k2 in kk:
        us.add(tuple(sorted(rng.choice(np.arange(1, 16), k2, replace=False))))
    allu = us | set(tuple(sorted(rng.choice(np.arange(1, 16), 1))) for _ in range(len([a for a in uniq if len(a) == 1])))
    c = 0
    for a in us:
        if any(len(b) == len(a) - 1 and set(b) < set(a) for b in allu):
            c += 1
    nul.append(c / max(len(us), 1))
print("  귀무(1-15 에서 무작위 집합) %.3f [%.3f-%.3f] → %s" % (
    np.mean(nul), np.percentile(nul, 2.5), np.percentile(nul, 97.5),
    "관측이 큼 = 포개짐" if n_nested / max(n_tot, 1) > np.percentile(nul, 97.5) else "차이 없음"))

g = io.open(os.path.join(OUT, "ZN_depthset_census.tsv"), "w", encoding="utf8", newline="\n")
g.write("depth_set\tsize\tn_hap\tn_indiv\t" + "\t".join(sups) + "\tsubset_parent\n")
print("\n=== D2. 집합 크기별 집단 분포 ===")
bysize = collections.defaultdict(lambda: collections.Counter())
for r in rec:
    if r["dset"] and r["sup"] != "?":
        bysize[len(r["dset"])][r["sup"]] += 1
print("  말단형 개수   " + "".join("%8s" % s for s in sups) + "     합")
for k2 in sorted(bysize):
    c = bysize[k2]
    print("  %6d        " % k2 + "".join("%8d" % c[s] for s in sups) + "%8d" % sum(c.values()))
for a in uniq:
    c = collections.Counter(r["sup"] for r in rec if r["dset"] == a)
    ind = len(set(r["sample"] for r in rec if r["dset"] == a))
    g.write("%s\t%d\t%d\t%d\t%s\t%s\n" % (",".join(map(str, a)), len(a), sets[a], ind,
                                          "\t".join(str(c[s]) for s in sups), ",".join(map(str, parent.get(a, ()))) or "-"))
g.close()

# ---------- D3 phylogenetic proximity ----------
print("\n=== D3. 같은 집합을 가진 하플로타입끼리 계통이 가까운가 (보통 유닛 변이로) ===")
prof = {}
for r in rec:
    k = r["key"]
    s = set()
    for u in U[k]:
        if u[1] <= 15 and u[1] <= u[0] and not TD.is_term(u[4]):
            s |= (set(u[4]) - DIAG)
    prof[k] = s
jac = lambda a, b: len(a & b) / max(len(a | b), 1)
car = [r for r in rec if r["dset"]]
same = []; diff = []
for i, a in enumerate(car):
    for b in car[i + 1:]:
        if a["sample"] == b["sample"]:
            continue
        v = jac(prof[a["key"]], prof[b["key"]])
        (same if a["dset"] == b["dset"] else diff).append(v)
print("  같은 집합 쌍 %d · 평균 %.4f   |   다른 집합 쌍 %d · 평균 %.4f" % (len(same), np.mean(same), len(diff), np.mean(diff)))
u, p = mannwhitneyu(same, diff)
print("  Mann-Whitney p = %.3g" % p)
# Containment pairs (one is a subset of the other)
sub = []; unrel = []
for i, a in enumerate(car):
    for b in car[i + 1:]:
        if a["sample"] == b["sample"] or a["dset"] == b["dset"]:
            continue
        v = jac(prof[a["key"]], prof[b["key"]])
        (sub if (set(a["dset"]) < set(b["dset"]) or set(b["dset"]) < set(a["dset"])) else unrel).append(v)
print("  포함 관계 쌍 %d · %.4f   |   무관 쌍 %d · %.4f · p = %.3g" % (
    len(sub), np.mean(sub), len(unrel), np.mean(unrel), mannwhitneyu(sub, unrel)[1] if sub and unrel else float("nan")))

# ---------- E1 ----------
print("\n=== E1. 집단 × 배치 패턴 (보유 하플로타입) ===")
g = io.open(os.path.join(OUT, "ZN_pattern_by_population_full.tsv"), "w", encoding="utf8", newline="\n")
g.write("depth_set\t" + "\t".join(sups) + "\ttotal\n")
tab2 = collections.defaultdict(collections.Counter)
for r in rec:
    if r["dset"] and r["sup"] != "?":
        tab2[",".join(map(str, r["dset"]))][r["sup"]] += 1
print("  패턴           " + "".join("%7s" % s for s in sups) + "   합")
for p in sorted(tab2, key=lambda x: -sum(tab2[x].values()))[:14]:
    c = tab2[p]
    print("  %-14s" % p + "".join("%7d" % c[s] for s in sups) + "%6d" % sum(c.values()))
for p in sorted(tab2, key=lambda x: -sum(tab2[x].values())):
    c = tab2[p]
    g.write("%s\t%s\t%d\n" % (p, "\t".join(str(c[s]) for s in sups), sum(c.values())))
g.close()
print("\n저장 :", OUT)
