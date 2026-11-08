# -*- coding: utf-8 -*-
r"""
Z_M - are terminal-type units (1) alike between people, and (2) linked to population (ancestry)?

M1. Terminal-type unit similarity: within the same person vs between people, control = ordinary units (same depth range).
    Also checks whether they are alike between people even AFTER removing the 12 diagnostic variants
    (before removal they are necessarily alike).
M2. Population association: join the 1kGP population table (superpopulation) and compare terminal-type carrier rates.
    - Test at the individual level (the 2 haplotypes are not independent), Fisher/chi-square, 2,000 permutations
    - Besides carrier rate, also whether the depth-set pattern (2,5,7 etc.) is linked to population
Output: $OUTDIR\ZM_*.tsv
"""
import io, os, sys, collections
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of the 12 diagnostic variants)

OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
POPF = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\v11_readvar\05_reference\1kGP_3202_ped_population_2026-09-23.txt")
rng = np.random.default_rng(20260930)
NPERM = 2000
MARK = "D1096:1"
DIAG = {MARK, "635CT", "1164GA", "740GA", "DB43:2", "DB16:17", "488GC", "DB46:10", "DB58:1", "R1752CG", "R405CG", "R347CG"}
U, A_, M_ = G.load("primary", verbose=False)
HK = sorted(U)
jac = lambda a, b: len(a & b) / max(len(a | b), 1)

# ---------- M1 ----------
term = []   # (sample, slot, depth, varset)
ordi = []
for k in HK:
    for r in U[k]:
        if not (r[1] <= 15 and r[1] <= r[0]):
            continue
        (term if TD.is_term(r[4]) else ordi).append((k[0], k[1], r[1], set(r[4])))
print("말단형 유닛 %d 개 · 같은 깊이대 보통 유닛 %d 개" % (len(term), len(ordi)))


def pairstats(pool, drop, same_person):
    """Draw 20,000 random pairs and return the mean Jaccard"""
    v = []
    n = len(pool)
    tries = 0
    while len(v) < 20000 and tries < 400000:
        tries += 1
        i, j = rng.integers(0, n, 2)
        if i == j:
            continue
        a, b = pool[i], pool[j]
        sp = (a[0] == b[0])
        if sp != same_person:
            continue
        A = a[3] - drop; B = b[3] - drop
        v.append(jac(A, B))
    return (np.mean(v), np.percentile(v, [2.5, 97.5]), len(v)) if v else (np.nan, (np.nan, np.nan), 0)


g = io.open(os.path.join(OUT, "ZM_terminal_similarity.tsv"), "w", encoding="utf8", newline="\n")
g.write("unit_class\tdiagnostic\tpair\tn_pairs\tmean_jaccard\tp2.5\tp97.5\n")
print("\n=== M1. 유닛 닮음 (원위 깊이 1-15) ===")
print("  유닛종류  진단변이   쌍         n      자카드 [2.5-97.5%]")
for cls, pool in (("말단형", term), ("보통", ordi)):
    for dname, drop in (("포함", set()), ("제외", DIAG)):
        for pname, sp in (("같은 사람", True), ("다른 사람", False)):
            m, pc, n = pairstats(pool, drop, sp)
            g.write("%s\t%s\t%s\t%d\t%.4f\t%.4f\t%.4f\n" % (cls, dname, pname, n, m, pc[0], pc[1]))
            print("  %-8s %-8s %-9s %6d   %.3f [%.3f-%.3f]" % (cls, dname, pname, n, m, pc[0], pc[1]))
g.close()

# ---------- M2 ----------
POP = {}
for i, Ln in enumerate(io.open(POPF, encoding="utf8")):
    f = Ln.split()
    if i == 0:
        continue
    POP[f[1]] = (f[5], f[6])
samples = sorted(set(k[0] for k in HK))
hit = [s for s in samples if s in POP]
print("\n=== M2. 집단 결합 ===")
print("  개체 %d 명 중 집단표에 있는 사람 %d 명 (%.1f%%)" % (len(samples), len(hit), 100.0 * len(hit) / len(samples)))
miss = [s for s in samples if s not in POP][:10]
if miss:
    print("  표에 없는 예 :", ", ".join(miss))

carrier_hap = collections.defaultdict(set)
for k in HK:
    if any(TD.is_term(r[4]) and r[1] <= 15 and r[1] <= r[0] for r in U[k]):
        carrier_hap[k[0]].add(k[1])
nhap = collections.Counter(k[0] for k in HK)

rows = []
for s in hit:
    rows.append((s, POP[s][1], POP[s][0], nhap[s], len(carrier_hap.get(s, ()))))

g = io.open(os.path.join(OUT, "ZM_population_carrier.tsv"), "w", encoding="utf8", newline="\n")
g.write("level\tgroup\tn_individuals\tn_haplotypes\tn_carrier_hap\tfrac_hap\tn_indiv_with_any\tfrac_indiv\n")
print("\n  [상위집단] 개체수  하플로타입  말단형 하플로  비율   보유 개체  개체비율")
sup = collections.defaultdict(lambda: [0, 0, 0, 0])
for s, sp, p, nh, nc in rows:
    a = sup[sp]; a[0] += 1; a[1] += nh; a[2] += nc; a[3] += (1 if nc else 0)
for spn in sorted(sup):
    a = sup[spn]
    g.write("super\t%s\t%d\t%d\t%d\t%.4f\t%d\t%.4f\n" % (spn, a[0], a[1], a[2], a[2] / max(a[1], 1), a[3], a[3] / max(a[0], 1)))
    print("  %-8s %6d %10d %13d %7.3f %10d %8.3f" % (spn, a[0], a[1], a[2], a[2] / max(a[1], 1), a[3], a[3] / max(a[0], 1)))

# Tests: haplotype-level 2 x k chi-square + individual-level permutation
tab = np.array([[sup[s][2], sup[s][1] - sup[s][2]] for s in sorted(sup)], float)
tab = tab[tab.sum(1) > 0]
chi2, pchi, dof, _ = chi2_contingency(tab)
print("\n  하플로타입 단위 카이제곱 : chi2 %.2f · dof %d · p = %.4g" % (chi2, dof, pchi))
labs = [sp for s, sp, p, nh, nc in rows]
frac = np.array([nc / max(nh, 1) for s, sp, p, nh, nc in rows])
uniq = sorted(set(labs))
obs = np.array([frac[[i for i, l in enumerate(labs) if l == u]].mean() for u in uniq])
stat = obs.max() - obs.min()
null = []
for _ in range(NPERM):
    pm = rng.permutation(frac)
    v = np.array([pm[[i for i, l in enumerate(labs) if l == u]].mean() for u in uniq])
    null.append(v.max() - v.min())
null = np.array(null)
pperm = (np.sum(null >= stat) + 1) / (NPERM + 1.0)
print("  개체 단위 순열(집단 딱지 섞기 %d 벌) : 관측 최대-최소 %.3f · 귀무 중앙 %.3f · p = %.4f" % (
    NPERM, stat, np.median(null), pperm))
g.write("#chi2_hap\t%.4f\tp\t%.6g\tdof\t%d\n" % (chi2, pchi, dof))
g.write("#perm_indiv\trange_obs\t%.4f\tp\t%.6g\tnperm\t%d\n" % (stat, pperm, NPERM))

# Subpopulations
print("\n  [하위집단 상위 12개] 개체  하플로  말단형  비율")
subp = collections.defaultdict(lambda: [0, 0, 0])
for s, sp, p, nh, nc in rows:
    a = subp[p]; a[0] += 1; a[1] += nh; a[2] += nc
for p in sorted(subp, key=lambda x: -subp[x][1])[:12]:
    a = subp[p]
    g.write("sub\t%s\t%d\t%d\t%d\t%.4f\t\t\n" % (p, a[0], a[1], a[2], a[2] / max(a[1], 1)))
    print("  %-8s %6d %8d %8d %7.3f" % (p, a[0], a[1], a[2], a[2] / max(a[1], 1)))
g.close()

# Pattern x population
print("\n=== M3. 깊이 집합 패턴 × 상위집단 (보유 하플로타입만) ===")
pat = collections.defaultdict(collections.Counter)
for k in HK:
    ds = tuple(sorted(r[1] for r in U[k] if r[1] <= 15 and r[1] <= r[0] and TD.is_term(r[4])))
    if not ds or k[0] not in POP:
        continue
    pat[POP[k[0]][1]][",".join(map(str, ds))] += 1
allp = collections.Counter()
for spn in pat:
    allp.update(pat[spn])
top = [p for p, _ in allp.most_common(6)]
print("  패턴        " + "".join("%8s" % s for s in sorted(pat)))
g = io.open(os.path.join(OUT, "ZM_pattern_by_population.tsv"), "w", encoding="utf8", newline="\n")
g.write("pattern\t" + "\t".join(sorted(pat)) + "\n")
for p in top:
    print("  %-11s" % p + "".join("%8d" % pat[s][p] for s in sorted(pat)))
    g.write("%s\t%s\n" % (p, "\t".join(str(pat[s][p]) for s in sorted(pat))))
g.close()
print("\n저장 :", OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('ZM_terminal_similarity.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
