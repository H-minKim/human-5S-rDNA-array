# -*- coding: utf-8 -*-
r"""
Z_B - ① distal group arrangement pattern of the 372 'no terminal type' haplotypes ② how different neighbouring units are in sequence (v11-fix2, primary)

① group pattern (grouped by variants of distal units only, Jaccard 0.85 single linkage, depth 1..15)
   - number of groups, share of the majority group, whether minority groups are adjacent (spacing 1) or apart (spacing >=2)
   - compared side by side with the 192 carriers
② sequence difference of neighbouring units
   - symmetric difference of the variant sets of two units converted to bp (SNV 1 bp, indel by length) -> mismatch rate over 2,240 bp
   - pair types: terminal-ordinary (neighbours), ordinary-ordinary (neighbours), terminal-terminal, by zone (proximal, middle)
   - percentiles are reported so it can be compared with recombination suppression thresholds (literature)
Output: $OUTDIR\ZB_nomark_groups.tsv, ZB_neighbour_divergence.tsv
"""
import io, os, sys, collections
import numpy as np
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of the 12 diagnostic variants)
OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
MARK = "D1096:1"; DMAX = 15; TH = 0.85; UNIT = 2240.0
U, A_, M_ = G.load("primary", verbose=False)
BLOCKS = G.BLOCKS
jac = lambda a, b: len(a & b) / max(len(a | b), 1)
def bp_of(t):
    if t.startswith("DB") or t.startswith("IB"): return int(t.split(":")[1])
    if t[0] in "DI":
        q = t[1:].split(":")
        return int(q[1]) if len(q) > 1 and q[1].isdigit() else 1
    return 1
def diff_bp(a, b): return sum(bp_of(t) for t in (a ^ b))
D = {k: {r[1]: r[4] for r in U[k] if r[1] <= DMAX} for k in sorted(U)}
P = {k: {r[0]: r[4] for r in U[k] if r[0] <= DMAX} for k in sorted(U)}
MIDU = {}
for k in sorted(U):
    m = sorted((r[0], r[4]) for r in U[k] if r[0] > 30 and r[1] > 30)
    MIDU[k] = {i + 1: v for i, (ap, v) in enumerate(m[:DMAX])}
HASMK = {k: any(TD.is_term(v) for v in D[k].values()) for k in D}
# ---------- ① group pattern ----------
def groups(dd):
    ds = sorted(dd); lab = list(range(len(ds)))
    for i in range(len(ds)):
        for j in range(i + 1, len(ds)):
            if jac(dd[ds[i]], dd[ds[j]]) >= TH:
                a, b = lab[i], lab[j]
                if a != b: lab = [a if x == b else x for x in lab]
    grp = collections.defaultdict(list)
    for i, c in enumerate(lab): grp[c].append(ds[i])
    return grp
g = io.open(os.path.join(OUT, "ZB_nomark_groups.tsv"), "w", encoding="utf8", newline="\n")
g.write("group\tsample\tslot\tn_units\tn_groups\tmajor_frac\tn_minor_groups\tminor_adjacent\tminor_spread\tpattern\n")
summ = collections.defaultdict(lambda: collections.Counter()); rows = []
for k in sorted(D):
    dd = D[k]
    if len(dd) < 12: continue
    grp = groups(dd)
    order = sorted(grp, key=lambda c: -len(grp[c]))
    ren = {c: chr(ord("A") + i) for i, c in enumerate(order)}
    pat = "".join(ren[c] for d in sorted(dd) for c, mem in grp.items() if d in mem)
    big = grp[order[0]]
    minors = [m for c, m in grp.items() if c != order[0] and len(m) >= 2]
    adj = sum(1 for m in minors if any(m[i + 1] - m[i] == 1 for i in range(len(m) - 1)))
    spr = sum(1 for m in minors if all(m[i + 1] - m[i] >= 2 for i in range(len(m) - 1)))
    gname = "말단형 보유" if HASMK[k] else "말단형 없음"
    g.write("%s\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%d\t%s\n" % (gname, k[0], k[1], len(dd), len(grp), len(big) / len(dd), len(minors), adj, spr, pat))
    rows.append((gname, len(grp), len(big) / len(dd), len(minors), adj, spr))
    s = summ[gname]
    s["하플로타입"] += 1; s["그룹수합"] += len(grp); s["다수몫합"] += len(big) / len(dd)
    s["소수그룹(2개이상 멤버)"] += len(minors); s["그중 붙은 것"] += adj; s["그중 떨어진 것"] += spr
g.close()
print("=== ① 원위 그룹 배치 : '말단형 없음' vs '보유' (원위 12 유닛 이상) ===")
print("  집단          하플로타입  평균 그룹수  다수그룹 몫  소수그룹(2+)  붙은 것  떨어진 것  떨어진 비율")
for gname in ("말단형 없음", "말단형 보유"):
    s = summ[gname]; n = s["하플로타입"]
    if not n: continue
    mn = s["소수그룹(2개이상 멤버)"]
    print("  %-12s %8d %11.2f %12.3f %12d %8d %9d %11.3f"
          % (gname, n, s["그룹수합"] / n, s["다수몫합"] / n, mn, s["그중 붙은 것"], s["그중 떨어진 것"],
             s["그중 떨어진 것"] / max(mn, 1)))
# ---------- ② sequence difference between neighbouring units ----------
acc = collections.defaultdict(list)
for k in sorted(D):
    dd = D[k]
    for d in sorted(dd):
        if d + 1 not in dd: continue
        a, b = dd[d], dd[d + 1]
        t1, t2 = TD.is_term(a), TD.is_term(b)
        kind = "말단형–말단형" if (t1 and t2) else ("말단형–보통(이웃)" if (t1 or t2) else "보통–보통(원위 이웃)")
        acc[kind].append(diff_bp(a, b) / UNIT)
    pp = P[k]
    for d in sorted(pp):
        if d + 1 in pp: acc["근위 이웃"].append(diff_bp(pp[d], pp[d + 1]) / UNIT)
    mm = MIDU[k]
    for d in sorted(mm):
        if d + 1 in mm: acc["가운데 이웃"].append(diff_bp(mm[d], mm[d + 1]) / UNIT)
g = io.open(os.path.join(OUT, "ZB_neighbour_divergence.tsv"), "w", encoding="utf8", newline="\n")
g.write("pair_kind\tn_pairs\tmedian_pct\tmean_pct\tq25\tq75\tp90\tfrac_over_1pct\tfrac_over_0.5pct\n")
print("\n=== ② 이웃한 두 유닛의 서열 불일치율 (변이 대칭차 bp / 2,240 bp) ===")
print("  짝 종류                 짝 수   중앙 %   평균 %   사분위 %        90% 지점   >1% 비율  >0.5% 비율")
for kind in ("말단형–말단형", "말단형–보통(이웃)", "보통–보통(원위 이웃)", "가운데 이웃", "근위 이웃"):
    v = np.array(acc[kind]) * 100
    if len(v) < 5: continue
    g.write("%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (kind, len(v), np.median(v), v.mean(),
            np.percentile(v, 25), np.percentile(v, 75), np.percentile(v, 90), np.mean(v > 1), np.mean(v > 0.5)))
    print("  %-22s %6d %8.3f %8.3f  %.3f–%.3f %10.3f %10.3f %10.3f"
          % (kind, len(v), np.median(v), v.mean(), np.percentile(v, 25), np.percentile(v, 75), np.percentile(v, 90),
             np.mean(v > 1), np.mean(v > 0.5)))
g.close()
print("\n저장 :", OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('ZB_neighbour_divergence.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
