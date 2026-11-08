# -*- coding: utf-8 -*-
# [archive copy] from the working copy lowmeth_changepoint.py (md5 18fc343b1add) · data and output roots via environment variables with no default (the script stops if one is unset) · module names changed to the archive names · no change to any computation
# WARNING — this script WRITES tables, into P11_TABLE_OUT_DIR (no default). Pointing it at the folder that holds the published
#   input tables (P11_TABLE_DIR) OVERWRITES them. To verify, set P11_TABLE_OUT_DIR to an empty temporary folder.
import os as _os_need, sys as _sys_need
def _need(k):
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set; it is required by %s (a module stops the script that imports it). The deposited copy has no default; see README.md section 2" % (k, _os_need.path.basename(globals().get("__file__", "an exec-loaded block"))))
    return v
# . The boundary where methylation goes from low to high, counted inward from each array end, is estimated as a change point from the gene body (2122 to 2240) CpG calls.
#   Model: units 1 to k from the end have methylation probability p1 and units k+1 to L have p2, with p1 < p2. Calls are binomial, and the (methylated calls, total calls) of each unit are summed.
#   k_hat is the maximum likelihood estimate, and the 95% interval is the profile likelihood set with delta log likelihood <= 1.92 times phi, where phi corrects for overdispersion between units, because reads tie calls together and the spread is wider than binomial.
#   The rule based method, which takes the units up to the first high one, cannot be used on hapB, where a unit carries only 1 to 3 reads. That is why this method is used. It is first matched against the rule based result on hapA and GM12878.
#   hapB: the methylation table is in the 98 unit frame and the 110 unit frame follows from it. The distal numbering is unchanged and proximal k in the 98 frame equals k+1 in the 110 frame (varlist sorted by DP, hapb_align.py: all 98 units pair up and 12 units exist only in the 110 frame).
import sys, io, os, gzip, csv, collections, math
import numpy as np
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
V = _need("V11PIPE_DIR")
T = _need("P11_TABLE_DIR"); TAB = chr(9)
G0, G1 = 2122, 2240
L = 40       # only this many units from the end are examined
def pooled(ru):
    U = collections.defaultdict(lambda: [0, 0, 0])
    with gzip.open(os.path.join(V, ru), "rt", encoding="utf-8") as g:
        for r in csv.DictReader(g, delimiter=TAB):
            if int(r["n_cpg"]) < 10: continue
            w = []
            for x in str(r["calls"]).split(","):
                if ":" not in x: continue
                p, m = x.split(":"); p = int(p)
                if G0 <= p <= G1: w.append(int(m))
            if len(w) < 3: continue
            k = (r["slot"], int(r["array_pos"]))
            U[k][0] += sum(1 for m in w if m >= 128); U[k][1] += len(w); U[k][2] += 1
    return U
def ll(m, c, p):
    p = min(max(p, 1e-6), 1 - 1e-6); return m*math.log(p) + (c - m)*math.log(1 - p)
def changepoint(units):
    """units is [(methylated, calls)] ordered as units 1 to L from the end. Returns k_hat, (lo, hi), p1, p2, phi, total calls, and the number of units with calls."""
    m = np.array([u[0] for u in units], float); c = np.array([u[1] for u in units], float)
    best = []
    for k in range(1, len(units)):
        m1, c1, m2, c2 = m[:k].sum(), c[:k].sum(), m[k:].sum(), c[k:].sum()
        if c1 == 0 or c2 == 0: best.append((-1e18, k, None, None)); continue
        p1, p2 = m1/c1, m2/c2
        if p1 >= p2: best.append((-1e18, k, p1, p2)); continue       # only the low to high direction is allowed
        best.append((ll(m1, c1, p1) + ll(m2, c2, p2), k, p1, p2))
    valid = [b for b in best if b[0] > -1e17]
    if not valid: return None
    top = max(valid, key=lambda b: b[0]); k0 = top[1]; p1, p2 = top[2], top[3]
    # overdispersion phi: the Pearson statistic per unit over the degrees of freedom (units with c >= 5)
    chi = 0.0; df = 0
    for i in range(len(units)):
        p = p1 if i < k0 else p2
        if c[i] >= 5 and 0 < p < 1:
            chi += (m[i] - c[i]*p)**2 / (c[i]*p*(1 - p)); df += 1
    phi = max(1.0, chi / max(df - 2, 1))
    inside = [b[1] for b in valid if top[0] - b[0] <= 1.92 * phi]
    return k0, (min(inside), max(inside)), p1, p2, phi, int(c.sum()), int(sum(1 for x in c if x > 0))

SETS = []   # (display name, U, slot, N which is the frame size, note)
UK = pooled("RU_K562r.tsv.gz"); UG = pooled("RU_GM12878FS.tsv.gz")
SETS = [("K562 hapA", UK, "1", 139), ("K562 hapB (98 틀)", UK, "2", 98),
        ("GM12878 hap1", UG, "1", 119), ("GM12878 hap2", UG, "2", 91)]
OUT = [TAB.join(["집합", "끝", "k_hat(해당 틀)", "구간_lo", "구간_hi", "p_저", "p_고", "φ", "전체call", "값 있는 유닛", "110틀_k_hat", "110틀_구간"])]
print("끝에서부터 저메틸 구간의 변화점 추정 (끝에서 1~%d 번 유닛 · 저→고 · 95%% 구간은 과분산 φ 로 넓힘)" % L)
print("%-18s %-3s %6s %-10s %6s %6s %5s %6s %5s" % ("집합", "끝", "k̂", "95% 구간", "p저", "p고", "φ", "call", "유닛"))
RES = {}
for nm, U, sl, N in SETS:
    for end in ("근위", "원위"):
        units = []
        for i in range(1, L+1):
            ap = (N - i + 1) if end == "근위" else i
            d = U.get((sl, ap)); units.append((d[0], d[1]) if d else (0, 0))
        r = changepoint(units)
        if r is None: print("%-18s %-3s 추정 불가" % (nm, end)); continue
        k0, (lo, hi), p1, p2, phi, tc, nu = r
        RES[(nm, end)] = r
        if "hapB" in nm:
            sh = 1 if end == "근위" else 0
            k110, i110 = k0 + sh, "%d–%d" % (lo + sh, hi + sh)
        else: k110, i110 = "", ""
        print("%-18s %-3s %6d %-10s %6.2f %6.2f %5.1f %6d %5d%s" % (nm, end, k0, "%d–%d" % (lo, hi), p1, p2, phi, tc, nu,
              ("   → 110 틀 : %s (%s)" % (k110, i110)) if k110 != "" else ""))
        OUT.append(TAB.join([nm, end, str(k0), str(lo), str(hi), "%.3f" % p1, "%.3f" % p2, "%.2f" % phi, str(tc), str(nu), str(k110), i110]))
print()
print("규칙 방식(첫 높음 앞까지) 결과 : K562 hapA 근위 14 · 원위 12 | GM12878 hap1 근위 1 · 원위 8 | hap2 근위 13 · 원위 6(MINN5)~10(MINN1)")
TO = _need("P11_TABLE_OUT_DIR"); os.makedirs(TO, exist_ok=True)   # writes go here (no default)
P = os.path.join(TO, "lowmeth_changepoint_4haplotypes_2026-10-12.tsv")
io.open(P, "w", encoding="utf-8", newline="\n").write("\n".join(OUT) + "\n")
print("→", P)
