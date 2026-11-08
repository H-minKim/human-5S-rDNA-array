# -*- coding: utf-8 -*-
# [archive copy] from the working copy fig6v5_step1_gmwgs.py (md5 d2067c977aff) · data and output roots via environment variables with no default (the script stops if one is unset) · module names changed to the archive names · no change to any computation
# WARNING — this script WRITES tables, into P11_TABLE_OUT_DIR (no default). Pointing it at the folder that holds the published
#   input tables (P11_TABLE_DIR) OVERWRITES them. To verify, set P11_TABLE_OUT_DIR to an empty temporary folder.
import os as _os_need, sys as _sys_need
def _need(k):
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set; it is required by %s (a module stops the script that imports it). The deposited copy has no default; see README.md section 2" % (k, _os_need.path.basename(globals().get("__file__", "an exec-loaded block"))))
    return v
# . Figure 6 version 5: the denominator of the GM12878 marker based values changes to whole genome NA12878 WGS, that is, through the same gate as the ChIP data (decided on 2026-10-14).
#   Denominator = 02_tables/denominator_comparison_2026-10-14/denom_GM_WGS_ERR3239334.tsv (gm_wgs_gate.py, 1000 Genomes ERR3239334 whole genome, same extraction and gate as K562 WGSNAT).
#   K562 is unchanged (ENCODE K562 WGS). The K562 rows are checked against the existing NEWZONE table and must agree in every cell, random number stream included, which is why K562 runs first.
#   rt_loop is the same code as in fig6v3_step1_marker_values.py (a copy). Only new files are written.
import sys, io, os, collections
import numpy as np, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
_HERE = os.path.dirname(os.path.abspath(__file__)); _ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "fig06"))
import fig06_K_N_zone_marker_lib as L
from scipy.stats import chi2
T, TAB = L.T, L.TAB
TO = _need("P11_TABLE_OUT_DIR"); os.makedirs(TO, exist_ok=True)   # writes go here (no default)
Z3 = ["근위", "가운데", "원위"]
def pci(k, e): return (chi2.ppf(0.025, 2*k)/2/e if k > 0 else 0.0), chi2.ppf(0.975, 2*k+2)/2/e
def rt_loop(specs, seed=20261008):                       # same code as fig6v3_step1
    rng = np.random.default_rng(seed); rows = []
    for cell, D, M, den, tgs in specs:
        Wd = D.get(den, {})
        for t in tgs:
            if t not in D: continue
            for z in Z3:
                vs = []
                for v, zz in M.items():
                    if zz != z: continue
                    a, dep = D[t].get(v, (0, 0)); wa, wn = Wd.get(v, (0, 0))
                    if wn == 0 or dep == 0 or wa == 0: continue
                    vs.append((a, dep*(wa/wn), dep))
                if not vs: continue
                o = sum(x[0] for x in vs); e = sum(x[1] for x in vs); d = sum(x[2] for x in vs)
                plo, phi = pci(o, e); ix = np.arange(len(vs)); bs = []
                for _ in range(2000):
                    p = rng.choice(ix, len(ix), replace=True)
                    oo = sum(vs[i][0] for i in p); ee = sum(vs[i][1] for i in p)
                    if ee > 0: bs.append(oo/ee)
                blo, bhi = (np.percentile(bs, [2.5, 97.5]) if len(bs) > 100 else (np.nan, np.nan))
                inc = "" if blo != blo else ("예" if blo <= 1 <= bhi else "아니오")
                rows.append(dict(cell=cell, target=t, zone=z, n_marker=len(vs), sum_alt=o, sum_depth=d, sum_exp=round(e, 2),
                                 ratio=round(o/e, 4) if e else "", poisson_lo=round(plo, 4), poisson_hi=round(phi, 4),
                                 boot_lo=round(blo, 4) if blo == blo else "", boot_hi=round(bhi, 4) if bhi == bhi else "",
                                 boot_includes_1=inc, denom=den))
    return pd.DataFrame(rows)
def read_gm(fn):
    D = collections.defaultdict(dict)
    for r in L.rd(fn):
        if r["filter"] == "filtered": D[r["label"]][r["variant"]] = (int(r["alt_n"]), int(r["depth"]))
    return D
def read_k562_pre(fn):
    D = collections.defaultdict(dict)
    for r in L.rd(fn): D[r["label"]][r["variant"]] = (int(r["alt_all"]), int(r["dep_all"]))
    return D
ZK = {"근위(합)": "근위", "가운데": "가운데", "원위(합)": "원위"}
M = {c: {v: r["zone"] for v, r in L.build_markers(c, "new").items()} for c in ("K562", "GM12878")}   # same marker order as step1, so that the bootstrap random number stream is identical and the result reproduces
DK = read_k562_pre("K562_QUALGATE_MIX.tsv"); DG = read_gm("READLEVEL_GM_MIX.tsv")
DG["WGS_ERR3239334"] = {r["variant"]: (int(r["alt_n"]), int(r["depth"])) for r in L.rd(os.path.join("denominator_comparison_2026-10-14", "denom_GM_WGS_ERR3239334.tsv"))}
miss = [v for v in M["GM12878"] if v not in DG["WGS_ERR3239334"] or DG["WGS_ERR3239334"][v][1] == 0]
print("GM 표지 %d · WGS 분모 없는 표지 %s" % (len(M["GM12878"]), miss or "없음"))
R = rt_loop([("K562", DK, M["K562"], "WGSNAT", ["TBP", "POLR3A", "GTF3C2", "ATAC", "SNAPC5"]),
             ("GM12878", DG, M["GM12878"], "WGS_ERR3239334", ["TBP", "POLR3G", "ATAC", "DNase", "INPUT"])])
OLD = pd.read_csv(os.path.join(T, "P11_percellline_NEWZONE_bootstrap_2026-10-13.tsv"), sep=TAB)
k_new = R[R.cell == "K562"].reset_index(drop=True); k_old = OLD[OLD.cell == "K562"].reset_index(drop=True)
def eq(a, b):
    try: return abs(float(a) - float(b)) < 1e-9
    except ValueError: return str(a) == str(b)
bad = [(c, i, a, b) for c in ("target", "zone", "n_marker", "sum_alt", "ratio", "boot_lo", "boot_hi") for i, (a, b) in enumerate(zip(k_new[c], k_old[c])) if not eq(a, b)]
print("   어긋난 칸", bad[:6])
same = len(k_new) == len(k_old) and not bad
print("검산 K562 행 = 기존 NEWZONE 표 :", "같음 (%d 행)" % len(k_new) if same else "다름")
assert same
fn = "P11_percellline_NEWZONE_GMWGS_bootstrap_2026-10-14.tsv"
R.to_csv(os.path.join(TO, fn), sep=TAB, index=False, encoding="utf-8"); print("→", fn)
for t in ("ATAC", "TBP", "DNase", "POLR3G", "INPUT"):
    a = R[(R.cell == "GM12878") & (R.target == t)].set_index("zone"); b = OLD[(OLD.cell == "GM12878") & (OLD.target == t)].set_index("zone")
    print("   GM %-6s " % t + " | ".join("%s %s→%s [%s-%s]" % (z, b.loc[z, "ratio"], a.loc[z, "ratio"], a.loc[z, "boot_lo"], a.loc[z, "boot_hi"]) for z in Z3 if z in a.index and z in b.index))
