# -*- coding: utf-8 -*-
# [archive copy] from the working copy fig6v3_step1_marker_values.py (md5 a6183bc54f37) · data and output roots via environment variables with no default (the script stops if one is unset) · module names changed to the archive names · no change to any computation
# WARNING — this script WRITES tables, into P11_TABLE_OUT_DIR (no default). Pointing it at the folder that holds the published
#   input tables (P11_TABLE_DIR) OVERWRITES them. To verify, set P11_TABLE_OUT_DIR to an empty temporary folder.
import os as _os_need, sys as _sys_need
def _need(k):
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set; it is required by %s (a module stops the script that imports it). The deposited copy has no default; see README.md section 2" % (k, _os_need.path.basename(globals().get("__file__", "an exec-loaded block"))))
    return v
# . Figure 6 v3 Step 1: the marker based values (e, f, g, h) are moved onto the new zones.
#   ① The new marker table and the new count table are saved. ② The zone table (RT) uses the same loop as fix_ci_and_figs.py, and it is first run on the old zones to reproduce the existing RT.
#   ③ The interaction against input uses the same code as inter() in p11_asb_bontable_v4_caveat.py, and it too is first run on the old zones to reproduce the existing BT.
#   Only new files are written, so the existing reference tables do not change.
import sys, io, os, re, csv, math, collections
import numpy as np, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
_HERE = os.path.dirname(os.path.abspath(__file__)); _ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "fig06"))
import fig06_K_N_zone_marker_lib as L
from scipy.stats import chi2
T, TAB = L.T, L.TAB
TO = _need("P11_TABLE_OUT_DIR"); os.makedirs(TO, exist_ok=True)   # writes go here (no default)
Z3 = ["근위", "가운데", "원위"]
def pci(k, e):
    return (chi2.ppf(0.025, 2*k)/2/e if k > 0 else 0.0), chi2.ppf(0.975, 2*k+2)/2/e
# ---------------- reading the counts
def read_gm(fn):
    D = collections.defaultdict(dict)
    for r in L.rd(fn):
        if r["filter"] == "filtered": D[r["label"]][r["variant"]] = (int(r["alt_n"]), int(r["depth"]))
    return D
def read_k562_pre(fn, kind):
    D = collections.defaultdict(dict)
    if kind == "READLEVEL":
        for r in L.rd(fn):
            if r["filter"] == "filtered": D[r["label"]][r["variant"]] = (int(r["alt_n"]), int(r["depth"]))
    else:   # QUALGATE: before the gate, which is alt_all over dep_all
        for r in L.rd(fn): D[r["label"]][r["variant"]] = (int(r["alt_all"]), int(r["dep_all"]))
    return D
# ---------------- ② the zone table, the same loop as the first block of fix_ci_and_figs.py
def rt_loop(specs, seed=20261008):
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
ZK3 = {"근위(합)": "근위", "가운데": "가운데", "원위(합)": "원위"}
def mtab(fn): return {r["variant"]: ZK3.get(r["zone"]) for r in L.rd(fn)}
# Reproducing the old zones: compared against the K562 and GM12878 rows of the existing RT, random number stream included, so K562 runs first and GM12878 second
OLDRT = pd.read_csv(os.path.join(T, "P11_percellline_latest_bootstrap_2026-10-08.tsv"), sep=TAB)
spec_old = [("K562", read_k562_pre("READLEVEL_3ZONE_V4.tsv", "READLEVEL"), mtab("P11_K562_MARKERS_3ZONE_v4_final.tsv"), "WGSNAT",
             ["TBP", "POLR3A", "BDP1", "BRF1", "BRF2", "GTF3C2", "POLR3G", "ATAC", "SNAPC5"]),
            ("GM12878", read_gm("READLEVEL_GM.tsv"), mtab("P11_GM12878_MARKERS_3ZONE.tsv"), "HIFI", ["TBP", "POLR3G", "ATAC", "DNase", "INPUT"])]
R_old = rt_loop(spec_old)
mm = 0; nn = 0
for _, r in R_old.iterrows():
    o = OLDRT[(OLDRT.cell == r.cell) & (OLDRT.target == r.target) & (OLDRT.zone == r.zone)]
    if len(o) != 1: mm += 1; continue
    o = o.iloc[0]; nn += 1
    same = (abs(float(o.ratio) - float(r.ratio)) < 1e-9) and (str(o.boot_lo) == str(r.boot_lo) or abs(float(o.boot_lo) - float(r.boot_lo)) < 1e-9) and (abs(float(o.boot_hi) - float(r.boot_hi)) < 1e-9) and int(o.n_marker) == int(r.n_marker)
    if not same: mm += 1; print("   불일치", r.cell, r.target, r.zone, o.ratio, r.ratio, o.boot_lo, r.boot_lo, o.boot_hi, r.boot_hi)
print("① 구역 표 재현 검산 (옛 구역 · 기존 RT 의 K562·GM12878 행) : 비교 %d 행 · 어긋남 %d" % (nn, mm))
# ---------------- ① the new marker table and the new count table
MKN = {c: L.build_markers(c, "new") for c in ("GM12878", "K562")}
L.write_marker_table("GM12878", MKN["GM12878"], "P11_GM12878_MARKERS_NEWZONE_2026-10-13.tsv")
L.write_marker_table("K562", MKN["K562"], "P11_K562_MARKERS_NEWZONE_2026-10-13.tsv")
GMmix = read_gm("READLEVEL_GM_MIX.tsv"); KQ = L.rd("K562_QUALGATE_MIX.tsv")
qd = collections.defaultdict(dict)
for r in KQ: qd[r["label"]][r["variant"]] = r
# the new count table, in which only the zone part of the tag is new, kept for tracing
with io.open(os.path.join(TO, "K562_QUALGATE_NEWZONE_2026-10-13.tsv"), "w", encoding="utf-8", newline="\n") as f:
    f.write(TAB.join(["label", "variant", "pos", "alt", "n_unit", "tag", "alt_all", "dep_all", "alt_hq", "dep_hq"]) + "\n")
    for lab in qd:
        for v, r in qd[lab].items():
            if v not in MKN["K562"]: continue
            m = MKN["K562"][v]
            f.write(TAB.join([lab, v, r["pos"], r["alt"], str(len(m["units"])), L.ZLAB[m["zone"]] + "|" + m["element"], r["alt_all"], r["dep_all"], r["alt_hq"], r["dep_hq"]]) + "\n")
miss = [v for v in MKN["K562"] if v not in qd["WGSNAT"]] + [v for v in MKN["GM12878"] if v not in GMmix["HIFI"]]
print("② 새 표지 중 계수가 없는 것 :", miss if miss else "없음 (GM12878 %d · K562 %d 모두 계수 있음)" % (len(MKN["GM12878"]), len(MKN["K562"])))
Mnew = {c: {v: r["zone"] for v, r in MKN[c].items()} for c in MKN}
spec_new = [("K562", read_k562_pre("K562_QUALGATE_MIX.tsv", "QUALGATE"), Mnew["K562"], "WGSNAT", ["TBP", "POLR3A", "GTF3C2", "ATAC", "SNAPC5"]),
            ("GM12878", GMmix, Mnew["GM12878"], "HIFI", ["TBP", "POLR3G", "ATAC", "DNase", "INPUT"])]
R_new = rt_loop(spec_new)
R_new.to_csv(os.path.join(TO, "P11_percellline_NEWZONE_bootstrap_2026-10-13.tsv"), sep=TAB, index=False, encoding="utf-8")
print("   새 구역 표 %d 행 → P11_percellline_NEWZONE_bootstrap_2026-10-13.tsv" % len(R_new))
for cell in ("K562", "GM12878"):
    for t in ("ATAC", "TBP", "POLR3A", "SNAPC5"):
        a = R_new[(R_new.cell == cell) & (R_new.target == t)].set_index("zone"); b = OLDRT[(OLDRT.cell == cell) & (OLDRT.target == t)].set_index("zone")
        if "근위" in a.index and "가운데" in a.index and "근위" in b.index:
            print("   %-8s %-7s 근위/가운데  옛 %6.2f → 새 %6.2f  (근위 표지 %d→%d · 가운데 %d→%d)" % (
                cell, t, float(b.loc["근위", "ratio"])/float(b.loc["가운데", "ratio"]), float(a.loc["근위", "ratio"])/float(a.loc["가운데", "ratio"]),
                int(b.loc["근위", "n_marker"]), int(a.loc["근위", "n_marker"]), int(b.loc["가운데", "n_marker"]), int(a.loc["가운데", "n_marker"])))
# ---------------- ③ the interaction against input, inter() unchanged
def build_qdf(qrows_src, zone_of):
    rows = []
    qrec = collections.defaultdict(dict)
    for r in qrows_src: qrec[r["label"]][r["variant"]] = (int(r["alt_hq"]), int(r["dep_hq"]))
    QW = qrec["WGSNAT"]
    for lab in qrec:
        if lab == "WGSNAT": continue
        for v, (a, d) in qrec[lab].items():
            wa, wd = QW.get(v, (0, 0))
            if wd == 0 or d == 0: continue
            q = wa / wd
            if not (0 < q < 1): continue
            if v not in zone_of: continue
            rows.append(dict(track=lab, variant=v, zone=zone_of[v], alt=a, exp=d * q))
    return pd.DataFrame(rows)
def inter_all(qdf, labs, seed=20261008, n=2000):
    R2 = np.random.default_rng(seed); out = {}
    for lab in labs:
        A = qdf[(qdf.track == lab) & (qdf.zone == "근위(합)")]; B = qdf[(qdf.track == lab) & (qdf.zone == "가운데")]
        C = qdf[(qdf.track == "INPUT") & (qdf.zone == "근위(합)")]; Dd = qdf[(qdf.track == "INPUT") & (qdf.zone == "가운데")]
        vz = sorted(set(A.variant) & set(C.variant)); vm = sorted(set(B.variant) & set(Dd.variant))
        if len(vz) < 2 or len(vm) < 2: out[lab] = (np.nan, np.nan, np.nan, len(vz), len(vm)); continue
        Ai, Bi, Ci, Di = A.set_index("variant"), B.set_index("variant"), C.set_index("variant"), Dd.set_index("variant")
        def one(pz, pm):
            a, b, c, d = Ai.loc[pz], Bi.loc[pm], Ci.loc[pz], Di.loc[pm]
            na = a.alt.sum() / max(a.exp.sum(), 1e-12); nb = b.alt.sum() / max(b.exp.sum(), 1e-12)
            nc = c.alt.sum() / max(c.exp.sum(), 1e-12); nd = d.alt.sum() / max(d.exp.sum(), 1e-12)
            if nb <= 0 or nc <= 0 or nd <= 0: return np.nan
            return (na / nb) / (nc / nd)
        pt = one(vz, vm)
        res = [one(R2.choice(vz, len(vz), replace=True), R2.choice(vm, len(vm), replace=True)) for _ in range(n)]
        res = [x for x in res if x == x and np.isfinite(x)]
        out[lab] = (pt,) + (tuple(np.percentile(res, [2.5, 97.5])) if len(res) >= 100 else (np.nan, np.nan)) + (len(vz), len(vm))
    return out
BT = pd.read_csv(os.path.join(T, "P11_K562_ASB_main_v4_with_caveats.tsv"), sep=TAB).set_index("target")
# old: K562_QUALGATE with 75 variants and the zone taken from the tag
qold = L.rd("K562_QUALGATE.tsv"); zold = {r["variant"]: r["tag"].split("|")[0] for r in qold}
LABS_OLD = ("POLR3A", "TBP", "BDP1", "BRF1", "GTF3C2", "ATAC", "SNAPC5")
I_old = inter_all(build_qdf(qold, zold), LABS_OLD)
bad = 0
print("③ 교호작용 재현 검산 (옛 구역 · 기존 본표 nc 열)")
for lab in LABS_OLD:
    pt, lo, hi, nz, nm_ = I_old[lab]; b = BT.loc[lab]
    same = abs(pt - float(b.nc_interaction_OR)) < 0.006 and abs(lo - float(b.nc_lo)) < 0.006 and abs(hi - float(b.nc_hi)) < 0.006
    bad += (not same)
    print("   %-7s 재현 %.2f [%.2f-%.2f] ↔ 본표 %.2f [%.2f-%.2f]  %s" % (lab, pt, lo, hi, float(b.nc_interaction_OR), float(b.nc_lo), float(b.nc_hi), "✓" if same else "불일치"))
print("   ⇒ 어긋남 %d" % bad)
# new: the new 70 markers on the new zones. BDP1 and BRF1 are absent because they were not counted, and the random number stream differs by that much.
znew = {v: L.ZLAB[r["zone"]] for v, r in MKN["K562"].items()}
LABS_NEW = ("POLR3A", "TBP", "GTF3C2", "ATAC", "SNAPC5")
I_new = inter_all(build_qdf(KQ, znew), LABS_NEW)
rows = []
for lab in LABS_NEW:
    pt, lo, hi, nz, nm_ = I_new[lab]; old = I_old[lab]
    print("   %-7s 새 %.2f [%.2f-%.2f] (근위 %d·가운데 %d 표지)   ← 옛 %.2f [%.2f-%.2f]" % (lab, pt, lo, hi, nz, nm_, old[0], old[1], old[2]))
    rows.append(TAB.join([lab, "%.2f" % pt, "%.2f" % lo, "%.2f" % hi, str(nz), str(nm_), "표적특이" if lo > 1 else ("음성대조" if lab == "SNAPC5" else "구분 안 됨")]))
io.open(os.path.join(TO, "K562_interaction_NEWZONE_2026-10-13.tsv"), "w", encoding="utf-8", newline="\n").write(
    TAB.join(["target", "nc_interaction_OR", "nc_lo", "nc_hi", "근위표지", "가운데표지", "판정"]) + "\n" + "\n".join(rows) + "\n")
print("→ K562_interaction_NEWZONE_2026-10-13.tsv")
# ---------------- a summary of the new marker table
for c in ("GM12878", "K562"):
    print("새 표지 %-8s %s" % (c, dict(collections.Counter(r["zone"] for r in MKN[c].values()))))
