# -*- coding: utf-8 -*-
# [archive copy] from the working copy s2_newzone_values.py (md5 15e5a60896bf) · data and output roots via environment variables with no default (the script stops if one is unset) · module names changed to the archive names · no change to any computation
# WARNING — this script WRITES tables, into P11_TABLE_OUT_DIR (no default). Pointing it at the folder that holds the published
#   input tables (P11_TABLE_DIR) OVERWRITES them. To verify, set P11_TABLE_OUT_DIR to an empty temporary folder.
import os as _os_need, sys as _sys_need
def _need(k):
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set; it is required by %s (a module stops the script that imports it). The deposited copy has no default; see README.md section 2" % (k, _os_need.path.basename(globals().get("__file__", "an exec-loaded block"))))
    return v
# . Supplementary S2, which holds the ChIP targets left out of the main figure, is recomputed on the new zones, that is, the low methylation span of the gene body. This is the value step.
#   ① K562 BDP1 and BRF1 come from the node table K562_QUALGATE_MIX_BDP1BRF1.tsv (mixk_bdp1brf1.py, same position list and same gate as K562_QUALGATE_MIX).
#      - The Q >= 20 zone ratio uses the same formula as agg() in fig6_main_v8.py. First the main figure h values are reproduced with POLR3A, TBP, GTF3C2 and SNAPC5.
#      - The interaction against input is build_qdf and inter_all of fig6v3_step1_marker_values.py, unchanged (a copy), with the same seed.
#        The order of the first five targets is kept and BDP1 and BRF1 are appended, so the first five must not differ in a single cell from K562_interaction_NEWZONE_2026-10-13 (checked).
#   ② GM12878 ChIP is read from the new zone table with the NA12878 WGS denominator, P11_percellline_NEWZONE_GMWGS_bootstrap_2026-10-14.
#      - TBP needs leave-one-out, because one cell carries a large share of the markers (the cancellation of 2056C and 2079C is recorded), and it uses the same denominator (WGS_ERR3239334).
#   Only new files are written.
import sys, io, os, collections
import numpy as np, pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
T = _need("P11_TABLE_DIR"); TAB = chr(9)
TO = _need("P11_TABLE_OUT_DIR"); os.makedirs(TO, exist_ok=True)   # writes go here (no default)
def rd(fn): return pd.read_csv(os.path.join(T, fn), sep=TAB)
def rdl(fn): return [dict(r) for _, r in rd(fn).astype(str).iterrows()]
MK = rd("P11_K562_MARKERS_NEWZONE_2026-10-13.tsv"); ZOF = dict(zip(MK.variant, MK.zone))          # the pooled proximal, the interior and the pooled distal
Z3L = ["근위(합)", "가운데", "원위(합)"]; ZS = {"근위(합)": "proximal", "가운데": "interior", "원위(합)": "distal"}
KQ = rdl("K562_QUALGATE_MIX.tsv"); KB = rdl("K562_QUALGATE_MIX_BDP1BRF1.tsv")
assert set(r["label"] for r in KB) == {"BDP1", "BRF1"}
assert sorted(r["variant"] for r in KB if r["label"] == "BDP1") == sorted(r["variant"] for r in KQ if r["label"] == "WGSNAT"), "자리표가 다르다"
ALL = KQ + KB
# ---------------- ① the Q >= 20 zone ratio, the same formula as fig6_main_v8.agg
W = {r["variant"]: r for r in ALL if r["label"] == "WGSNAT"}
def agg(lab, z):
    O = Ex = 0.0; n = 0
    for r in ALL:
        if r["label"] != lab or ZOF.get(r["variant"]) != z or r["variant"] not in W: continue
        w = W[r["variant"]]; a, d, wa, wd = int(r["alt_hq"]), int(r["dep_hq"]), int(w["alt_hq"]), int(w["dep_hq"])
        if wd == 0 or d == 0 or wa == 0: continue
        O += a; Ex += d * (wa / wd); n += 1
    return O, Ex, n
H = {"POLR3A": (5.418, 0.113), "TBP": (5.911, 0.212), "GTF3C2": (5.976, 0.188), "SNAPC5": (0.899, 0.336)}   # the main figure h values, used for the reproduction check only and not drawn
bad = 0
for lab, (hn, hm) in H.items():
    o1, e1, _ = agg(lab, "근위(합)"); o2, e2, _ = agg(lab, "가운데")
    ok = abs(o1 / e1 - hn) < 0.0006 and abs(o2 / e2 - hm) < 0.0006; bad += not ok
    print("  재현 %-7s 근위 %.4f 가운데 %.4f  ↔ 본그림 h %.3f %.3f  %s" % (lab, o1 / e1, o2 / e2, hn, hm, "✓" if ok else "불일치"))
print("  ⇒ h 재현 어긋남 %d" % bad); assert bad == 0
zrows = []
for lab in ("BDP1", "BRF1", "POLR3A", "TBP", "GTF3C2", "SNAPC5", "INPUT"):
    for z in Z3L:
        o, e, n = agg(lab, z)
        zrows.append(dict(target=lab, zone=ZS[z], n_marker=n, sum_alt_hq=int(o), sum_exp_hq=round(e, 3), ratio=round(o / e, 4) if e else ""))
ZR = pd.DataFrame(zrows)
# ---------------- ① the interaction, a copy of build_qdf and inter_all from fig6v3_step1_marker_values.py
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
LABS = ("POLR3A", "TBP", "GTF3C2", "ATAC", "SNAPC5", "BDP1", "BRF1")
I = inter_all(build_qdf(ALL, ZOF), LABS)
BT = rd("K562_interaction_NEWZONE_2026-10-13.tsv").set_index("target"); bad = 0
for lab in LABS[:5]:
    pt, lo, hi, nz, nm = I[lab]; b = BT.loc[lab]
    ok = ("%.2f %.2f %.2f %d %d" % (pt, lo, hi, nz, nm)) == ("%.2f %.2f %.2f %d %d" % (b.nc_interaction_OR, b.nc_lo, b.nc_hi, b["근위표지"], b["가운데표지"])); bad += not ok
    print("  교호작용 재현 %-7s %.2f [%.2f-%.2f] (%d·%d) ↔ 본표 %.2f [%.2f-%.2f]  %s" % (lab, pt, lo, hi, nz, nm, b.nc_interaction_OR, b.nc_lo, b.nc_hi, "✓" if ok else "불일치"))
print("  ⇒ 교호작용 재현 어긋남 %d" % bad); assert bad == 0
rows = []
for lab in LABS:
    pt, lo, hi, nz, nm = I[lab]
    rows.append(dict(target=lab, nc_interaction_OR=round(pt, 2), nc_lo=round(lo, 2), nc_hi=round(hi, 2), 근위표지=nz, 가운데표지=nm,
                     판정=("표적특이" if lo > 1 else ("음성대조" if lab == "SNAPC5" else "판정 불가(CI 가 1 을 문다)"))))
IT = pd.DataFrame(rows)
IT.to_csv(os.path.join(TO, "K562_interaction_NEWZONE_withBDP1BRF1_2026-10-14.tsv"), sep=TAB, index=False, encoding="utf-8")
ZR.to_csv(os.path.join(TO, "K562_zoneratio_Q20_NEWZONE_withBDP1BRF1_2026-10-14.tsv"), sep=TAB, index=False, encoding="utf-8")
for lab in ("BDP1", "BRF1"):
    s = ZR[ZR.target == lab].set_index("zone"); r = IT.set_index("target").loc[lab]
    print("  %-5s 근위 %.3f (표지 %d) · 가운데 %.3f (표지 %d) · 원위 %s (표지 %d) · 근위/가운데 %.1f배 · 교호작용 %.2f [%.2f-%.2f] → %s" % (
        lab, s.ratio["근위"], s.n_marker["근위"], s.ratio["가운데"], s.n_marker["가운데"], s.ratio["원위"], s.n_marker["원위"],
        s.ratio["근위"] / s.ratio["가운데"], r.nc_interaction_OR, r.nc_lo, r.nc_hi, r.판정))
# ---------------- ② GM12878 TBP leave-one-out on the new zones with the WGS denominator, using the same counts as the RT, that is, READLEVEL_GM_MIX filtered alt_n over depth
# Markers that drop out of the WGS denominator, meaning that genome carries no alternate allele and q is 0, are recorded as well. 2079C of the earlier HiFi denominator version is such a case.
GM = rd("P11_GM12878_MARKERS_NEWZONE_2026-10-13.tsv"); GZ = dict(zip(GM.variant, GM.zone))
RL = rd("READLEVEL_GM_MIX.tsv"); RL = RL[RL["filter"] == "filtered"]
DEN = pd.read_csv(os.path.join(T, "denominator_comparison_2026-10-14", "denom_GM_WGS_ERR3239334.tsv"), sep=TAB)
WD = {v: (int(a), int(d)) for v, a, d in zip(DEN.variant, DEN.alt_n, DEN.depth)}
RT = rd("P11_percellline_NEWZONE_GMWGS_bootstrap_2026-10-14.tsv")
loo = []
for z in Z3L:
    vs = []
    for _, r in RL[RL.label == "TBP"].iterrows():
        if GZ.get(r.variant) != z: continue
        wa, wn = WD.get(r.variant, (0, 0))
        if wn == 0 or int(r.depth) == 0 or wa == 0: continue
        vs.append((r.variant, int(r.alt_n), int(r.depth) * wa / wn))
    O = sum(x[1] for x in vs); E = sum(x[2] for x in vs)
    rt = RT[(RT.cell == "GM12878") & (RT.target == "TBP") & (RT.zone == ZS[z])]
    assert len(rt) == 1 and abs(float(rt.ratio.iloc[0]) - O / E) < 0.0006 and int(rt.n_marker.iloc[0]) == len(vs), ("RT 재현 실패", z)
    pick = {x[0] for x in sorted(vs, key=lambda x: -x[1])[:3]} | {x[0] for x in sorted(vs, key=lambda x: -x[2])[:3]}   # the top 3 by observed count joined with the top 3 by expected count
    for v, o, e in sorted([x for x in vs if x[0] in pick], key=lambda x: -x[1]):
        loo.append(dict(zone=ZS[z], 전체비=round(O / E, 4), 표지수=len(vs), 뺀표지=v, 그표지_관측=o, 그표지_기대=round(e, 2), 기대중몫=round(e / E, 3),
                        그표지_비=round(o / e, 3) if e else "", 관측중몫=round(o / O, 3) if O else "", 빼면=round((O - o) / (E - e), 4) if E - e > 0 else ""))
LO = pd.DataFrame(loo)
for v in ("2056C", "2079C", "2061G"):
    t_ = RL[(RL.label == "TBP") & (RL.variant == v)]; print("  %s : 구역 %s · TBP %s/%s · WGS 분모 %s/%s" % (v, GZ.get(v), t_.alt_n.iloc[0] if len(t_) else "-", t_.depth.iloc[0] if len(t_) else "-", *WD.get(v, ("-", "-"))))
LO.to_csv(os.path.join(TO, "GM12878_TBP_LOO_NEWZONE_GMWGS_2026-10-14.tsv"), sep=TAB, index=False, encoding="utf-8")
print("  GM12878 TBP LOO (RT 재현 ✓)"); print(LO.to_string(index=False))
print("  → K562_interaction_NEWZONE_withBDP1BRF1_2026-10-14.tsv · K562_zoneratio_Q20_NEWZONE_withBDP1BRF1_2026-10-14.tsv · GM12878_TBP_LOO_NEWZONE_GMWGS_2026-10-14.tsv")
# ---------------- ③ the GM12878 interaction against input, inter_all unchanged, counts = READLEVEL_GM_MIX filtered, q = NA12878 WGS (ERR3239334)
#   GM12878 has no negative control ChIP of the SNAPC5 kind. The INPUT here is an untargeted control, and the reason it is not put in the main figure is unchanged.
gq = []
for _, r in RL.iterrows():
    if r.label not in ("TBP", "POLR3G", "INPUT", "ATAC") or r.variant not in GZ: continue
    wa, wn = WD.get(r.variant, (0, 0))
    if wn == 0 or int(r.depth) == 0: continue
    q = wa / wn
    if not (0 < q < 1): continue
    gq.append(dict(track=r.label, variant=r.variant, zone=GZ[r.variant], alt=int(r.alt_n), exp=int(r.depth) * q))
IG = inter_all(pd.DataFrame(gq), ("TBP", "POLR3G", "ATAC"))
grow = []
for lab in ("TBP", "POLR3G", "ATAC"):
    pt, lo, hi, nz, nm = IG[lab]
    grow.append(dict(target=lab, nc_interaction_OR=round(pt, 2), nc_lo=round(lo, 2), nc_hi=round(hi, 2), 근위표지=nz, 가운데표지=nm))
    print("  GM12878 %-6s 교호작용(INPUT 기준) %.2f [%.2f-%.2f] (근위 %d · 가운데 %d 표지)" % (lab, pt, lo, hi, nz, nm))
pd.DataFrame(grow).to_csv(os.path.join(TO, "GM12878_interaction_NEWZONE_GMWGS_2026-10-14.tsv"), sep=TAB, index=False, encoding="utf-8")
gi = RT[(RT.cell == "GM12878") & (RT.target == "INPUT")].set_index("zone")
print("  GM12878 INPUT 근위 관측 %d · 기대 %.2f · 가운데 관측 %d · 기대 %.2f" % (gi.sum_alt["근위"], gi.sum_exp["근위"], gi.sum_alt["가운데"], gi.sum_exp["가운데"]))
print("  → GM12878_interaction_NEWZONE_GMWGS_2026-10-14.tsv")
