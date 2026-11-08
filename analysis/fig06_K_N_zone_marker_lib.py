# -*- coding: utf-8 -*-
# [archive copy] from the working copy fig6v3_lib.py (md5 4f0e1a56926a) · data and output roots via environment variables with no default (the script stops if one is unset) · module names changed to the archive names · no change to any computation
import os as _os_need, sys as _sys_need
def _need(k):
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set; it is required by %s (a module stops the script that imports it). The deposited copy has no default; see README.md section 2" % (k, _os_need.path.basename(globals().get("__file__", "an exec-loaded block"))))
    return v
# . Figure 6 main figure version 3: the shared module for the new zone definition, which is the low methylation span of the gene body.
#   New zones = for each haplotype and each end, "units 1 to k counted from that end", where k is the change point estimate in 02_tables/lowmeth_changepoint_4haplotypes_2026-10-12.tsv. The rest is the interior.
#   There are two frames. The marker table uses the canonical 110 unit frame for K562 hapB, while the methylation, accessibility and occupancy measurements use the earlier 98 unit frame. The boundary is set separately in each frame.
#     K562 hapB: marker frame proximal 19 and distal 14, measurement frame proximal 18 and distal 14. The distal number is the same in both and the proximal differs because the 110 unit frame has one more unit at that end (hapb_frame_align.py).
import sys, io, os, re, csv, collections
import numpy as np
T = _need("P11_TABLE_DIR"); TAB = chr(9)
V = _need("V11PIPE_DIR")
NBD = os.environ.get("NB81_DIR")   # not used by any archived code path; kept as a name only, no default        # the read level and position level tables of , read only
def rd(fn, d=T): return list(csv.DictReader(io.open(os.path.join(d, fn), encoding="utf-8"), delimiter=TAB))
_CP = {(r["집합"], r["끝"]): r for r in rd("lowmeth_changepoint_4haplotypes_2026-10-12.tsv")}
def _k(nm, end, key="k_hat(해당 틀)"): return int(_CP[(nm, end)][key])
# (cell line, slot) -> {frame: (proximal k, distal k, N)}
KB = {("GM12878", "1"): dict(marker=(_k("GM12878 hap1", "근위"), _k("GM12878 hap1", "원위"), 119), meas=(_k("GM12878 hap1", "근위"), _k("GM12878 hap1", "원위"), 119)),
      ("GM12878", "2"): dict(marker=(_k("GM12878 hap2", "근위"), _k("GM12878 hap2", "원위"), 91), meas=(_k("GM12878 hap2", "근위"), _k("GM12878 hap2", "원위"), 91)),
      ("K562", "1"): dict(marker=(_k("K562 hapA", "근위"), _k("K562 hapA", "원위"), 139), meas=(_k("K562 hapA", "근위"), _k("K562 hapA", "원위"), 139)),
      ("K562", "2"): dict(marker=(_k("K562 hapB (98 틀)", "근위", ) + 1, _k("K562 hapB (98 틀)", "원위"), 110), meas=(_k("K562 hapB (98 틀)", "근위"), _k("K562 hapB (98 틀)", "원위"), 98))}
# Check: the 110 unit frame values equal the 110 frame k_hat column of the change point table
assert KB[("K562", "2")]["marker"][0] == int(_CP[("K562 hapB (98 틀)", "근위")]["110틀_k_hat"]) and KB[("K562", "2")]["marker"][1] == int(_CP[("K562 hapB (98 틀)", "원위")]["110틀_k_hat"])
HAPNAME = {("K562", "1"): "hapA", ("K562", "2"): "hapB", ("GM12878", "1"): "hap1", ("GM12878", "2"): "hap2"}
def zone_new(cell, slot, ap, frame="marker"):
    kp, kd, N = KB[(cell, str(slot))][frame]
    if ap <= kd: return "원위"
    if N - ap + 1 <= kp: return "근위"
    return "가운데"
def zone_old(cell, slot, ap, N):
    return "원위" if ap <= 20 else ("근위" if N - ap + 1 <= 20 else "가운데")
isS1 = lambda v: bool(re.fullmatch(r"\d+[ACGT]", v or ""))
def elem(p):
    return ("AluY" if p <= 318 else "NTS1" if p <= 1178 else "GA-rich" if p <= 1551 else
            "NTS2" if p <= 2031 else "상류(5'-flank)" if p <= 2121 else "유전자")
def _rows():
    gm = [dict(slot=r["slot"], N=int(r["N"]), ap=int(r["array_pos"]), varlist=r["varlist"]) for r in rd("MAP_G8_GM12878FS.tsv")
          if r["sample"] == "GM12878FS" and r["seg_type"] == "UNIT" and str(r["array_pos"]).isdigit()]
    g8 = [r for r in csv.DictReader(io.open(os.path.join(V, "MAP_G8_K562r.tsv"), encoding="utf-8"), delimiter=TAB)]
    ka = [dict(slot="1", N=int(r["N"]), ap=int(r["array_pos"]), varlist=r["varlist"]) for r in g8
          if r["sample"] == "K562r" and r["seg_type"] == "UNIT" and str(r["array_pos"]).isdigit() and r["slot"] == "1"]
    nb = [r for r in csv.DictReader(io.open(os.path.join(_need("P11_ASSEMBLY_DIR"), "MAP_hapB_v4_110units.tsv"), encoding="utf-8"), delimiter=TAB)
          if r["sample"] == "K562Ex" and r["seg_type"] == "UNIT" and str(r["array_pos"]).isdigit()]
    kb = [dict(slot="2", N=len(nb), ap=int(r["array_pos"]), varlist=r["varlist"]) for r in nb]
    return {"GM12878": gm, "K562": ka + kb}
ROWS = _rows()
DROP = {"GM12878": set(), "K562": {"1922A", "1927T", "809A", "875T"}}      # the four dropped by the gate when K562 v4 became v4_final
def build_markers(cell, mode):
    """mode 'old' = proximal and distal units 1 to 20. mode 'new' = the low methylation span. Returns {variant: dict(zone, element, units[list of 'slot:ap'])}"""
    zf = (lambda sl, ap, N: zone_old(cell, sl, ap, N)) if mode == "old" else (lambda sl, ap, N: zone_new(cell, sl, ap, "marker"))
    car = collections.defaultdict(list)
    for r in ROWS[cell]:
        z = zf(r["slot"], r["ap"], r["N"]); key = "%s:%d" % (r["slot"], r["ap"])
        for v in str(r["varlist"]).split(","):
            if not isS1(v): continue
            p = int(re.match(r"(\d+)", v).group(1))
            if p <= 318 or 1179 <= p <= 1551: continue
            car[v].append((key, z))
    nm = collections.Counter(u for c in car.values() for u, _ in c)
    per = [nm.get("%s:%d" % (r["slot"], r["ap"]), 0) for r in ROWS[cell]]
    cut = max(12, int(np.percentile(per, 99)) + 1); HOT = {u for u, n in nm.items() if n >= cut}
    out = {}
    for v, lst in car.items():
        if v in DROP[cell]: continue
        keep = [(u, z) for u, z in lst if u not in HOT]
        if not keep: continue
        zs = {z for _, z in keep}
        if len(zs) != 1: continue
        p = int(re.match(r"(\d+)", v).group(1))
        out[v] = dict(zone=list(zs)[0], element=elem(p), pos=p, alt=v[-1], units=[u for u, _ in keep])
    return out
ZLAB = {"근위": "근위(합)", "가운데": "가운데", "원위": "원위(합)"}          # the zone labels used by the existing marker table
def write_marker_table(cell, mk, fn):
    rows = sorted(mk.items(), key=lambda kv: (["근위", "가운데", "원위"].index(kv[1]["zone"]), -len(kv[1]["units"]), kv[1]["pos"]))
    TO = _need("P11_TABLE_OUT_DIR"); os.makedirs(TO, exist_ok=True)   # writes go here (no default)
    with io.open(os.path.join(TO, fn), "w", encoding="utf-8", newline="\n") as f:
        f.write(TAB.join(["zone", "variant", "pos", "alt", "element", "n_units", "units"]) + "\n")
        for v, r in rows:
            f.write(TAB.join([ZLAB[r["zone"]], v, str(r["pos"]), r["alt"], r["element"], str(len(r["units"])), ";".join(r["units"])]) + "\n")
