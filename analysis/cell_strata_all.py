# -*- coding: utf-8 -*-
"""층(유전자부 β · 고정 구역 · 실측 저메틸 구역)별 두 자리 footprint (2026-10-12)
usage: python cell_strata.py HepG2|K562|GM12878 gene|zonefix|zoneemp
층마다 : 블라인드(10bp · 격자 전체) · 호출기 없음 확인검정 · 두 자리 배수 · 곡선
"""
import io
import os
import sys
import json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cell_sig_lib as S

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
_SELF = os.path.basename(os.path.abspath(__file__))


def _need(k):
    # No defaults are set: if a variable is missing the script stops with a readable message
    v = os.environ.get(k, "")
    if not v:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s." % (k, _SELF)
            + chr(10)
            + "The deposited copy has no default; see README.md section 2")
    return v


DS = {"HepG2": (_need("RDNA_HEPG2_DIR"), "HEPG2"),
      "K562": (_need("RDNA_CELLS_DIR"), "K562r"),
      "GM12878": (_need("RDNA_CELLS_DIR"), "GM12878FS"),
      "K562all": (_need("RDNA_CELLS_DIR"), "K562all"),
      "GM12878all": (_need("RDNA_CELLS_DIR"), "GM12878all")}
NAME, SCHEME = sys.argv[1], sys.argv[2]
DIR, LAB = DS[NAME]
HALF, NP, PER = S.HALF, S.NP, S.PER
rng = np.random.default_rng(20261013)
NB_BLIND, NB_CONF, NB_CURVE = 500, 1000, 200
MINEV = 15

keys = []
with io.open(os.path.join(DIR, "UNIT_%s.tsv" % LAB), encoding="utf8") as f:
    f.readline()
    for l in f:
        keys.append(l.split(TB, 3)[2])
N = len(keys)
ix = dict((k, i) for i, k in enumerate(keys))
mi, ce, sz = [], [], []
with io.open(os.path.join(DIR, "FHMM_UNITS_%s.tsv" % LAB), encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(NL).split(TB)
        if x[4] != "tf":
            continue
        c = (int(x[5]) + int(x[6])) // 2
        if -HALF <= c < HALF and x[2] in ix:
            mi.append(ix[x[2]]); ce.append(c + HALF); sz.append(int(x[7]))
mi, ce, sz = np.array(mi), np.array(ce), np.array(sz)
A = np.zeros((N, NP), np.int16)
M = np.zeros((N, NP), np.int16)
with io.open(os.path.join(DIR, "SITE_M6A_%s.tsv" % LAB), encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(NL).split(TB)
        i = ix.get(x[1])
        if i is None:
            continue
        for v in (x[3].split(",") if len(x) > 3 and x[3] else []):
            r = int(v)
            if S.LO <= r <= S.HI:
                A[i, r - S.LO] += 1
        for v in (x[4].split(",") if len(x) > 4 and x[4] else []):
            r = int(v)
            if S.LO <= r <= S.HI:
                M[i, r - S.LO] += 1

lab = {}
with io.open(os.path.join(DIR, "strata", "%s_%s.tsv" % (NAME, SCHEME)), encoding="utf8") as f:
    f.readline()
    for l in f:
        k, v = l.rstrip(NL).split(TB)
        if k in ix:
            lab[ix[k]] = v
names = sorted(set(lab.values()))
ORDER = {"gene": ["low", "high"], "zonefix": ["원위", "가운데", "근위"]}
if SCHEME in ORDER:
    names = [n for n in ORDER[SCHEME] if n in names]
print("%s %s : 층 %s" % (NAME, SCHEME, [(n, sum(1 for v in lab.values() if v == n)) for n in names]))
R = {"name": NAME, "scheme": SCHEME, "strata": []}
arr = {"rel": S.REL}
for si, nm in enumerate(names):
    idx = np.array(sorted(i for i, v in lab.items() if v == nm))
    n = len(idx)
    remap = -np.ones(N, int)
    remap[idx] = np.arange(n)
    sel = remap[mi] >= 0
    smi, sce, ssz = remap[mi[sel]], ce[sel], sz[sel]
    res = dict(label=nm, n=int(n), n_events=int(sel.sum()))
    flo = 1000.0 / (n * S.SM)
    # 곡선
    for k, (a0, a1) in enumerate([(20, 29), (30, 39), (40, 49), (50, 59)]):
        s = (ssz >= a0) & (ssz <= a1)
        o, b = S.dens_fold(smi[s], sce[s], n, NB_CURVE, rng)
        arr["s%d_c%d_obs" % (si, k)] = o
        arr["s%d_c%d_base" % (si, k)] = b
    for bnm, (a0, a1), cpos in (("b35", (30, 39), -11), ("b40", (40, 49), 62)):
        s = (ssz >= a0) & (ssz <= a1)
        o, b = S.dens_fold(smi[s], sce[s], n, NB_CURVE, rng)
        arr["s%d_%s_obs" % (si, bnm)] = o
        arr["s%d_%s_base" % (si, bnm)] = b
        j = int(np.where(S.REL == cpos)[0][0])
        res[bnm] = dict(center=cpos, fold=float(o[j] / max(b[j], flo)),
                        n_events=int(((sce[s] - HALF >= cpos - 5) & (sce[s] - HALF <= cpos + 5)).sum()))
    sA = A[idx].sum(0).astype(float)
    sM = M[idx].sum(0).astype(float)
    res["dub_m11"] = S.dub_paper(sA, sM, -11)
    res["dub_p58"] = S.dub_paper(sA, sM, 58)
    # 블라인드
    bl = S.blind(smi, sce, ssz, n, NB_BLIND, rng, minev=MINEV, sbin=10)
    if bl is None:
        res["blind"] = None
    else:
        res["blind"] = dict(nk=bl["nk"], thr=bl["thr"], NB=bl["NB"], sites=bl["sites"],
                            peaks=[list(p[:2]) + list(p[2]) for p in bl["peaks"][:6]])
    # 호출기 없음 확인검정
    okm = A[idx].sum(1) >= 20
    if okm.sum() >= 100:
        res["confirm"] = S.confirm(A[idx][okm], M[idx][okm], NB_CONF, rng)
        res["confirm_n"] = int(okm.sum())
    else:
        res["confirm"] = None
    R["strata"].append(res)
    bs = res["blind"]["sites"] if res["blind"] else None
    print("  [%s] n=%d 사건 %d · TFIIIB %.1f×(사건 %d) ICR %.1f×(사건 %d) · flank/core %.2f/%.2f"
          % (nm, n, res["n_events"], res["b35"]["fold"], res["b35"]["n_events"],
             res["b40"]["fold"], res["b40"]["n_events"], res["dub_m11"], res["dub_p58"]))
    if bs:
        print("      블라인드 칸 %d · p_격자 TFIIIB %.4f ICR %.4f · p_창 %.4f / %.4f"
              % (res["blind"]["nk"], bs["TFIIIB"]["p_grid"], bs["ICR"]["p_grid"],
                 bs["TFIIIB"]["p_win"], bs["ICR"]["p_win"]))
    else:
        print("      블라인드 : 사건 %d 문턱을 넘는 칸 없음 -> 검정 불가" % MINEV)
    if res["confirm"]:
        print("      확인(호출기 없음) p TFIIIB %.4f ICR %.4f (분자 %d)"
              % (res["confirm"]["TFIIIB"]["p"], res["confirm"]["ICR"]["p"], res["confirm_n"]))
out = os.path.join(DIR, "strata_%s_%s" % (NAME, SCHEME))
np.savez_compressed(out + ".npz", **arr)
io.open(out + ".json", "w", encoding="utf8").write(
    json.dumps(R, ensure_ascii=False, indent=1,
               default=lambda o: float(o) if hasattr(o, "item") else str(o)))
print("저장 ->", out)
