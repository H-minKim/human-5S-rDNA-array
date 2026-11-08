# -*- coding: utf-8 -*-
r"""Ranks per zone separately: layer (all / shared type removed / single-unit type only) x 4 zones

  For each zone compute the null z of that zone, and rank by |z| within the zone
  (earlier only the single 'most distant zone' was used. Here every zone is reported)
Output: REPRO/VTZONE4_byzone.tsv
"""
import io, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CODE = __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir))
sys.path.insert(0, CODE)
import v11_io as G
from vgroup import GROUPS, GI, vgroup, iscpg_factory
import origin_filter as OF

OUT = os.path.join((os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"]), "REPRO")
rng = np.random.default_rng(20260913)
NN = 200
iscpg = iscpg_factory(G.CAN)
U, _, _ = G.load("primary", verbose=False)
HK = sorted(U)
NH = len(HK)
NG = len(GROUPS)
hapf, kcnt, NS, _ = OF.build(U, HK)
LAY = {"전부": {"다개체 공유형", "배열내 다중형", "단일유닛형"},
       "공유형 뺌": {"배열내 다중형", "단일유닛형"},
       "단일유닛형만": {"단일유닛형"}}
LAYS = ["전부", "공유형 뺌", "단일유닛형만"]
ZN = ["근위 1-3", "근위 4-17", "가운데", "원위 1-10"]

AP, AD, HI = [], [], []
V = {l: [] for l in LAYS}
for hi, k in enumerate(HK):
    for row in U[k]:
        AP.append(row[0]); AD.append(row[1]); HI.append(hi)
        c = {l: np.zeros(NG) for l in LAYS}
        for x in row[4]:
            j = GI[vgroup(x, iscpg)]
            kl = OF.klass(x, k, hapf, kcnt)
            for l in LAYS:
                if kl in LAY[l]:
                    c[l][j] += 1
        for l in LAYS:
            V[l].append(c[l])
AP = np.asarray(AP); AD = np.asarray(AD); HI = np.asarray(HI)
V = {l: np.asarray(v) for l, v in V.items()}


def zone_of(ap, ad):
    z = np.full(len(ap), 2, int)
    pz = np.where(ap <= 3, 0, np.where(ap <= 17, 1, -1))
    dz = ad <= 10
    both = (pz >= 0) & dz
    z = np.where(pz >= 0, pz, z)
    z = np.where(dz, 3, z)
    z = np.where(both, np.where(ap <= ad, pz, 3), z)
    return z


def zm(X, zz):
    out = np.zeros((len(ZN), X.shape[1]))
    for i in range(len(ZN)):
        acc = [X[(zz == i) & (HI == hi)].mean(0) for hi in range(NH)
               if ((zz == i) & (HI == hi)).any()]
        out[i] = np.mean(np.asarray(acc), 0)
    return out


Z0 = zone_of(AP, AD)
ORD = {hi: np.where(HI == hi)[0] for hi in range(NH)}
PERM = []
for b in range(NN):
    ap2 = AP.copy(); ad2 = AD.copy()
    for hi, ix in ORD.items():
        q = rng.permutation(len(ix))
        ap2[ix] = AP[ix][q]; ad2[ix] = AD[ix][q]
    PERM.append(zone_of(ap2, ad2))

# produce two rank rulers. The code builds the table to the end (no manual pasting later)
#   within-zone rank = by |null z| (how unexpected)   ,   ratio rank = by |log2 ratio to middle| (size of change)
MIDI4 = ZN.index("가운데")
OUTZ = {}
REC = []
for l in LAYS:
    O = zm(V[l], Z0)
    NUL = np.array([zm(V[l], p) for p in PERM])
    Zz = (O - NUL.mean(0)) / np.where(NUL.std(0) > 0, NUL.std(0), np.nan)
    OUTZ[l] = (O, Zz)
    GG4 = [gg for gg in GROUPS if gg != "기타"]
    for i, zn in enumerate(ZN):
        zr = {gg: rk for rk, gg in enumerate(
            sorted(GG4, key=lambda gg: -abs(np.nan_to_num(Zz[i, GI[gg]]))), 1)}
        lg = {}
        for gg in GG4:
            m = O[MIDI4, GI[gg]]
            v = O[i, GI[gg]]
            lg[gg] = np.log2(v / m) if (m > 0 and v > 0) else np.nan
        rr = {gg: rk for rk, gg in enumerate(
            sorted(GG4, key=lambda gg: -(abs(lg[gg]) if np.isfinite(lg[gg]) else -1)), 1)}
        for gg in GG4:
            REC.append((l, zn, gg, O[i, GI[gg]], Zz[i, GI[gg]], zr[gg], lg[gg], rr[gg]))

g = io.open(os.path.join(OUT, "VTZONE4_byzone.tsv"), "w", encoding="utf8", newline="\n")
g.write("층\t구역\t무리\t유닛당\t귀무z\t구역안_순위\t가운데대비_log2비\t비_순위\n")
for l in LAYS:
    for zn in ZN:
        for r0 in [x for x in REC if x[0] == l and x[1] == zn]:
            g.write("%s\t%s\t%s\t%.5f\t%+.2f\t%d\t%s\t%d\n"
                    % (r0[0], r0[1], r0[2], r0[3], r0[4], r0[5],
                       ("%.4f" % r0[6]) if np.isfinite(r0[6]) else "", r0[7]))
g.close()

for l in ("전부", "단일유닛형만"):
    O, Zz = OUTZ[l]
    print("\n■ %s — 구역마다 상위 6 (그 구역의 귀무 z)" % l)
    for i, zn in enumerate(ZN):
        order = sorted([gg for gg in GROUPS if gg != "기타"],
                       key=lambda gg: -abs(np.nan_to_num(Zz[i, GI[gg]])))[:6]
        print("   %-10s %s" % (zn, " · ".join("%s(%+.0f)" % (gg, Zz[i, GI[gg]])
                                              for gg in order)))
print("\n저장 %s/VTZONE4_byzone.tsv" % OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('VTZONE4_byzone.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
