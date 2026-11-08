# -*- coding: utf-8 -*-
r"""After removing inherited variants (multi-individual shared type), which group changes the most between zones

Layer: environment variable ORIGIN_KEEP
   notshared  = multi-individual shared type removed (within-array multi-copy type plus single-unit type)   <- default
   single     = single-unit type only
   all        = everything
Output: REPRO/VTZONE3_<layer>_zones.tsv
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
os.makedirs(OUT, exist_ok=True)
KEEP = os.environ.get("ORIGIN_KEEP", "notshared")
rng = np.random.default_rng(20260913)
NN = 200
iscpg = iscpg_factory(G.CAN)

U, _, _ = G.load("primary", verbose=False)
HK = sorted(U)
NH = len(HK)
NG = len(GROUPS)
hapf, kcnt, NS, _ = OF.build(U, HK)
print("기원 코호트 : 유닛>=%d 슬롯 %d · fT %.2f · 층 = %s" % (OF.MINU, NS, OF.FT, KEEP))

OKC = {"notshared": {"배열내 다중형", "단일유닛형"},
       "single": {"단일유닛형"},
       "all": {"다개체 공유형", "배열내 다중형", "단일유닛형"}}[KEEP]

VEC, AP, AD, HI = [], [], [], []
drop = kept = 0
for hi, k in enumerate(HK):
    for row in U[k]:
        v = np.zeros(NG)
        for x in row[4]:
            if OF.klass(x, k, hapf, kcnt) in OKC:
                v[GI[vgroup(x, iscpg)]] += 1; kept += 1
            else:
                drop += 1
        VEC.append(v); AP.append(row[0]); AD.append(row[1]); HI.append(hi)
VEC = np.asarray(VEC); AP = np.asarray(AP); AD = np.asarray(AD); HI = np.asarray(HI)
print("변이 출현 : 남김 %d · 뺌 %d (남은 몫 %.4f) · 유닛당 %.3f"
      % (kept, drop, kept / (kept + drop), VEC.sum() / len(VEC)))

ZN = ["근위 1-3", "근위 4-17", "가운데", "원위 1-10"]


def zone_of(ap, ad):
    z = np.full(len(ap), 2, int)
    pz = np.where(ap <= 3, 0, np.where(ap <= 17, 1, -1))
    dz = ad <= 10
    both = (pz >= 0) & dz
    z = np.where(pz >= 0, pz, z)
    z = np.where(dz, 3, z)
    z = np.where(both, np.where(ap <= ad, pz, 3), z)
    return z


def zmeans(zz):
    out = np.full((len(ZN), NG), np.nan); nh = np.zeros(len(ZN), int)
    for i in range(len(ZN)):
        acc = [VEC[(zz == i) & (HI == hi)].mean(0)
               for hi in range(NH) if ((zz == i) & (HI == hi)).any()]
        out[i] = np.mean(np.asarray(acc), 0); nh[i] = len(acc)
    return out, nh


Z0 = zone_of(AP, AD)
OBS, NHZ = zmeans(Z0)
MIDI = ZN.index("가운데")
ORD = {hi: np.where(HI == hi)[0] for hi in range(NH)}
NUL = np.empty((NN, len(ZN), NG))
for b in range(NN):
    ap2 = AP.copy(); ad2 = AD.copy()
    for hi, ix in ORD.items():
        q = rng.permutation(len(ix))
        ap2[ix] = AP[ix][q]; ad2[ix] = AD[ix][q]
    NUL[b] = zmeans(zone_of(ap2, ad2))[0]

g = io.open(os.path.join(OUT, "VTZONE3_%s_zones.tsv" % KEEP), "w", encoding="utf8", newline="\n")
g.write("무리\t" + "\t".join(ZN) + "\t최대_log2비\t최대_절대차\t최대_귀무z\t가장_먼_구역\n")
rows = []
for t in GROUPS:
    j = GI[t]
    mid = OBS[MIDI, j]
    lg = np.array([np.log2(OBS[i, j] / mid) if (mid > 0 and OBS[i, j] > 0) else np.nan
                   for i in range(len(ZN))])
    zz = np.array([(OBS[i, j] - NUL[:, i, j].mean()) / (NUL[:, i, j].std() or np.nan)
                   for i in range(len(ZN))])
    k = int(np.nanargmax(np.abs(np.nan_to_num(lg)))) if np.isfinite(lg).any() else MIDI
    kz = int(np.nanargmax(np.abs(np.nan_to_num(zz)))) if np.isfinite(zz).any() else MIDI
    rows.append((t, OBS[:, j], lg[k], OBS[k, j] - mid, zz[kz], ZN[kz]))
    g.write("%s\t%s\t%s\t%+.5f\t%+.2f\t%s\n"
            % (t, "\t".join("%.5f" % x for x in OBS[:, j]),
               ("%.4f" % lg[k]) if np.isfinite(lg[k]) else "",
               OBS[k, j] - mid, zz[kz] if np.isfinite(zz[kz]) else 0.0, ZN[kz]))
g.close()
print("\n■ 층 '%s' · 구역 사이 변화 (|귀무 z| 로 정렬)" % KEEP)
print("   %-14s %8s %8s %8s %8s | %8s %9s %s"
      % ("무리", ZN[0], ZN[1], ZN[2], ZN[3], "log2비", "귀무z", "가장 먼 구역"))
for r in sorted(rows, key=lambda r: -abs(r[4]) if np.isfinite(r[4]) else 0)[:14]:
    print("   %-14s %8.4f %8.4f %8.4f %8.4f | %8s %9.1f %s"
          % (r[0], r[1][0], r[1][1], r[1][2], r[1][3],
             ("%.2f" % r[2]) if np.isfinite(r[2]) else "-", r[4], r[5]))
print("\n   하플로타입 : " + " · ".join("%s %d" % (z, n) for z, n in zip(ZN, NHZ)))
print("저장 %s/VTZONE3_%s_zones.tsv" % (OUT, KEEP))
