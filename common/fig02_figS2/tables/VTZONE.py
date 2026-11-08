# -*- coding: utf-8 -*-
r"""Which variant type changes the most across the zones defined by the S24 boundaries

Zones (S24 boundaries: proximal 3|4, 17|18, distal 10|11)
   proximal 1-3, proximal 4-17, middle, distal 1-10
   The distal side has only one boundary, so there is no 'distal 2'. There are 4 zones.
   Overlapping units (short arrays that are shallow on both the proximal and distal side) go to the shallower side.

What is measured
   - two-stage mean per zone (mean per unit within haplotype -> mean across haplotypes)
   - change = max |log2(zone / middle)|, absolute change = max|zone - middle|
   - null = units shuffled within haplotype, 200 sets -> z
Output: REPRO/VTZONE_zones.tsv
"""
import io, os, sys
import numpy as np

CODE = __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir))
sys.path.insert(0, CODE)
import v11_io as G

OUT = os.path.join((os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"]), "REPRO")
os.makedirs(OUT, exist_ok=True)
CAN = G.CAN
rng = np.random.default_rng(20260913)
NN = 200

U, _, _ = G.load("primary", verbose=False)
HK = sorted(U)


def iscpg(p):
    c = CAN[p - 1]
    return (c == "C" and p < 2240 and CAN[p] == "G") or (c == "G" and p >= 2 and CAN[p - 2] == "C")


def vtype(x):
    if x.startswith("DB"):
        return "반복블록 결실"
    if x.startswith("IB"):
        return "반복블록 삽입"
    if x[0] == "D":
        p, n = x[1:].split(":")
        return "1bp 결실" if int(n) == 1 else "긴 결실"
    if x[0] == "I":
        q = x[1:].split(":")
        return "1bp 삽입" if int(q[1]) == 1 else "긴 삽입"
    p = int(x[1:-2]) if x[0] == "R" else int(x[:-2])
    s = x[-2] + ">" + x[-1]
    if s in ("C>T", "G>A"):
        return "탈아미노 CpG" if iscpg(p) else "탈아미노 비CpG"
    if s in ("C>A", "G>T"):
        return "산화 C>A·G>T"
    if s in ("T>C", "A>G"):
        return "그 밖 전이 T>C·A>G"
    return "그 밖 전환"


TY = ["1bp 결실", "반복블록 결실", "산화 C>A·G>T", "탈아미노 비CpG", "탈아미노 CpG",
      "그 밖 전이 T>C·A>G", "그 밖 전환", "긴 결실", "1bp 삽입", "긴 삽입", "반복블록 삽입"]
TI = {t: i for i, t in enumerate(TY)}
NT = len(TY)
ZN = ["근위 1-3", "근위 4-17", "가운데", "원위 1-10"]

VEC, AP, AD, HI = [], [], [], []
for hi, k in enumerate(HK):
    for row in U[k]:
        v = np.zeros(NT)
        for x in row[4]:
            v[TI[vtype(x)]] += 1
        VEC.append(v); AP.append(row[0]); AD.append(row[1]); HI.append(hi)
VEC = np.asarray(VEC); AP = np.asarray(AP); AD = np.asarray(AD); HI = np.asarray(HI)
NH = len(HK)


def zone_of(ap, ad):
    """If overlapping, assign to the shallower side"""
    z = np.full(len(ap), 2, int)                       # middle
    pz = np.where(ap <= 3, 0, np.where(ap <= 17, 1, -1))
    dz = np.where(ad <= 10, 3, -1)
    both = (pz >= 0) & (dz >= 0)
    z = np.where(pz >= 0, pz, z)
    z = np.where(dz >= 0, 3, z)
    z = np.where(both, np.where(ap <= ad, pz, 3), z)   # the shallower side
    return z, int(both.sum())


Z0, NOV = zone_of(AP, AD)
print("유닛 %d · 하플로타입 %d · 두 구역에 다 걸리는 유닛 %d (%.2f%%) -> 더 얕은 쪽"
      % (len(VEC), NH, NOV, 100 * NOV / len(VEC)))
for i, z in enumerate(ZN):
    print("   %-10s 유닛 %6d" % (z, int((Z0 == i).sum())))


def zmeans(zz):
    out = np.full((len(ZN), NT), np.nan)
    nh = np.zeros(len(ZN), int)
    for i in range(len(ZN)):
        acc = []
        for hi in range(NH):
            s = (zz == i) & (HI == hi)
            n = int(s.sum())
            if n:
                acc.append(VEC[s].sum(0) / n)
        out[i] = np.mean(np.asarray(acc), 0)
        nh[i] = len(acc)
    return out, nh


OBS, NHZ = zmeans(Z0)
MIDI = ZN.index("가운데")

ORD = {hi: np.where(HI == hi)[0] for hi in range(NH)}
NUL = np.empty((NN, len(ZN), NT))
for b in range(NN):
    ap2 = AP.copy(); ad2 = AD.copy()
    for hi, ix in ORD.items():
        q = rng.permutation(len(ix))
        ap2[ix] = AP[ix][q]; ad2[ix] = AD[ix][q]
    NUL[b] = zmeans(zone_of(ap2, ad2)[0])[0]

g = io.open(os.path.join(OUT, "VTZONE_zones.tsv"), "w", encoding="utf8", newline="\n")
g.write("유형\t" + "\t".join(ZN) + "\t최대_log2비\t최대_절대차\t최대_귀무z\t가장_먼_구역\n")
rank = []
for t in TY:
    j = TI[t]
    mid = OBS[MIDI, j]
    lg = np.array([np.log2(OBS[i, j] / mid) if (mid > 0 and OBS[i, j] > 0) else
                   (-np.inf if mid > 0 else np.nan) for i in range(len(ZN))])
    zz = np.array([(OBS[i, j] - NUL[:, i, j].mean()) / (NUL[:, i, j].std() or np.nan)
                   for i in range(len(ZN))])
    k = int(np.nanargmax(np.abs(np.where(np.isfinite(lg), lg, 0))))
    kz = int(np.nanargmax(np.abs(np.where(np.isfinite(zz), zz, 0))))
    rank.append((t, OBS[:, j], lg[k], OBS[k, j] - mid, zz[kz], ZN[kz]))
    g.write("%s\t%s\t%.4f\t%+.5f\t%+.2f\t%s\n"
            % (t, "\t".join("%.5f" % x for x in OBS[:, j]),
               lg[k] if np.isfinite(lg[k]) else float("nan"),
               OBS[k, j] - mid, zz[kz], ZN[kz]))
g.close()

print("\n■ 구역 사이에서 가장 많이 변하는 유형 (|귀무 z| 로 정렬)")
if 0:
    print("x"
      % (["유형"] + ZN + ["최대 log2비", "최대 절대차", "어느 구역"])[0],
      end="")
print("")
print("   %-20s %8s %8s %8s %8s | %9s %9s %s"
      % ("유형", ZN[0], ZN[1], ZN[2], ZN[3], "log2비", "귀무z", "가장 먼 구역"))
for r in sorted(rank, key=lambda r: -abs(r[4]) if np.isfinite(r[4]) else 0):
    print("   %-20s %8.4f %8.4f %8.4f %8.4f | %9.2f %9.1f %s"
          % (r[0], r[1][0], r[1][1], r[1][2], r[1][3], r[2], r[4], r[5]))
print("\n   하플로타입 수 : " + " · ".join("%s %d" % (z, n) for z, n in zip(ZN, NHZ)))
print("저장 %s/VTZONE_zones.tsv" % OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('VTZONE_zones.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
