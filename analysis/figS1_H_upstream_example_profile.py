# -*- coding: utf-8 -*-
# → EXAMPLE_PROFILE.tsv
#   단계 : S1 H 상류
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 그림 C 패널용 실측 예시 : 한 하플로에서 '손대지 않은 정본' vs '가운데 6 유닛 역위' 의 자리별 국소최고율
# Output : EXAMPLE_PROFILE.tsv (the figure reads this table only)
import io, os, glob, collections
K = os.environ["RDNA_HANDOFF_AD6468DD"]
OUT = os.path.join(K, "result_inversion_positive_control_20260930")
S = os.path.dirname(os.path.abspath(__file__))
MAPF = os.path.join(S, "pos_in", "CONSENSUS_MAP_G8_full.tsv")
V = os.path.join(os.environ["RDNA_V11"], "06_units_final")
TAB = chr(9); NL = chr(10); W = 5
SMP, SLOT, BLK = "NA19118", "1", 6
def rows(p):
    f = io.open(p, encoding="utf8"); h = f.readline().rstrip(NL).split(TAB)
    for L in f:
        if L.strip(): yield dict(zip(h, L.rstrip(NL).split(TAB)))
CAN = {}; N = 0
for r in rows(MAPF):
    if r["sample"] != SMP or r["slot"] != SLOT or r["varlist"] in ("", "-"): continue
    try: ap = int(r["array_pos"])
    except ValueError: continue
    if ap < 1: continue
    CAN[ap] = set(r["varlist"].split(",")); N = max(N, ap)
c = collections.Counter()
for ts in CAN.values():
    for t in ts: c[t] += 1
FIX = {t for t, v in c.items() if v >= 0.9 * len(CAN)}
CANF = {ap: (ts - FIX) for ap, ts in CAN.items()}
U = []
for p in sorted(glob.glob(V + "/AO_G8_UNITS_*.mid.tsv")):
    for L in io.open(p, encoding="utf8"):
        a = L.rstrip(NL).split(TAB)
        if len(a) > 14 and a[0] == SMP and a[3] == SLOT and a[14] not in ("", "-"):
            try: U.append((int(a[7]), set(a[14].split(",")) - FIX))
            except ValueError: pass
for p in sorted(glob.glob(V + "/JX_G8_UNITS_*.mid.tsv")):
    for L in io.open(p, encoding="utf8"):
        a = L.rstrip(NL).split(TAB)
        if len(a) > 12 and a[0] == SMP and a[2] == SLOT and a[12] not in ("", "-"):
            try: U.append((int(a[6]), set(a[12].split(",")) - FIX))
            except ValueError: pass
def jac(x, y):
    if not x or not y: return 0.0
    i = len(x & y); return i / float(len(x) + len(y) - i)
lo = (N - BLK) // 2 + 1; hi = lo + BLK - 1
def profile(inv):
    can = dict(CANF)
    if inv:
        for ap in range(lo, hi + 1):
            if (lo + hi - ap) in CANF: can[ap] = CANF[lo + hi - ap]
    hit = collections.defaultdict(list)
    for ap, ts in U:
        if ap not in can: continue
        best = jac(ts, can[ap]); bap = ap
        for d in range(-W, W + 1):
            if d == 0 or (ap + d) not in can: continue
            v = jac(ts, can[ap + d])
            if v > best: best, bap = v, ap + d
        hit[ap].append(1 if bap == ap else 0)
    return {ap: (sum(v) / float(len(v)), len(v)) for ap, v in hit.items()}
p0 = profile(False); p1 = profile(True)
g = io.open(os.path.join(OUT, "EXAMPLE_PROFILE.tsv"), "w", encoding="utf8", newline=NL)
g.write(TAB.join(["sample", "slot", "array_pos", "n_units", "rate_real", "rate_inverted", "in_block"]) + NL)
for ap in sorted(set(list(p0.keys()) + list(p1.keys()))):
    a0 = p0.get(ap); a1 = p1.get(ap)
    if not a0: continue
    g.write(TAB.join([SMP, SLOT, str(ap), str(a0[1]), "%.3f" % a0[0], "%.3f" % (a1[0] if a1 else 0.0),
                      "Y" if lo <= ap <= hi else "N"]) + NL)
g.close()
inb0 = [p0[ap][0] for ap in range(lo, hi + 1) if ap in p0]; inb1 = [p1[ap][0] for ap in range(lo, hi + 1) if ap in p1]
print("%s h%s N=%d · 심은 구간 %d-%d" % (SMP, SLOT, N, lo, hi))
print("  구간 안 국소최고율 : 진짜 %.2f → 역위 %.2f" % (sum(inb0) / len(inb0), sum(inb1) / len(inb1)))
out0 = [p0[ap][0] for ap in p0 if not (lo <= ap <= hi)]; out1 = [p1[ap][0] for ap in p1 if not (lo <= ap <= hi)]
print("  구간 밖              : 진짜 %.2f → 역위 %.2f" % (sum(out0) / len(out0), sum(out1) / len(out1)))
