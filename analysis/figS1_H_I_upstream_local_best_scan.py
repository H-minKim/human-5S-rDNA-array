# -*- coding: utf-8 -*-
# → ALIGNED_POS_RATE.tsv
#   단계 : S1 H·I 상류
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 배치된 자리에서 본 앵커 — 전량 스캔 + 역위 검출력
#   (A) 손대지 않은 정본 : 자리마다 '배치된 자리가 ±5 이웃 중 최고인 read 유닛 비율'(국소 최고율)
#       → 낮은 자리 · 낮은 자리가 이어진 구간 = 정본 이상 후보 (귀무 없이도 해석 가능하나 분포로 문턱을 정한다)
#   (B) 역위 검출력 : 구간 크기 3·6·12 × 위치를 옮기며 정본 유닛 내용을 뒤집고, 그 구간 국소 최고율이
#       문턱 아래로 떨어지는지 · 구간 밖은 그대로인지
import io, os, glob, collections, statistics

S = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(S, "inv_out"); os.makedirs(OUT, exist_ok=True)
MAPF = os.path.join(S, "pos_in", "CONSENSUS_MAP_G8_full.tsv")
V = os.path.join(os.environ["RDNA_V11"], "06_units_final")
TAB = chr(9); NL = chr(10); W = 5

# 정본 지도
CAN = {}; NUN = collections.Counter()
f = io.open(MAPF, encoding="utf8"); h = f.readline().rstrip(NL).split(TAB)
for L in f:
    r = dict(zip(h, L.rstrip(NL).split(TAB)))
    if r["varlist"] in ("", "-"): continue
    try: ap = int(r["array_pos"])
    except ValueError: continue
    if ap < 1: continue
    k = (r["sample"], r["slot"])
    CAN[(k[0], k[1], ap)] = set(r["varlist"].split(","))
    NUN[k] = max(NUN[k], ap)
cnt = collections.defaultdict(collections.Counter); tot = collections.Counter()
for (s, sl, ap), ts in CAN.items():
    tot[(s, sl)] += 1
    for t in ts: cnt[(s, sl)][t] += 1
FIX = {k: {t for t, v in cnt[k].items() if v >= 0.9 * tot[k]} for k in tot}
CANF = {k: (v - FIX[(k[0], k[1])]) for k, v in CAN.items()}

# 배치 read 유닛
def load(path, ci_slot, ci_pos, ci_var, has_read16):
    out = []
    for L in io.open(path, encoding="utf8"):
        a = L.rstrip(NL).split(TAB)
        if len(a) <= max(ci_slot, ci_pos, ci_var): continue
        rid = a[2] if has_read16 else a[1][:16]
        try: ap = int(a[ci_pos])
        except ValueError: continue
        if a[ci_var] in ("", "-"): continue
        out.append((a[0], a[ci_slot], rid, ap, set(a[ci_var].split(","))))
    return out
U = []
for p in sorted(glob.glob(V + "/AO_G8_UNITS_*.mid.tsv")): U += load(p, 3, 7, 14, True)
for p in sorted(glob.glob(V + "/JX_G8_UNITS_*.mid.tsv")): U += load(p, 2, 6, 12, False)
U = [(s, sl, rid, ap, ts - FIX.get((s, sl), set())) for (s, sl, rid, ap, ts) in U if (s, sl) in FIX]

def jac(a, b):
    if not a or not b: return 0.0
    i = len(a & b); return i / float(len(a) + len(b) - i)

# (A) 손대지 않은 정본
pos_hit = collections.defaultdict(list)
for s, sl, rid, ap, ts in U:
    c0 = CANF.get((s, sl, ap))
    if c0 is None: continue
    j0 = jac(ts, c0); best = j0; bap = ap
    for d in range(-W, W + 1):
        if d == 0: continue
        c = CANF.get((s, sl, ap + d))
        if c is None: continue
        v = jac(ts, c)
        if v > best: best, bap = v, ap + d
    pos_hit[(s, sl, ap)].append(1 if bap == ap else 0)
rate = {k: (sum(v) / float(len(v)), len(v)) for k, v in pos_hit.items()}
deep = {k: v for k, v in rate.items() if v[1] >= 5}
allr = [v[0] for v in deep.values()]
L = ["# 배치된 자리에서 본 앵커 — 전량 (2026-09-30)", "",
     "- 배치표 AO+JX 3 노드 전부 · 배치 유닛 %d · 정본 자리(깊이>=5) %d · 하플로 %d" % (
         len(U), len(deep), len(set((k[0], k[1]) for k in deep))), ""]
L.append("- 국소 최고율(배치 자리가 ±%d 중 최고인 read 유닛 비율) 분위(5/25/50/75/95) : %s" % (
    W, " / ".join("%.2f" % x for x in [sorted(allr)[min(len(allr) - 1, int(p * len(allr)))] for p in (0.05, 0.25, 0.5, 0.75, 0.95)])))
for th in (0.5, 0.3, 0.1):
    low = [k for k, v in deep.items() if v[0] <= th]
    runs = collections.Counter(); byh = collections.defaultdict(set)
    for s, sl, ap in low: byh[(s, sl)].add(ap)
    mx = 0; nrun = 0
    for k, aps in byh.items():
        for ap in sorted(aps):
            if ap - 1 not in aps:
                n = 1
                while ap + n in aps: n += 1
                if n >= 2: nrun += 1
                mx = max(mx, n)
    L.append("- 국소 최고율 <= %.1f 인 자리 : %d (%.2f%%) · 2 칸 이상 이어진 구간 %d · 최장 %d 칸" % (
        th, len(low), 100.0 * len(low) / len(deep), nrun, mx))
g = io.open(os.path.join(OUT, "ALIGNED_POS_RATE.tsv"), "w", encoding="utf8", newline=NL)
g.write(TAB.join(["sample", "slot", "array_pos", "local_best_rate", "n_units"]) + NL)
for k, v in sorted(deep.items()): g.write(TAB.join([k[0], k[1], str(k[2]), "%.3f" % v[0], str(v[1])]) + NL)
g.close()
low_rows = sorted([(v[0], v[1], k) for k, v in deep.items() if v[0] <= 0.3], key=lambda x: (x[0], -x[1]))
L.append("")
L.append("## 손대지 않은 정본에서 국소 최고율이 낮은 자리 (깊이>=5 · <=0.30) 상위")
L.append("| 표본 | 슬롯 | 자리 | 국소최고율 | 유닛 수 |")
L.append("|---|---|---|---|---|")
for r0, n0, k in low_rows[:20]:
    L.append("| %s | h%s | %d | %.2f | %d |" % (k[0], k[1], k[2], r0, n0))

# (B) 역위 검출력
byhap = collections.defaultdict(list)
for s, sl, rid, ap, ts in U: byhap[(s, sl)].append((ap, ts))
def block_rate(key, lo, hi, invert):
    s, sl = key; inb = []; outb = []
    can = CANF
    if invert:
        can = dict(CANF)
        for ap in range(lo, hi + 1):
            src = (s, sl, lo + hi - ap)
            if src in CANF: can[(s, sl, ap)] = CANF[src]
    for ap, ts in byhap[key]:
        c0 = can.get((s, sl, ap))
        if c0 is None: continue
        j0 = jac(ts, c0); best = j0; bap = ap
        for d in range(-W, W + 1):
            if d == 0: continue
            c = can.get((s, sl, ap + d))
            if c is None: continue
            v = jac(ts, c)
            if v > best: best, bap = v, ap + d
        (inb if lo <= ap <= hi else outb).append(1 if bap == ap else 0)
    return inb, outb
res = []
for key in sorted(byhap):
    n = NUN.get(key, 0)
    if n < 25: continue
    for blk in (3, 6, 12):
        for lo in range(6, n - blk - 5, 5):
            hi = lo + blk - 1
            i0, o0 = block_rate(key, lo, hi, False)
            if len(i0) < 5: continue
            i1, o1 = block_rate(key, lo, hi, True)
            r0 = sum(i0) / float(len(i0)); r1 = sum(i1) / float(len(i1))
            ro0 = sum(o0) / float(len(o0)) if o0 else 0; ro1 = sum(o1) / float(len(o1)) if o1 else 0
            res.append((key[0], key[1], blk, lo, hi, len(i0), r0, r1, ro0, ro1))
g = io.open(os.path.join(OUT, "ALIGNED_INV_POWER.tsv"), "w", encoding="utf8", newline=NL)
g.write(TAB.join(["sample", "slot", "block_units", "ap_lo", "ap_hi", "n_units_in_block",
                  "rate_in_real", "rate_in_inverted", "rate_out_real", "rate_out_inverted"]) + NL)
for r in res: g.write(TAB.join(str(x) if not isinstance(x, float) else "%.3f" % x for x in r) + NL)
g.close()
L.append("")
L.append("## 역위 검출력 — 규칙 : 구간 국소최고율 <= 0.5 (구간 밖은 그대로여야 함)")
L.append("| 구간 크기 | 시험 | 검출(역위판) | 거짓양성(진짜판) | 구간 안 최고율 중앙 진짜→역위 | 구간 밖 중앙 진짜→역위 |")
L.append("|---|---|---|---|---|---|")
for blk in (3, 6, 12):
    a = [r for r in res if r[2] == blk]
    if not a: continue
    L.append("| %d 유닛 | %d | **%.1f%%** | %.1f%% | %.2f → **%.2f** | %.2f → %.2f |" % (
        blk, len(a), 100.0 * sum(1 for r in a if r[7] <= 0.5) / len(a),
        100.0 * sum(1 for r in a if r[6] <= 0.5) / len(a),
        statistics.median([r[6] for r in a]), statistics.median([r[7] for r in a]),
        statistics.median([r[8] for r in a]), statistics.median([r[9] for r in a])))
io.open(os.path.join(OUT, "ALIGNED_SCAN_ALL.txt"), "w", encoding="utf8", newline=NL).write(NL.join(L) + NL)
print(NL.join(L))
