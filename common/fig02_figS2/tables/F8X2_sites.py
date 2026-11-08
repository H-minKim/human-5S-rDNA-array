# -*- coding: utf-8 -*-
r"""Neighbour exchange rate, choosing positions by a principle (instead of removing the 12 diagnostic positions by hand)

User definition
  positions used = only variant positions that satisfy both conditions at the same time
    (1) pooling the terminal-type units of all individuals and all positions, the position is polymorphic
       (among terminal-type units there are both units that carry it and units that do not)
    (2) pooling the non-terminal-type units in the same way, it is also polymorphic
  This automatically removes positions 'fixed in the terminal type' (= diagnostic variants) and positions 'fixed in the ordinary type'.

What is measured
  Neighbour pairs (gap 1 and 2) in the distal region within one haplotype are split into three classes,
  and the exchange rate (Jaccard = |A∩B| / |A∪B|) is looked at only at the chosen positions.
    terminal-terminal, non-terminal-non-terminal, terminal-non-terminal
  CI = haplotype bootstrap 1,000, with depth matching as well.
  Null = shuffling the terminal-type labels within haplotypes, 1,000 times.
Output: REPRO/F8X2_sites.tsv, F8X2_pairs.tsv, F8X2_depth.tsv, F8X2_null.tsv, F8X2_bytoken.tsv
"""
import io, os, sys, collections
import numpy as np

sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of the 12 diagnostic variants)

OUT = os.path.join((os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"]), "REPRO")
rng = np.random.default_rng(20260927)
NB, NPERM = 1000, 1000
GAPS = (1, 2)
DEPTH = int(os.environ.get("F8DEPTH", "15"))
MINC = int(os.environ.get("F8MINC", "5"))        # minimum number of carrier units and non-carrier units in each group
MARK = "D1096:1"

U, A_, M_ = G.load("primary", verbose=False)
HK = sorted(U)
# length floor: drop arrays that are too short and check sensitivity (exclude below F8MINN units)
MINN = int(os.environ.get("F8MINN", "0"))
_drop = [h for h in HK if len(U[h]) < MINN]
HK = [h for h in HK if len(U[h]) >= MINN]
if MINN:
    print("길이 바닥 N>=%d — 뺀 하플로타입 %d : %s"
          % (MINN, len(_drop), ", ".join("_".join(h) for h in _drop[:8])))


# ---- (1) split all units into terminal-type and non-terminal-type and count carriers per position ----
nT = nO = 0
cT, cO = collections.Counter(), collections.Counter()
for h in HK:
    for r in U[h]:
        isT = TD.is_term(r[4])
        if isT:
            nT += 1
            for t in r[4]:
                cT[t] += 1
        else:
            nO += 1
            for t in r[4]:
                cO[t] += 1
print("전체 유닛 — 말단형 %d · 비말단형 %d" % (nT, nO))

ALL = set(cT) | set(cO)
SITES, why = set(), collections.Counter()
for t in ALL:
    a, b = cT[t], cO[t]
    polyT = min(a, nT - a) >= MINC
    polyO = min(b, nO - b) >= MINC
    if polyT and polyO:
        SITES.add(t)
    else:
        why["말단형에서 다형성 아님" if not polyT else "비말단형에서 다형성 아님"] += 1
print("자리 — 전체 %d · 둘 다 다형성 %d (최소 %d)" % (len(ALL), len(SITES), MINC))
for k, v in why.most_common():
    print("   빠짐 : %-22s %d" % (k, v))

g = io.open(os.path.join(OUT, "F8X2_sites.tsv"), "w", encoding="utf8", newline="\n")
g.write("토큰\t말단형보유\t말단형전체\t말단형몫\t비말단형보유\t비말단형전체\t비말단형몫\t채택\n")
for t in sorted(ALL):
    g.write("%s\t%d\t%d\t%.5f\t%d\t%d\t%.5f\t%d\n"
            % (t, cT[t], nT, cT[t] / float(nT), cO[t], nO, cO[t] / float(nO), int(t in SITES)))
g.close()
# check that the 12 diagnostic positions were actually removed
DIAG = {MARK, "635CT", "1164GA", "740GA", "DB43:2", "DB16:17", "488GC", "DB46:10",
        "DB58:1", "R1752CG", "R405CG", "R347CG"}
print("   옛 진단 12 개 중 이 규칙으로도 빠진 것 : %d / %d  (남은 것 %s)"
      % (len(DIAG - SITES), len(DIAG), ",".join(sorted(DIAG & SITES)) or "없음"))

# ---- (2) neighbour pairs in the distal region within a haplotype ----
HAP = {}
for h in HK:
    rows = [(r[1], TD.is_term(r[4]), set(t for t in r[4] if t in SITES))
            for r in U[h] if r[1] <= DEPTH and r[1] <= r[0]]
    if len(rows) >= 3:
        HAP[h] = sorted(rows, key=lambda z: z[0])
CARR = [h for h in HAP if any(z[1] for z in HAP[h])]
print("원위 1-%d 쓸 수 있는 하플로타입 %d · 그중 말단형 보유 %d" % (DEPTH, len(HAP), len(CARR)))

KIND = {(True, True): "말단형-말단형", (False, False): "비말단형-비말단형"}


def pairs(h, flags=None):
    rows = HAP[h]
    isT = flags if flags is not None else [z[1] for z in rows]
    out = []
    for i in range(len(rows)):
        for gp in GAPS:
            j = i + gp
            if j >= len(rows):
                continue
            a, b = rows[i][2], rows[j][2]
            if not (a | b):
                continue
            k = KIND.get((isT[i], isT[j]), "말단형-비말단형")
            d = int(round((rows[i][0] + rows[j][0]) / 2.0))
            out.append((k, gp, d, len(a & b) / float(len(a | b)), len(a & b), len(a | b)))
    return out


OBS = collections.defaultdict(list)
CELL = collections.defaultdict(list)
for h in HAP:
    for k, gp, d, v, sh, un in pairs(h):
        OBS[(k, gp)].append((h, v, sh, un))
        CELL[(k, gp, d)].append(v)
        if gp <= 2:                      # gaps up to 2 are also accumulated as one group (position table 0)
            OBS[(k, 0)].append((h, v, sh, un))
            CELL[(k, 0, d)].append(v)
KS = ["말단형-말단형", "비말단형-비말단형", "말단형-비말단형"]


def hboot(rows):
    byh = collections.defaultdict(list)
    for h, v, _s, _u in rows:
        byh[h].append(v)
    hs = list(byh)
    out = []
    for _ in range(NB):
        s = rng.choice(len(hs), len(hs))
        vv = [x for i in s for x in byh[hs[i]]]
        out.append(np.mean(vv) if vv else np.nan)
    return np.nanpercentile(out, [2.5, 97.5])


g = io.open(os.path.join(OUT, "F8X2_pairs.tsv"), "w", encoding="utf8", newline="\n")
g.write("부류\t간격\t짝수\t하플로수\t자카드평균\tci_lo\tci_hi\t공유자리평균\t합집합평균\n")
print("\n=== 둘 다 다형성인 자리에서만 본 이웃 교환율 (원위 1-%d) ===" % DEPTH)
print("  %-20s %s" % ("부류", "간격≤2                 간격1                   간격2"))
for k in KS:
    line = "  %-20s" % k
    for gp in (0,) + GAPS:
        rows = OBS[(k, gp)]
        if len(rows) < 20:
            line += " %-23s" % ("(짝 %d · 모자람)" % len(rows))
            continue
        v = [x[1] for x in rows]
        lo, hi = hboot(rows)
        nh = len({x[0] for x in rows})
        g.write("%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n"
                % (k, gp, len(v), nh, np.mean(v), lo, hi,
                   np.mean([x[2] for x in rows]), np.mean([x[3] for x in rows])))
        line += " %.3f[%.3f-%.3f] n=%-5d" % (np.mean(v), lo, hi, len(v))
    print(line)
g.close()

# ---- depth matching ----
g = io.open(os.path.join(OUT, "F8X2_depth.tsv"), "w", encoding="utf8", newline="\n")
g.write("견줌\t간격\t맞춘칸\t말단형-비말단형평균\t상대평균\t차\t차_lo\t차_hi\n")
print("\n=== 깊이를 맞춘 차 ===")
for other in ("비말단형-비말단형", "말단형-말단형"):
    for gp in (0,) + GAPS:
        ds = [d for d in range(1, DEPTH + 1)
              if len(CELL[("말단형-비말단형", gp, d)]) >= 5 and len(CELL[(other, gp, d)]) >= 5]
        if len(ds) < 3:
            continue
        a = np.array([np.mean(CELL[("말단형-비말단형", gp, d)]) for d in ds])
        b = np.array([np.mean(CELL[(other, gp, d)]) for d in ds])
        dif = a - b
        bs = [np.mean(dif[rng.choice(len(ds), len(ds))]) for _ in range(NB)]
        lo, hi = np.percentile(bs, [2.5, 97.5])
        g.write("말단형-비말단형 vs %s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n"
                % (other, gp, len(ds), a.mean(), b.mean(), dif.mean(), lo, hi))
        print("  vs %-18s 간격%d · 칸 %2d · %.3f vs %.3f · 차 %+.3f [%+.3f,%+.3f]"
              % (other, gp, len(ds), a.mean(), b.mean(), dif.mean(), lo, hi))
g.close()

# ---- null ----
obs1 = np.mean([x[1] for x in OBS[("말단형-비말단형", 0)]])   # gap <=2
nul = []
for _ in range(NPERM):
    v = []
    for h in CARR:
        fl = [z[1] for z in HAP[h]]
        rng.shuffle(fl)
        for k, gp, d, j, sh, un in pairs(h, fl):
            if gp <= 2 and k == "말단형-비말단형":
                v.append(j)
    if v:
        nul.append(np.mean(v))
nul = np.array(nul)
p = (np.sum(nul <= obs1) + 1) / (len(nul) + 1.0)
print("\n=== 귀무 (딱지 섞기 %d 벌 · 간격≤2 말단형-비말단형) ===" % NPERM)
print("  관측 %.4f · 귀무 %.4f [%.4f-%.4f] · p %.4f"
      % (obs1, nul.mean(), np.percentile(nul, 2.5), np.percentile(nul, 97.5), p))
io.open(os.path.join(OUT, "F8X2_null.tsv"), "w", encoding="utf8", newline="\n").write(
    "stat\tvalue\n간격\t<=2\nobs_jaccard\t%.4f\nnull_mean\t%.4f\nnull_lo\t%.4f\nnull_hi\t%.4f\n"
    "p_one_sided\t%.4f\nnperm\t%d\nn_sites\t%d\n"
    % (obs1, nul.mean(), np.percentile(nul, 2.5), np.percentile(nul, 97.5), p, len(nul), len(SITES)))

# ---- per token ----
SHARED = collections.defaultdict(collections.Counter)
NPAIR = collections.Counter()
for h in HAP:
    rows = HAP[h]
    isT = [z[1] for z in rows]
    for i in range(len(rows)):
        for gp in (1, 2):                 # gaps up to 2
            j = i + gp
            if j >= len(rows):
                continue
            a, b = rows[i][2], rows[j][2]
            if not (a | b):
                continue
            k = KIND.get((isT[i], isT[j]), "말단형-비말단형")
            NPAIR[k] += 1
            for t in (a & b):
                SHARED[k][t] += 1
g = io.open(os.path.join(OUT, "F8X2_bytoken.tsv"), "w", encoding="utf8", newline="\n")
g.write("부류\t토큰\t공유짝수\t그부류짝수\t몫\n")
print("\n=== 간격 2 이내에서 함께 지닌 자리 상위 6 ===")
for k in KS:
    tot = NPAIR[k]
    if not tot:
        continue
    print("  %-20s 짝 %5d · %s" % (k, tot, " · ".join(
        "%s %.2f" % (t, c / float(tot)) for t, c in SHARED[k].most_common(6))))
    for t, c in SHARED[k].most_common(30):
        g.write("%s\t%s\t%d\t%d\t%.4f\n" % (k, t, c, tot, c / float(tot)))
g.close()
print("\n저장 REPRO/F8X2_*.tsv")

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('F8X2_pairs.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
