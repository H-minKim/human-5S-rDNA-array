# -*- coding: utf-8 -*-
r"""Recurrent mutation vs copying (2026-10-16: compute (1) and (2) and make the tables first)
   Set: same as P08_X7 (middle window 0.2-0.8, haplotypes with 23 or more units in the window, variants = substitutions and deletions, insertions excluded)
   (1) Co-sharing: unit pairs (i,j) that BOTH carry a rare variant A (cohort unit frequency < 1%) vs pairs (i,j') at the same distance where only ONE carries A.
      The fraction of other variants B, carried by at least one of the two units, that are present in BOTH units, by |pos_B - pos_A| bin.
      Copying = a large excess in near bins that shrinks with distance. Recurrence plus overall pair similarity = an excess independent of distance (flat).
      Interval estimate = haplotype cluster bootstrap, 500 rounds
   (2) Variant kinds: cohort-rare (< 1%) variants are split into those in only 1 unit of the array vs those shared by 2 or more units, and the composition of kinds is compared.
      If recurrent (hotspot), CpG C>T would concentrate on the shared side
   Read only. Output goes to the scratchpad"""
import os, sys, io, collections, pickle, time
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
SRC = __import__("os").environ.get("P04_SRC", __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
OUT = os.path.dirname(os.path.abspath(__file__))
BSTRAT = [("B 희귀 (<1%)", 0.0, 0.01), ("B 중간 (1-10%)", 0.01, 0.10), ("B 흔함 (>=10%)", 0.10, 1.01)]
SEED, NB = 20261016, 500
rng = np.random.default_rng(SEED)
RARE = 0.01
BINS = [(1, 10), (11, 100), (101, 250), (251, 500), (501, 1000), (1001, 2240)]   # 1-10 kept separate, a single event spanning several positions (deletion and the like)
CAN = G.CAN

U, _, _ = G.load("h1", verbose=False)
NCOPY = {("%s" % k[0], "%s" % k[1]): len(v) for k, v in U.items()}
with open(os.path.join(SRC, "OUT_P04", "Q0_alignment.pkl"), "rb") as f:
    DIFF = pickle.load(f)["diff"]
with open(os.path.join(SRC, "OUT_P04", "S3_matrices.pkl"), "rb") as f:
    S3 = pickle.load(f)
MATS, APS = S3["mats"], S3["aps"]
keys = [k for k in DIFF if k not in S3["bad"]]
A_OF = {}
for k in keys:
    sub, dele, ins, ul = DIFF[k]
    d = dict(sub)
    for p in dele:
        d[p] = "-"
    A_OF[k] = d
CF = collections.Counter()
for k in keys:
    for q, b in A_OF[k].items():
        CF[(q, b)] += 1
NT = len(keys)
# ---- check: coordinates are 1-based on the reference sequence, and substitution alleles differ from the reference
pos_all = [q for d in A_OF.values() for q in d]
subs = [(q, b) for d in A_OF.values() for q, b in d.items() if b != "-"]
same = sum(1 for q, b in subs if CAN[q - 1].upper() == b.upper())
print("참조 길이 %d · 자리 범위 %d-%d · 치환 %d 중 대립=참조 %d" % (len(CAN), min(pos_all), max(pos_all), len(subs), same))
assert len(CAN) == 2240 and 1 <= min(pos_all) and max(pos_all) <= 2240 and same == 0

PACK = []
for s in sorted(MATS):
    aps = APS[s]; N = NCOPY.get(s, 0)
    sel = [i for i, ap in enumerate(aps) if 0.2 * N < ap <= 0.8 * N]
    if len(sel) < 23:
        continue
    PACK.append((s, [A_OF[(s[0], s[1], aps[i])] for i in sel]))
print("하플로타입 %d · 유닛 %d" % (len(PACK), sum(len(x[1]) for x in PACK)))


def binof(d):
    for k, (a, b) in enumerate(BINS):
        if a <= d <= b:
            return k
    return None


# ---------------- (1) co-sharing ----------------
# per haplotype, (B shared count, B candidate count) for [bin][class (0 both, 1 one side)]
H = []
npair = [0, 0]
for s, UA in PACK:
    n = len(UA)
    acc = np.zeros((len(BSTRAT), len(BINS), 2, 2))
    carry = collections.defaultdict(list)
    for i in range(n):
        for q, b in UA[i].items():
            carry[(q, b)].append(i)
    for (qa, ba), car in carry.items():
        if len(car) < 2 or CF[(qa, ba)] / NT >= RARE:
            continue
        cs = set(car)
        for x in range(len(car)):
            for y in range(x + 1, len(car)):
                i, j = car[x], car[y]
                L = j - i
                # one-sided control: the unit on the opposite side of i at the same distance, or on the opposite side of j at the same distance, that lacks A
                ctrl = [(i, c) for c in (i - L,) if 0 <= c < n and c not in cs] + \
                       [(j, c) for c in (j + L,) if 0 <= c < n and c not in cs]
                for cls, pairs in ((0, [(i, j)]), (1, ctrl)):
                    for u, v in pairs:
                        npair[cls] += 1
                        Au, Av = UA[u], UA[v]
                        for q in set(Au) | set(Av):
                            if q == qa:
                                continue
                            kb = binof(abs(q - qa))
                            if kb is None:
                                continue
                            for bb in ([(q, Au[q])] if q in Au else []) + ([(q, Av[q])] if q in Av and (q not in Au or Av[q] != Au[q]) else []):
                                fr = CF[bb] / NT
                                si = [i_ for i_, (_, lo_, hi_) in enumerate(BSTRAT) if lo_ <= fr < hi_][0]
                                acc[si, kb, cls, 1] += 1
                                if q in Au and q in Av and Au[q] == Av[q]:
                                    acc[si, kb, cls, 0] += 1
    H.append(acc)
H = np.array(H)     # hap x strat x bin x cls x (shared, cand)
np.save(os.path.join(OUT, "COSHARE_BSTRAT_H.npy"), H)   # for the figure draft (haplotype-level bootstrap redone in the figure code)
# input table for the final figure (2026-10-16), (shared count, candidate count) per haplotype, B layer, distance bin and class
with io.open(os.path.join(OUT, "COSHARE_PERHAP.tsv"), "w", encoding="utf8", newline="\n") as gp:
    gp.write("# 희귀 변이 A(코호트 유닛 빈도 < 1%) 공동 공유 · H1 557 가운데 창(0.2-0.8 · 창 유닛 23 이상) · 생산 COSHARE_BSTRAT.py · "
             "cls 0 = A 를 함께 가진 쌍 · 1 = 같은 거리에서 한쪽만 가진 쌍 · bin 1-10·11-100 은 같은 사건 혼입이라 그림에서 뺀다\n")
    gp.write("sample\tslot\tB_stratum\tbin_bp\tcls\tn_shared\tn_cand\n")
    for hi_, (s_, _) in enumerate(PACK):
        for si_, (nm_, _, _) in enumerate(BSTRAT):
            for k_, (a_, b_) in enumerate(BINS):
                for c_ in (0, 1):
                    gp.write("%s\t%s\t%s\t%d-%d\t%d\t%d\t%d\n" % (s_[0], s_[1], nm_, a_, b_, c_,
                                                             H[hi_, si_, k_, c_, 0], H[hi_, si_, k_, c_, 1]))
print("\n① B 를 코호트 빈도로 나눔 · 쌍 : 함께 %d · 한쪽 %d" % tuple(npair))
g = io.open(os.path.join(OUT, "COSHARE_BSTRAT_out.tsv"), "w", encoding="utf8", newline="\n")
g.write("# B 코호트 빈도 층별 공동 공유 · 같은 집합 · cluster bootstrap %d · seed %d\n" % (NB, SEED))
g.write("B_stratum\tbin_bp\tn_together\tshare_together\tn_onesided\tshare_onesided\texcess\tlo\thi\n")
for si, (nm, _, _) in enumerate(BSTRAT):
    T = H[:, si].sum(0)
    r = T[:, :, 0] / T[:, :, 1]
    BS = []
    for b in range(NB):
        idx = rng.integers(0, len(H), len(H))
        Tb = H[idx, si].sum(0); rb = Tb[:, :, 0] / Tb[:, :, 1]; BS.append(rb[:, 0] - rb[:, 1])
    BS = np.array(BS)
    print("  ■ %s" % nm)
    for k, (a, b) in enumerate(BINS):
        lo, hi = np.percentile(BS[:, k], [2.5, 97.5])
        print("     %5d-%-5d  함께 %8d %.4f · 한쪽 %8d %.4f · 초과 %+.4f [%+.4f, %+.4f]" % (a, b, T[k, 0, 1], r[k, 0], T[k, 1, 1], r[k, 1], r[k, 0] - r[k, 1], lo, hi))
        g.write("%s\t%d-%d\t%d\t%.5f\t%d\t%.5f\t%.5f\t%.5f\t%.5f\n" % (nm, a, b, T[k, 0, 1], r[k, 0], T[k, 1, 1], r[k, 1], r[k, 0] - r[k, 1], lo, hi))
g.close()
# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('COSHARE_PERHAP.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
import sys; sys.exit(0)
# ---------------- (2) variant kinds ----------------
def vtype(q, b):
    if b == "-":
        return "결실"
    r = CAN[q - 1].upper(); b = b.upper()
    nxt = CAN[q].upper() if q < len(CAN) else ""
    prv = CAN[q - 2].upper() if q >= 2 else ""
    if (r == "C" and b == "T" and nxt == "G") or (r == "G" and b == "A" and prv == "C"):
        return "C>T CpG"
    if (r == "C" and b == "T") or (r == "G" and b == "A"):
        return "C>T 비CpG"
    if (r == "T" and b == "C") or (r == "A" and b == "G"):
        return "T>C"
    return "전환 아닌 치환(transversion)"


TY = ["C>T CpG", "C>T 비CpG", "T>C", "전환 아닌 치환(transversion)", "결실"]
HC = []    # per haplotype, [single (1 unit) / shared (2+)] x kind
for s, UA in PACK:
    c = np.zeros((2, len(TY)))
    carry = collections.Counter()
    for d in UA:
        for q, b in d.items():
            carry[(q, b)] += 1
    for (q, b), m in carry.items():
        if CF[(q, b)] / NT >= RARE:
            continue
        c[0 if m == 1 else 1, TY.index(vtype(q, b))] += 1
    HC.append(c)
HC = np.array(HC)
C = HC.sum(0)
print("\n② 코호트 희귀(< 1%%) 변이 · array 안 1 유닛만 %d ↔ 2 유닛 이상 공유 %d" % (C[0].sum(), C[1].sum()))
g.write("[B]\ttype\tsingle_n\tsingle_frac\tshared_n\tshared_frac\tdiff\tlo\thi\n")
BSC = []
for b in range(NB):
    idx = rng.integers(0, len(HC), len(HC))
    t = HC[idx].sum(0); BSC.append(t / t.sum(1, keepdims=True))
BSC = np.array(BSC)
for k, t in enumerate(TY):
    f0, f1 = C[0, k] / C[0].sum(), C[1, k] / C[1].sum()
    lo, hi = np.percentile(BSC[:, 1, k] - BSC[:, 0, k], [2.5, 97.5])
    print("   %-28s 홑 %6d (%.4f) · 공유 %6d (%.4f) · 차 %+.4f [%+.4f, %+.4f]" % (t, C[0, k], f0, C[1, k], f1, f1 - f0, lo, hi))
    g.write("\t%s\t%d\t%.5f\t%d\t%.5f\t%.5f\t%.5f\t%.5f\n" % (t, C[0, k], f0, C[1, k], f1, f1 - f0, lo, hi))
chi = stats.chi2_contingency(C)
print("   종류 구성 카이제곱 %.1f · 자유도 %d · p %.3g" % (chi[0], chi[2], chi[1]))
g.write("[B2]\tchi2\t%.3f\tdf\t%d\tp\t%.4g\n" % (chi[0], chi[2], chi[1]))
g.close()
print("\n산출 :", os.path.join(OUT, "COSHARE_TEST_out.tsv"))
