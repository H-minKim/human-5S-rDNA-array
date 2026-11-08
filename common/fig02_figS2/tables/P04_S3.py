# -*- coding: utf-8 -*-
r"""S3 - alignment-based pairwise distance (the main scalar of Q2), built by projection onto CAN
   Section 2.3 consensus profile projection: each unit was aligned to CAN (2,240), so
   the distance between two units = number of CAN positions with different alleles + bp difference in inserted segments
   Construction
     A_i  = positions where unit i differs from CAN (substitution or deletion)
     diff = (A_i Δ A_j) ∪ {p ∈ A_i∩A_j : alleles differ}
     ins  = Hamming distance after padding the insertion sequences at the same anchor to the maximum length
   An n×n matrix is built for each slot and saved as npz (reused by Q2 and Q4)
   QC-failed units (69 units concentrated in one slot) are excluded and their number is written"""
import os, sys, io, re, time, collections, pickle
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "OUT_P04")
CAN = G.CAN
LU = len(CAN)
with open(os.path.join(OUT, "Q0_alignment.pkl"), "rb") as f:
    Q0 = pickle.load(f)
DIFF = Q0["diff"]
print("정렬 적재 유닛 %d" % len(DIFF))
# ---------- filter out QC-failed units ----------
SNV4 = re.compile(r"^R?(\d+)([ACGT])([ACGT])$")
BAD = set()
for f, Hx in G._rows(G.MAPF):
    if f[Hx["seg_type"]] != "UNIT":
        continue
    k = (f[0], f[Hx["slot"]], int(f[Hx["array_pos"]]))
    if k not in DIFF:
        continue
    vl = f[Hx["varlist"]]
    raw = [] if vl == "-" else vl.split(",")
    wsub, wd, wi = set(), 0, 0
    for t in raw:
        v = G.to_v4(t)
        if not v:
            continue
        m = SNV4.match(v)
        if m:
            wsub.add(v[1:] if v[0] == "R" else v)
            continue
        m = re.match(r"^D(\d+):(\d+)$", v)
        if m:
            wd += int(m.group(2))
            continue
        m = re.match(r"^I(\d+):(\d+)", v)
        if m:
            wi += int(m.group(2))
            continue
        m = re.match(r"^([DI])B(\d+):(\d+)$", v)
        if m:
            if m.group(1) == "D":
                wd += int(m.group(3))
            else:
                wi += int(m.group(3))
    sub, dele, ins, ul = DIFF[k]
    gsub = {"%d%s%s" % (p, CAN[p - 1], b) for p, b in sub.items()}
    gi = sum(len(s) for _, s in ins)
    if gsub != wsub or len(dele) != wd or gi != wi:
        BAD.add(k)
print("QC 실패 유닛 %d (%.5f) · 슬롯 분포 상위 %s"
      % (len(BAD), len(BAD) / len(DIFF),
         collections.Counter(("%s|%s" % (k[0], k[1])) for k in BAD).most_common(4)))

# ---------- S3 matrix per slot ----------
bySlot = collections.defaultdict(list)
for k in DIFF:
    if k in BAD:
        continue
    bySlot[(k[0], k[1])].append(k[2])
for s in bySlot:
    bySlot[s].sort()
print("슬롯 %d · 쓸 유닛 %d" % (len(bySlot), sum(len(v) for v in bySlot.values())))


def unit_rep(k):
    sub, dele, ins, ul = DIFF[k]
    # position -> allele code (0=CAN, 1-4=ACGT, 5=gap)
    A = {}
    for p, b in sub.items():
        A[p] = b
    for p in dele:
        A[p] = "-"
    ia = {}
    for p, s in ins:
        ia[p] = ia.get(p, "") + s
    return A, ia, ul


def pair_d(r1, r2):
    A1, I1, _ = r1
    A2, I2, _ = r2
    d = 0
    for p, b in A1.items():
        if A2.get(p, None) != b:
            d += 1
    for p in A2:
        if p not in A1:
            d += 1
    for p in set(I1) | set(I2):
        a, b = I1.get(p, ""), I2.get(p, "")
        if a != b:
            L = max(len(a), len(b))
            a = a.ljust(L, "\0")
            b = b.ljust(L, "\0")
            d += sum(1 for x, y in zip(a, b) if x != y)
    return d


t0 = time.time()
MATS = {}
for si, s in enumerate(sorted(bySlot)):
    aps = bySlot[s]
    reps = [unit_rep((s[0], s[1], ap)) for ap in aps]
    n = len(reps)
    M = np.zeros((n, n), np.int16)
    for i in range(n):
        for j in range(i + 1, n):
            M[i, j] = M[j, i] = pair_d(reps[i], reps[j])
    MATS[s] = M
    if si % 100 == 99:
        print("   ... 슬롯 %d/%d · %.0fs" % (si + 1, len(bySlot), time.time() - t0), flush=True)
print("■ S3 행렬 완료 %.0fs · 슬롯 %d" % (time.time() - t0, len(MATS)))
allv = np.concatenate([MATS[s][np.triu_indices(len(MATS[s]), 1)] for s in MATS])
print("   같은 array 쌍 %d · 편집거리 median %d · 평균 %.3f · 0 인 몫 %.6f · max %d"
      % (len(allv), int(np.median(allv)), allv.mean(), (allv == 0).mean(), allv.max()))
ulen = np.array([DIFF[k][3] for k in DIFF if k not in BAD])
np.savez_compressed(os.path.join(OUT, "S3_matrices.npz"),
                    slots=np.array(["%s|%s" % s for s in sorted(MATS)], dtype=object),
                    **{"%s|%s" % s: MATS[s] for s in MATS},
                    aps=np.array([np.array(bySlot[s]) for s in sorted(bySlot)], dtype=object))
with open(os.path.join(OUT, "S3_matrices.pkl"), "wb") as f:
    pickle.dump({"mats": MATS, "aps": {s: bySlot[s] for s in bySlot},
                 "bad": BAD, "ulen": {k: DIFF[k][3] for k in DIFF if k not in BAD}}, f, protocol=4)
g = io.open(os.path.join(OUT, "S3_summary.tsv"), "w", encoding="utf8", newline="\n")
g.write("# S3 정렬 기반 쌍별 거리 (CAN 투영) · parasail 정렬 · 2026-10-10\n")
g.write("# QC 실패 유닛 %d 제외 · 쓸 유닛 %d · 슬롯 %d\n"
        % (len(BAD), sum(len(v) for v in bySlot.values()), len(bySlot)))
g.write("item\tvalue\n")
for k_, v_ in (("n_unit_used", sum(len(v) for v in bySlot.values())), ("n_unit_bad", len(BAD)),
               ("n_slot", len(bySlot)), ("n_pair_within", len(allv)),
               ("median_edit", float(np.median(allv))), ("mean_edit", float(allv.mean())),
               ("frac_identical", float((allv == 0).mean())), ("max_edit", int(allv.max())),
               ("unit_len_median", float(np.median(ulen))), ("unit_len_min", int(ulen.min())),
               ("unit_len_max", int(ulen.max()))):
    g.write("%s\t%s\n" % (k_, v_))
g.close()
print("\n산출 :", os.path.join(OUT, "S3_matrices.pkl"), "·", os.path.join(OUT, "S3_summary.tsv"))
