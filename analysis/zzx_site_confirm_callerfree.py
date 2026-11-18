# -*- coding: utf-8 -*-
"""Confirmation test: the two sites chosen by the HMM, tested again with a caller-free method.

Logic
  (1) main method = HMM footprint calls (FiberHMM). No site is fixed in advance; the whole grid is scanned.
     -> fold over the null: first the TFIIIB site, second the ICR site
  (2) only those two sites are carried over, and a method that uses no caller (flank-core m6A frequency ratio)
     tests whether the same sites exceed the null.

Why this is a sounder test than the rank
  The rank depends on how many other entries the grid holds (a finer width grid lowers the rank).
  In the confirmation test the window is fixed in advance, so multiplicity needs controlling only inside it :
    p = (permutations whose maximum inside the window under the null >= the observed maximum, + 1) / (permutations + 1)
  The maximum is taken over every position x width inside the window, so multiplicity inside the window is controlled.
  The window was chosen by another method (HMM), so the test itself is not circular.

Null : the (A, m6A) position table of each molecule is circularly shifted as a whole. Sequence composition and the per-molecule labelling rate are preserved.
"""
import io
import os
import sys
import gzip
import collections
import numpy as np

_SELF = os.path.basename(os.path.abspath(__file__))


def _need(k):
    # No defaults are set: if a variable is missing the script stops with a readable message
    v = os.environ.get(k, "")
    if not v:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s." % (k, _SELF)
            + chr(10)
            + "The deposited copy has no default; see README.md, section 'Environment variables'")
    return v


sys.path.insert(0, _need("RDNA_LIB_DIR"))
import a_lib as L

NL = chr(10)
TB = chr(9)
DD = _need("RDNA_DATA_ROOT")
OUTD = _need("RDNA_OUT_DIR")
MLTHR = 128
LO, HI = -150, 118
NP = HI - LO + 1
NB = 2000                                  # this is a confirmatory test, so the number of permutations is raised
GAP, FLANKW = 8, 17
WIDTHS = list(range(9, 72, 2))
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
SITES = [("TFIIIB 자리", -28, -10), ("ICR 자리", 57, 68)]
rel = np.arange(LO, HI + 1)
rng = np.random.default_rng(20261012)

UNIT, _ = L.load_units()
TSSU = {k: (v[2] - v[0] + 1) for k, v in UNIT.items()}
A, M = [], []
for node in ("w3", "iB", "iA"):
    p = os.path.join(L.R, "08_methylation_s5", node, L.METH["PACBIO_FIBERSEQ"],
                     "M6A_READUNIT_%s.tsv.gz" % L.HOST[node])
    if not os.path.exists(p):
        continue
    mol = collections.defaultdict(list)
    with gzip.open(p, "rt", encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if len(x) < 14 or x[1] != "PACBIO_FIBERSEQ" or x[12] != "Y" or not x[13]:
                continue
            smp, slot, ap, rid = x[0], x[4], x[5], x[6]
            if smp in L.PB_EXCL or (smp, slot) in SHIFT:
                continue
            mol[(smp, rid, slot, ap)].append(x[13])
    for k, blocks in mol.items():
        tu = TSSU.get((k[0], k[2], k[3]))
        if tu is None:
            continue
        a = np.zeros(NP, np.int32)
        m = np.zeros(NP, np.int32)
        for blk in blocks:
            for t in blk.split(","):
                if not t:
                    continue
                p0, v0 = t.split(":")
                r = int(p0) - tu
                if LO <= r <= HI:
                    a[r - LO] += 1
                    if int(v0) >= MLTHR:
                        m[r - LO] += 1
        if a.sum() >= 20:
            A.append(a)
            M.append(m)
A = np.array(A, np.int64)
M = np.array(M, np.int64)
nmol = len(A)
print("molecules %s, widths %d, permutations %d" % (format(nmol, ","), len(WIDTHS), NB))


def W(cs, c0, c1):
    o = np.full(NP, np.nan)
    i = np.arange(NP)
    x0, x1 = i + c0, i + c1 + 1
    ok = (x0 >= 0) & (x1 <= NP)
    o[ok] = cs[x1[ok]] - cs[x0[ok]]
    return o


def dprof(sa, sm, w):
    h = w // 2
    ca = np.concatenate([[0.0], np.cumsum(sa)])
    cm = np.concatenate([[0.0], np.cumsum(sm)])
    cA, cM = W(ca, -h, h), W(cm, -h, h)
    fA = W(ca, -h - GAP - FLANKW, -h - GAP - 1) + W(ca, h + GAP + 1, h + GAP + FLANKW)
    fM = W(cm, -h - GAP - FLANKW, -h - GAP - 1) + W(cm, h + GAP + 1, h + GAP + FLANKW)
    s = np.full(NP, np.nan)
    ok = np.isfinite(cA) & np.isfinite(fA) & (cA >= 200) & (fA >= 200)
    s[ok] = (fM[ok] / fA[ok]) / ((cM[ok] + 0.5) / (cA[ok] + 1.0))
    return s


MASK = {nm: ((rel >= a) & (rel <= b)) for nm, a, b in SITES}
sA, sM = A.sum(0).astype(float), M.sum(0).astype(float)

OBS = {nm: {} for nm, _, _ in SITES}
OBSMAX = {}
for nm, a, b in SITES:
    best = (-9, None, None)
    for w in WIDTHS:
        o = dprof(sA, sM, w)
        k = MASK[nm] & np.isfinite(o)
        if not k.any():
            continue
        j = int(np.nanargmax(np.where(k, o, np.nan)))
        OBS[nm][w] = (o[j], rel[j])
        if o[j] > best[0]:
            best = (o[j], rel[j], w)
    OBSMAX[nm] = best
    print("  %-10s observed max %.4f (position %+d, core %d bp)" % (nm, best[0], best[1], best[2]))

NMAX = {nm: np.empty(NB) for nm, _, _ in SITES}
NMAXW = {nm: {w: np.empty(NB) for w in WIDTHS} for nm, _, _ in SITES}
for t in range(NB):
    d = rng.integers(0, NP, nmol)
    ra = np.zeros(NP)
    rm = np.zeros(NP)
    for o in range(NP):
        sel = d == o
        if sel.any():
            ra += np.roll(A[sel].sum(0), o)
            rm += np.roll(M[sel].sum(0), o)
    pr = {w: dprof(ra, rm, w) for w in WIDTHS}
    for nm, _, _ in SITES:
        vs = []
        for w in WIDTHS:
            k = MASK[nm] & np.isfinite(pr[w])
            v = np.nanmax(np.where(k, pr[w], np.nan)) if k.any() else np.nan
            NMAXW[nm][w][t] = v
            if np.isfinite(v):
                vs.append(v)
        NMAX[nm][t] = max(vs) if vs else np.nan
    if (t + 1) % 500 == 0:
        print("   permutation %d / %d" % (t + 1, NB))

print()
print("=" * 94)
print("## confirmation test: the two positions chosen by the HMM, by a caller-free method (multiplicity within the window controlled)")
print("  %-10s %9s %8s %6s | %9s %9s %9s | %s"
      % ("자리", "관측최대", "자리", "폭", "귀무중앙", "귀무95%", "배수", "p (치환 %d)" % NB))
ROWS = []
for nm, a, b in SITES:
    ob, q, w = OBSMAX[nm]
    nd = NMAX[nm][np.isfinite(NMAX[nm])]
    med, p95 = float(np.median(nd)), float(np.percentile(nd, 95))
    pv = (int((nd >= ob).sum()) + 1) / float(len(nd) + 1)
    print("  %-10s %9.4f %+8d %6d | %9.4f %9.4f %9.3f | %.4g"
          % (nm, ob, q, w, med, p95, ob / med, pv))
    ROWS.append((nm, a, b, ob, q, w, med, p95, ob / med, pv, len(nd)))

print()
print("## each width separately (only position multiplicity within the window controlled)")
print("  %-10s %6s %9s %8s %9s %9s  %s" % ("자리", "폭", "관측", "자리", "귀무중앙", "배수", "p"))
PW = []
for nm, a, b in SITES:
    for w in WIDTHS:
        if w not in OBS[nm]:
            continue
        ob, q = OBS[nm][w]
        nd = NMAXW[nm][w]
        nd = nd[np.isfinite(nd)]
        if len(nd) < 50:
            continue
        pv = (int((nd >= ob).sum()) + 1) / float(len(nd) + 1)
        med = float(np.median(nd))
        PW.append((nm, w, ob, q, med, ob / med, pv))
        if w in (9, 15, 21, 25, 27, 33, 35, 37, 39, 41, 47, 57, 71):
            print("  %-10s %6d %9.4f %+8d %9.4f %9.3f  %.4g"
                  % (nm, w, ob, q, med, ob / med, pv))

out = os.path.join(OUTD, "TableS5_confirmation_test_no_caller.tsv")
with io.open(out, "w", encoding="utf8", newline=NL) as g:
    g.write("# Table S5, confirmation test: the two positions chosen by the HMM, retested by a caller-free method" + NL)
    g.write("# the windows come from the HMM side (TFIIIB -28..-10, ICR +57..+68), so this test is not circular." + NL)
    g.write("# p = (permutations where the null maximum within the window was at least the observed maximum + 1) / (permutations + 1)" + NL)
    g.write("#   the maximum is over all positions x widths within the window, so multiplicity within the window is controlled." + NL)
    g.write("# null : circular shift of the whole (A, m6A) position table per molecule, %d permutations" % NB + NL)
    g.write(TB.join(["자리", "창lo", "창hi", "관측최대", "그자리", "그폭",
                     "귀무중앙", "귀무95", "귀무대비배수", "p", "치환수"]) + NL)
    for r in ROWS:
        g.write(TB.join(["%s" % r[0], "%d" % r[1], "%d" % r[2], "%.4f" % r[3],
                         "%+d" % r[4], "%d" % r[5], "%.4f" % r[6], "%.4f" % r[7],
                         "%.4f" % r[8], "%.5g" % r[9], "%d" % r[10]]) + NL)
    g.write("#" + NL + "# each width separately (only position multiplicity within the window controlled)" + NL)
    g.write(TB.join(["자리", "core폭", "관측", "그자리", "귀무중앙", "배수", "p"]) + NL)
    for r in PW:
        g.write(TB.join(["%s" % r[0], "%d" % r[1], "%.4f" % r[2], "%+d" % r[3],
                         "%.4f" % r[4], "%.4f" % r[5], "%.5g" % r[6]]) + NL)
print()
print("written ->", out)
np.savez_compressed(os.path.join(DD, "site_confirm_callerfree.npz"),
                    nmax_T=NMAX[SITES[0][0]], nmax_I=NMAX[SITES[1][0]],
                    obs_T=OBSMAX[SITES[0][0]][0], obs_I=OBSMAX[SITES[1][0]][0],
                    wid=np.array(WIDTHS),
                    pw=np.array([[r[1], r[2], r[4], r[6]] for r in PW
                                 if r[0] == SITES[0][0]]),
                    pw_i=np.array([[r[1], r[2], r[4], r[6]] for r in PW
                                   if r[0] == SITES[1][0]]))
print("saved -> site_confirm_callerfree.npz")
