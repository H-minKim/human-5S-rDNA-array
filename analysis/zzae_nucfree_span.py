# -*- coding: utf-8 -*-
"""flank/core again, on the molecules without a nucleosome only.

Why
  The flank/core of the main table (MAIN f) pools every molecule.
  On one molecule a nucleosome at that position excludes a TF: the two are exclusive.
  Whether the ICR ratio is due to nucleosomes is therefore decided by measuring it again
  on the molecules without a nucleosome. (The caller-free version of SUPP (h)(i).)

Written : nucfree_flankcore_span.npz and a table on screen
  nA/nM  (site 2, stratum 2, nucleosome-free 2, rel 401)
  nmol   (site 2, stratum 2, nucleosome-free 2)
  ntf    (site 2, stratum 2, nucleosome-free 2)   <- occupancy inside that group (for the ceiling)

Definitions (those of zzab)
  'no nucleosome' = no `nuc` segment overlaps the core window
  'occupied'      = a `tf` call whose centre is inside the zone window and whose size is in range
"""
import io
import os
import sys
import gzip
import glob
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
FD = os.path.join(DD, "fhmm")
CUT = 0.508
MLTHR = 128
SRC = "PACBIO_FIBERSEQ"
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
LO, HI = -200, 200
NR = HI - LO + 1

# the positions: (centre, zone window, TF size range), the same values as the main table
SITE = [("TFIIIB", -11, (-28, -10), (30, 59)),
        ("ICR", 58, (57, 68), (20, 49))]
NS = len(SITE)

BETA = {(s, r, sl, a): v for (s, sr, r, sl, a), v in
        L.load_beta(min_cpg=10).items() if sr == SRC}
print("β n=%s" % format(len(BETA), ","))

NUC = collections.defaultdict(list)
TFS = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio":
                continue
            k = (x[0], x[4], x[2], x[3])
            if x[5] == "nuc":
                NUC[k].append((int(x[6]), int(x[7])))
            elif x[5] == "tf":
                TFS[k].append((int(x[6]), int(x[7]), int(x[8])))
print("molecules with nucleosome %s, molecules with TF %s"
      % (format(len(NUC), ","), format(len(TFS), ",")))

UNIT, _ = L.load_units()
TSSU = {k: (v[2] - v[0] + 1) for k, v in UNIT.items()}

nA = np.zeros((NS, 2, 2, NR))
nM = np.zeros((NS, 2, 2, NR))
nmol = np.zeros((NS, 2, 2))
ntf = np.zeros((NS, 2, 2))

for node in ("w3", "iB", "iA"):
    p = os.path.join(L.R, "08_methylation_s5", node, L.METH[SRC],
                     "M6A_READUNIT_%s.tsv.gz" % L.HOST[node])
    if not os.path.exists(p):
        continue
    MOL = collections.defaultdict(list)
    with gzip.open(p, "rt", encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if len(x) < 14 or x[1] != SRC or x[12] != "Y" or not x[13]:
                continue
            smp, slot, ap, rid = x[0], x[4], x[5], x[6]
            if smp in L.PB_EXCL or (smp, slot) in SHIFT:
                continue
            MOL[(smp, rid, slot, ap)].append(x[13])
    print("  %s : molecules %s" % (node, format(len(MOL), ",")))
    for k, blocks in MOL.items():
        b = BETA.get(k)
        tu = TSSU.get((k[0], k[2], k[3]))
        if b is None or tu is None:
            continue
        si = 0 if b < CUT else 1
        nu = NUC.get(k, ())
        tf = TFS.get(k, ())
        # nucleosome present or absent, and TF present or absent, per position
        flag = []
        for _, c0, zw, szr in SITE:
            # The symmetry condition: gating on the core alone leaves the flank ungated and tilts the ratio downwards,
            #   so the whole interval covering both the core and the flank, c0-35 to c0+27, is walked.
            a0, a1 = c0 - 35, c0 + 27
            cov = any(e >= a0 and s <= a1 for s, e in nu)
            hit = any(szr[0] <= sz <= szr[1] and zw[0] <= (s + e) // 2 <= zw[1]
                      for s, e, sz in tf)
            flag.append((0 if cov else 1, hit))     # 0 = a nucleosome lies in the interval; 1 = the whole interval is clear
        for j, (fr, hit) in enumerate(flag):
            nmol[j, si, fr] += 1
            ntf[j, si, fr] += hit
        # accumulating m6A
        for blk in blocks:
            for t in blk.split(","):
                if not t:
                    continue
                a_, b_ = t.split(":")
                r = int(a_) - tu
                if LO <= r <= HI:
                    m = int(b_) >= MLTHR
                    i = r - LO
                    for j, (fr, _h) in enumerate(flag):
                        nA[j, si, fr, i] += 1
                        nM[j, si, fr, i] += m

np.savez_compressed(os.path.join(DD, "nucfree_flankcore_span.npz"),
                    nA=nA, nM=nM, nmol=nmol, ntf=ntf, lo=LO, hi=HI, cut=CUT,
                    sites=np.array([s[0] for s in SITE]))

rel = np.arange(LO, HI + 1)


def ratio(j, si, fr):
    _, c0, _, _ = SITE[j]
    mc = (rel >= c0 - 7) & (rel <= c0 + 1)
    mf = ((rel >= c0 - 35) & (rel <= c0 - 27)) | ((rel >= c0 + 22) & (rel <= c0 + 27))
    cA, cM = nA[j, si, fr][mc].sum(), nM[j, si, fr][mc].sum()
    fA, fM = nA[j, si, fr][mf].sum(), nM[j, si, fr][mf].sum()
    cr = (cM + 0.5) / (cA + 0.5)
    fr_ = (fM + 0.5) / (fA + 0.5)
    return fr_ / cr, int(cA), int(fA)


print()
print("## flank/core: molecules with nucleosome vs without")
print("%-8s %-4s %-10s %9s %9s %9s %9s %9s"
      % ("자리", "층", "무리", "분자", "점유", "ceiling", "관측 비", "관측÷천장"))
for j, (nm, c0, _, _) in enumerate(SITE):
    for si, sl in ((0, "저"), (1, "고")):
        for fr, fn in ((0, "구간에 뉴클"), (1, "구간 깨끗")):
            n = nmol[j, si, fr]
            f = ntf[j, si, fr] / max(n, 1)
            ceil = 1.0 / max(1.0 - f, 1e-9)
            r, cA, fA = ratio(j, si, fr)
            print("%-8s %-4s %-10s %9s %8.4f %9.3f %9.3f %9.3f"
                  % (nm, sl, fn, format(int(n), ","), f, ceil, r, r / ceil))
print()
print("saved -> nucfree_flankcore_span.npz")
