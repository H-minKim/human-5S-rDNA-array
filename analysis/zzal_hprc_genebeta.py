# -*- coding: utf-8 -*-
"""HPRC split by the stratum label used for the cell lines.

Two points
  (1) The 'whole pool' axis is not available in HPRC.
     Units that are not final in the placement (in_final=N) number 3,632 in FiberHMM,
     but only 3 of them span the window: they were filtered for being short.
     This axis is left out, and that count is written in the figure and the table.
  (2) The ruler follows the published version. `fhmm_peaks_bin10_both.npz` (the 10 bp row of table 3)
     removes only msp from the segments (nuc+tf). With tf alone the value is 49.85, with the published rule 52.77.
     The published rule is followed.

What is compared is therefore the stratum label alone
    whole-unit β 0.508 (main table)  against  gene body (last 119 bp, CpG>=5) β 0.5 (the cell-line way)

Written : hprc_genebeta.npz
  R   (label 2, site 2, stratum 2)  fold over the null at that site
  N   (label 2, stratum 2)          number of molecules
  CF  (2 x 2)                       confusion table of the two labels
  GC                                distribution of the gene-body CpG count
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

NL, TB = chr(10), chr(9)
DD = _need("RDNA_DATA_ROOT")
FD = os.path.join(DD, "fhmm")
SRC = "PACBIO_FIBERSEQ"
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
PER, HALF = 2240, 1120
NB, SM, MLTHR = 200, 11, 128
GENE_FROM, MIN_GENE_CPG = 2122, 5
CUT_UNIT, CUT_GENE = 0.508, 0.5
BINS = [("TFIIIB", 30, 39, -16), ("ICR", 40, 49, 63)]
NC = len(BINS)

# -- whole-unit beta: load_beta, exactly as in the published version
UB = {(s, r, sl, a): v for (s, sr, r, sl, a), v in
      L.load_beta(min_cpg=10).items() if sr == SRC}

# -- gene-region beta: from the calls, the last 119 bp only, with the two strand rows merged
raw = collections.defaultdict(lambda: [0, 0])
for node in ("w3", "iB", "iA"):
    p = os.path.join(L.R, "08_methylation_s5", node, L.METH[SRC],
                     "METH_READUNIT_%s.tsv.gz" % L.HOST[node])
    if not os.path.exists(p):
        continue
    with gzip.open(p, "rt", encoding="utf8") as f:
        ix = {k: i for i, k in enumerate(f.readline().rstrip(NL).split(TB))}
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[ix["data_source"]] != SRC or x[ix["in_final"]] != "Y":
                continue
            smp, sl = x[ix["sample"]], x[ix["slot"]]
            if smp in L.PB_EXCL or (smp, sl) in SHIFT:
                continue
            k = (smp, x[ix["read_id"]], sl, x[ix["array_pos"]])
            for t in (x[ix["calls"]] or "").split(","):
                if not t:
                    continue
                q, m = t.split(":")
                if int(q) >= GENE_FROM:
                    raw[k][0] += 1
                    raw[k][1] += (int(m) >= MLTHR)
GC = collections.Counter(v[0] for v in raw.values())
GB = {k: v[1] / float(v[0]) for k, v in raw.items() if v[0] >= MIN_GENE_CPG}
print("whole-unit β %s, gene-region β %s (CpG>=%d gate dropped %s)"
      % (format(len(UB), ","), format(len(GB), ","), MIN_GENE_CPG,
         format(len(raw) - len(GB), ",")))

# -- FiberHMM segments: the published rule, removing only the MSPs
SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio" or x[5] == "msp":
                continue
            SEG[(x[0], x[4], x[2], x[3])].append(
                ((int(x[6]) + int(x[7])) // 2, int(x[8])))

OBS = np.zeros((2, 2, NC, PER))
NUL = np.zeros((2, 2, NC, PER))
NMO = np.zeros((2, 2))
CF = np.zeros((2, 2))
rng = np.random.default_rng(20261012)   # the same seed as the published version (zzac), to match the ruler
for p in sorted(glob.glob(os.path.join(FD, "FHMM_READCOV_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio":
                continue
            if int(x[5]) > -HALF or int(x[6]) < HALF - 1:
                continue
            k = (x[0], x[4], x[2], x[3])
            labs = []
            if k in UB:
                labs.append((0, 0 if UB[k] < CUT_UNIT else 1))
            if k in GB:
                labs.append((1, 0 if GB[k] < CUT_GENE else 1))
            if not labs:
                continue
            if k in UB and k in GB:
                CF[0 if UB[k] < CUT_UNIT else 1, 0 if GB[k] < CUT_GENE else 1] += 1
            cen, cls = [], []
            for c, sz in SEG.get(k, []):
                if not (-HALF <= c < HALF):
                    continue
                for bi, (_n, a0, a1, _p) in enumerate(BINS):
                    if a0 <= sz <= a1:
                        cen.append(c + HALF)
                        cls.append(bi)
                        break
            cen = np.array(cen, dtype=int)
            cls = np.array(cls, dtype=int)
            d = rng.integers(0, PER, NB) if len(cen) else None
            for li, si in labs:
                NMO[li, si] += 1
                if not len(cen):
                    continue
                np.add.at(OBS[li, si], (cls, cen), 1)
                for t in range(NB):
                    np.add.at(NUL[li, si], (cls, (cen + d[t]) % PER), 1)
NUL /= float(NB)

K = np.ones(SM) / SM


def sm(v):
    return np.convolve(np.r_[v[-SM:], v, v[:SM]], K, "same")[SM:-SM]


R = np.zeros((2, NC, 2))
for li in (0, 1):
    for bi in range(NC):
        j = BINS[bi][3] + HALF
        for si in (0, 1):
            o, e = sm(OBS[li, si, bi]), sm(NUL[li, si, bi])
            R[li, bi, si] = o[j] / e[j] if e[j] > 0 else np.nan

np.savez_compressed(os.path.join(DD, "hprc_genebeta.npz"),
                    R=R, N=NMO, CF=CF, obs=OBS, nul=NUL,
                    gc=np.array(sorted(GC.items())),
                    sites=np.array([b[0] for b in BINS]),
                    pos=np.array([b[3] for b in BINS]))

LN = ["유닛 전체 β 0.508 (본표)", "유전자부 β 0.5 (세포주 방식)"]
print()
print("## fold over the null at that position, positions fixed (TFIIIB -16, ICR +63, 10 bp bins)")
print("%-26s %-8s %9s %9s %9s %14s" % ("stratum label", "자리", "저", "고", "low÷high", "분자 저/고"))
for li in (0, 1):
    for bi in range(NC):
        v = R[li, bi]
        print("%-26s %-8s %9.2f %9.2f %9.2f %14s"
              % (LN[li], BINS[bi][0], v[0], v[1],
                 v[0] / v[1] if v[1] > 0 else float("nan"),
                 "%s / %s" % (format(int(NMO[li, 0]), ","), format(int(NMO[li, 1]), ","))))
print()
print("## confusion table of the two labels (units with both)")
print("                      gene-region low    gene-region high")
for i, nm in enumerate(("유닛 저", "유닛 고")):
    print("   %-16s %11s %13s" % (nm, format(int(CF[i, 0]), ","), format(int(CF[i, 1]), ",")))
tot = CF.sum()
print("   agreement %.1f%% (n %s)" % (100 * (CF[0, 0] + CF[1, 1]) / tot, format(int(tot), ",")))
print()
print("## distribution of gene-region CpG count")
for k in sorted(GC):
    if GC[k] > 50:
        print("   CpG %2d : units %s" % (k, format(GC[k], ",")))
print()
print("saved -> hprc_genebeta.npz")
