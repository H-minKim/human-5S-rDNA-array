# -*- coding: utf-8 -*-
"""The calculation of zz7 with the null not thinned.

Here the null is not thinned (SUB=1), so nothing is multiplied.
The output is written as fhmm_peaks_v2.npz.
"""
import io
import os
import sys
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

DD = _need("RDNA_DATA_ROOT")
FD = os.path.join(DD, "fhmm")
NL = chr(10)
TB = chr(9)
PER = 2240                      # the unit period
HALF = PER // 2                 # 1120
NB = 200                        # null replicates
SUB = 1                         # no thinning: the scale factor is 1
SM = 11                         # smoothing by position, in bp
CLS = [(10, 19, "10~19"), (20, 29, "20~29"), (30, 39, "30~39"), (40, 49, "40~49"),
       (50, 59, "50~59"), (60, 89, "60~89"), (90, 140, "90~140"), (141, 200, "141~200")]
NC = len(CLS)
PLAT = {"pacbio": 1, "nanopore": 0}
SRCN = {"pacbio": "PACBIO_FIBERSEQ", "nanopore": "UCSC_ONT_FIBERSEQ"}

old = np.load(os.path.join(DD, "perbase_strat.npz"), allow_pickle=True)
CUTB = old["cuts"][:, 0, :]
BETA = {}
for plat, src in SRCN.items():
    for (s, sr, r, sl, a), v in L.load_beta(min_cpg=10).items():
        if sr == src:
            BETA[(plat, s, r, sl, a)] = v

SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[5] == "msp":
                continue
            SEG[(x[1], x[0], x[4], x[2], x[3])].append((int(x[6]), int(x[7]), int(x[8])))
print("read x unit with intervals %s" % format(len(SEG), ","))

OBS = np.zeros((2, 3, NC, PER))         # platform x stratum x size group x position
NUL = np.zeros((2, 3, NC, PER))
NMOL = np.zeros((2, 3))
NSUB = np.zeros((2, 3))
rng = np.random.default_rng(20261011)
nskip = 0

for p in sorted(glob.glob(os.path.join(FD, "FHMM_READCOV_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            plat, smp, sl, ap, rid = x[1], x[0], x[2], x[3], x[4]
            lo, hi = int(x[5]), int(x[6])
            if lo > -HALF or hi < HALF - 1:      # only molecules covering the window completely
                nskip += 1
                continue
            b = BETA.get((plat, smp, rid, sl, ap))
            if b is None:
                continue
            pi = PLAT[plat]
            si = 0 if b < CUTB[pi, 0] else (1 if b < CUTB[pi, 1] else 2)
            fp = SEG.get((plat, smp, rid, sl, ap), [])
            cen, cls = [], []
            for (s0, s1, sz) in fp:
                c = (s0 + s1) // 2
                if not (-HALF <= c < HALF):
                    continue
                for k, (a0, a1, _) in enumerate(CLS):
                    if a0 <= sz <= a1:
                        cen.append(c + HALF)
                        cls.append(k)
                        break
            NMOL[pi, si] += 1
            if not cen:
                continue
            cen = np.array(cen)
            cls = np.array(cls)
            np.add.at(OBS[pi, si], (cls, cen), 1)
            if int(NMOL[pi, si]) % SUB:
                continue
            NSUB[pi, si] += 1
            d = rng.integers(0, PER, NB)
            for t in range(NB):
                np.add.at(NUL[pi, si], (cls, (cen + d[t]) % PER), 1)

NUL /= float(NB)
np.savez_compressed(os.path.join(DD, "fhmm_peaks_v2.npz"),
                    obs=OBS, nul=NUL, nmol=NMOL, nsub=NSUB,
                    per=PER, half=HALF, cls=np.array([c[2] for c in CLS]))
print("molecules dropped for not covering the whole window %s" % format(nskip, ","))
print("molecule count  low / mid / high")
for nm, pi in (("PacBio", 1), ("ONT", 0)):
    print("  %-8s %s   (null samples %s)"
          % (nm, " / ".join("%d" % NMOL[pi, s] for s in range(3)),
             " / ".join("%d" % NSUB[pi, s] for s in range(3))))

# -- the peak table: observed over null, after smoothing
K = np.ones(SM) / SM


def sm(v):
    return np.convolve(np.r_[v[-SM:], v, v[:SM]], K, "same")[SM:-SM]


print()
print("=" * 86)
print("## peaks: observed ÷ positional null (smoothing %d bp, scaled to molecule count)" % SM)
for nm, pi in (("PacBio", 1), ("ONT", 0)):
    print()
    print("  %s" % nm)
    print("    %-9s %10s %10s %10s %10s %10s"
          % ("크기", "저 자리", "low fold", "중 배수", "high fold", "low observed/1000"))
    for k, (a0, a1, lab) in enumerate(CLS):
        o0 = sm(OBS[1 if nm == "PacBio" else 0, 0, k])
        n0 = sm(NUL[1 if nm == "PacBio" else 0, 0, k])     # scale 1
        r0 = np.where(n0 > 0, o0 / np.maximum(n0, 1e-9), np.nan)
        j = int(np.nanargmax(np.where(np.arange(PER) < PER, r0, np.nan)))
        pos = j - HALF
        row = [r0[j]]
        for si in (1, 2):
            oo = sm(OBS[pi, si, k])
            nn = sm(NUL[pi, si, k])                            # scale 1
            row.append(oo[j] / nn[j] if nn[j] > 0 else np.nan)
        print("    %-9s %+10d %10.2f %10.2f %10.2f %10.2f"
              % (lab, pos, row[0], row[1], row[2], 1000.0 * o0[j] / max(NMOL[pi, 0], 1)))
print()
print("saved -> fhmm_peaks_v2.npz")
