# -*- coding: utf-8 -*-
"""Per-position tally in one pass, local negative controls (sliding), and a verdict over both platforms.

Design
  (1) The four combinations are not independent replicates: ONT with unit β and ONT with gene-body β
     are the same reads (91,819 against 91,774), so only ONT and PacBio are independent.
  (2) A pass/fail at 5% alone would discard the fact that both platforms fall in the lowest 15%.
  (3) The elements A-box, IE and C-box lie inside the gene. The gene is transcribed and differs
     from the NTS (400 to 2,000 bp away from the gene) to begin with,
     so a control for the elements must also lie inside the gene.

What this script does
  - one pass stacks per-position (rel) nA and nM for each stratum (every later calculation uses this one array)
  - the control is built by sliding in each of three regions
      inside the gene (elements excluded) / upstream 90 / far NTS
  - the verdict combines the two platforms (both in the same tail, not a product of one-sided percentiles)
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

OUTD = _need("RDNA_OUT_DIR")
DD = _need("RDNA_DATA_ROOT")
NPZ = os.path.join(DD, "perbase_strat.npz")
NL = chr(10)
MLTHR = 128
LO, HI = -2100, 170
BGW = (-2000, -100)
MINCPG_GENE = 5
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
ST = ["저", "중", "고"]
AX = ["유닛 β", "유전자부 β"]
NW = HI - LO + 1
SRCS = [("UCSC_ONT_FIBERSEQ", "ONT"), ("PACBIO_FIBERSEQ", "PacBio")]
ELEM = [("TFIIIB", -30, -10), ("A-box", 49, 63), ("IE", 66, 71),
        ("C-box", 79, 96), ("ICR", 49, 96), ("유전자", 0, 118)]
# control regions: (name, start, end, and the window within which no element may overlap)
CTLZ = [("유전자안", 0, 118), ("상류90", -90, -1), ("먼 NTS", -2000, -400)]

if os.path.exists(NPZ) and os.environ.get("REBUILD", "") != "1":
    z = np.load(NPZ, allow_pickle=True)
    NAr, NMr, CUT, NU = z["nA"], z["nM"], z["cuts"], z["nunit"]
    print("rewriting the per-position summary ->", NPZ)
else:
    UNIT, _ = L.load_units()
    TSSU = {k: (v[2] - v[0] + 1) for k, v in UNIT.items()}
    NAr = np.zeros((2, 2, 3, NW))          # platform x stratification axis x stratum x position
    NMr = np.zeros((2, 2, 3, NW))
    NU = np.zeros((2, 2, 3))
    CUT = np.zeros((2, 2, 2))
    for pi, (SRC, nm_src) in enumerate(SRCS):
        UB = {(s, r, sl, a): v for (s, sr, r, sl, a), v in L.load_beta(min_cpg=10).items() if sr == SRC}
        GB = {}
        for node in ("w3", "iB", "iA"):
            p = os.path.join(L.R, "08_methylation_s5", node, L.METH[SRC],
                             "METH_READUNIT_%s.tsv.gz" % L.HOST[node])
            if not os.path.exists(p):
                continue
            with gzip.open(p, "rt", encoding="utf8") as f:
                f.readline()
                for l in f:
                    x = l.rstrip(NL).split(chr(9))
                    if len(x) < 13 or x[1] != SRC or x[11] != "Y" or not x[12]:
                        continue
                    tu = TSSU.get((x[0], x[4], x[5]))
                    if tu is None:
                        continue
                    v = np.array(x[12].replace(":", ",").split(","), dtype=np.int64).reshape(-1, 2)
                    rel = v[:, 0] - tu
                    m = (rel >= 0) & (rel <= 118)
                    n = int(m.sum())
                    if n >= MINCPG_GENE:
                        GB[(x[0], x[6], x[4], x[5])] = (v[m, 1] >= MLTHR).sum() / float(n)
        for ax, d in enumerate((UB, GB)):
            CUT[pi, ax] = np.percentile(np.array(list(d.values())), [33.333, 66.667])
        print("%-7s unit β n %d cuts %.4f/%.4f, gene-region β n %d cuts %.4f/%.4f"
              % (nm_src, len(UB), CUT[pi, 0, 0], CUT[pi, 0, 1], len(GB), CUT[pi, 1, 0], CUT[pi, 1, 1]))
        used = 0
        for node in ("w3", "iB", "iA"):
            p = os.path.join(L.R, "08_methylation_s5", node, L.METH[SRC],
                             "M6A_READUNIT_%s.tsv.gz" % L.HOST[node])
            if not os.path.exists(p):
                continue
            with gzip.open(p, "rt", encoding="utf8") as f:
                f.readline()
                for l in f:
                    x = l.rstrip(NL).split(chr(9))
                    if len(x) < 14 or x[1] != SRC or x[12] != "Y" or not x[13]:
                        continue
                    smp, slot, ap, rid = x[0], x[4], x[5], x[6]
                    if SRC == "PACBIO_FIBERSEQ" and (smp in L.PB_EXCL or (smp, slot) in SHIFT):
                        continue
                    tu = TSSU.get((smp, slot, ap))
                    if tu is None:
                        continue
                    bs = (UB.get((smp, rid, slot, ap)), GB.get((smp, rid, slot, ap)))
                    if bs[0] is None and bs[1] is None:
                        continue
                    v = np.array(x[13].replace(":", ",").split(","), dtype=np.int64).reshape(-1, 2)
                    rel = v[:, 0] - tu
                    m = (rel >= LO) & (rel <= HI)
                    if not m.any():
                        continue
                    idx = rel[m] - LO
                    hit = (v[m, 1] >= MLTHR)
                    for ax in (0, 1):
                        b = bs[ax]
                        if b is None:
                            continue
                        si = 0 if b < CUT[pi, ax, 0] else (1 if b < CUT[pi, ax, 1] else 2)
                        np.add.at(NAr[pi, ax, si], idx, 1)
                        if hit.any():
                            np.add.at(NMr[pi, ax, si], idx[hit], 1)
                        NU[pi, ax, si] += 1
                    used += 1
            print("   read %s, cumulative %d" % (node, used))
    np.savez_compressed(NPZ, nA=NAr, nM=NMr, cuts=CUT, nunit=NU, lo=LO, hi=HI)
    print("saved ->", NPZ)


def inter(pi, ax, a, b):
    """interaction of window [a,b] = [rate(high)/rate(low)] / [background rate(high)/rate(low)], and the minimum observed A"""
    i0, i1 = a - LO, b - LO + 1
    A = NAr[pi, ax, :, i0:i1].sum(1)
    M = NMr[pi, ax, :, i0:i1].sum(1)
    j0, j1 = BGW[0] - LO, BGW[1] - LO + 1
    bA = NAr[pi, ax, :, j0:j1].sum(1)
    bM = NMr[pi, ax, :, j0:j1].sum(1)
    if A[0] < 300 or A[2] < 300:
        return float("nan"), min(A[0], A[2])
    r = (M / A)
    br = (bM / bA)
    return (r[2] / r[0]) / (br[2] / br[0]), min(A[0], A[2])


def nulls(pi, ax, size, zname, z0, z1, excl):
    """slide a window of size over region [z0,z1] and give the distribution of the interaction; windows overlapping excl are left out"""
    out = []
    for s in range(z0, z1 - size + 2):
        e = s + size - 1
        if any(not (e < ea or s > eb) for ea, eb in excl):
            continue
        v, n = inter(pi, ax, s, e)
        if v == v:
            out.append(v)
    return np.array(out)


EXCL = [(-30, -10), (49, 96)]        # positions carrying an element are excluded from the control
print()
print("=== element interaction and in-place negative control (sliding, two platforms, stratification axes 2)")
rows = []
for ax in (0, 1):
    print()
    print("---- stratification axis: %s" % AX[ax])
    print("  %-8s %5s | %-26s | %-26s" % ("요소", "크기", "ONT  값 / 대조중앙 / 백분위", "PacBio  값 / 대조중앙 / 백분위"))
    for nm, a, b in ELEM:
        size = b - a + 1
        zname, z0, z1 = ("상류90", -90, -1) if nm == "TFIIIB" else ("유전자안", 0, 118)
        if size > (z1 - z0 + 1) - 20:
            zname, z0, z1 = "먼 NTS", -2000, -400
        cells = []
        for pi, (_, src) in enumerate(SRCS):
            v, n = inter(pi, ax, a, b)
            cn = nulls(pi, ax, size, zname, z0, z1, EXCL)
            if v != v or len(cn) < 8:
                cells.append("  값없음 (대조 %d)" % len(cn)); continue
            pc = 100.0 * (cn < v).mean()
            cells.append("%6.3f / %6.3f / %5.1f%% (n%d)" % (v, np.median(cn), pc, len(cn)))
            rows.append((AX[ax], src, nm, size, zname, v, float(np.median(cn)),
                         float(np.percentile(cn, 5)), float(np.percentile(cn, 95)), pc, len(cn), int(n)))
        print("  %-8s %5d | %-26s | %-26s" % (nm, size, cells[0], cells[1]))
        print("           %5s   control zone = %s" % ("", zname))

print()
print("=== verdict with the two platforms combined (the independent replicates are the two, ONT and PacBio)")
print("  %-8s %-10s %10s %10s %10s" % ("요소", "층화축", "ONT 백분위", "PB percentile", "판정"))
for ax in (0, 1):
    for nm, a, b in ELEM:
        pcs = [r[9] for r in rows if r[0] == AX[ax] and r[2] == nm]
        if len(pcs) < 2:
            continue
        lo_both = all(p <= 15 for p in pcs)
        hi_both = all(p >= 85 for p in pcs)
        vd = "양쪽 하위 15%" if lo_both else ("양쪽 상위 15%" if hi_both else "어긋남/안")
        print("  %-8s %-10s %9.1f%% %9.1f%% %10s" % (nm, AX[ax], pcs[0], pcs[1], vd))

out = os.path.join(OUTD, "TF_20261010_in_place_negative_control_sliding.tsv")
with io.open(out, "w", encoding="utf8", newline=NL) as g:
    g.write(chr(9).join(["층화축", "플랫폼", "요소", "창크기", "대조구역", "교호작용",
                         "대조중앙", "대조5%", "대조95%", "백분위", "대조창수", "A관측최소"]) + NL)
    for r in rows:
        g.write(chr(9).join([r[0], r[1], r[2], "%d" % r[3], r[4]]
                            + ["%.4f" % r[i] for i in range(5, 10)]
                            + ["%d" % r[10], "%d" % r[11]]) + NL)
print()
print("written ->", out)
