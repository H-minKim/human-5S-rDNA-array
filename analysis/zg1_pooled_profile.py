# -*- coding: utf-8 -*-
"""Per-position profile over all reads pooled.

All reads are pooled; the proximal and distal ends are in absolute positions, the middle in relative position.
Why : one slot alone has only 5 to 20 molecules per position, and the di-nuc and patch tracks are too spiky.

Position frame
   array_pos 1 = genomic distal (DJ), array_pos N = proximal (centromeric)
   - distal 1..E      -> absolute number d = ap            (counted from the end)
   - proximal 1..E    -> absolute number p = N - ap + 1    (counted from the end)
   - middle           -> relative r = (ap - E) / (N - 2E)  in NB bins
   so the x axis = [E distal bins] [NB middle bins] [E proximal bins]

Measures (as in the reference, Dubocanin 2025)
   - di-nucleosome density = (bases covered by nucleosomes >210 bp) / (called bases)
   - accessible patch density = (bases covered by large MSPs >150 bp) / (called bases)
   - the count density (per kb) is also given
Pooling : per position bin all molecules are gathered, sum of numerators (bases) / sum of denominators (bases)  = pooled
   The spikiness goes away. The number of molecules (n) is written alongside.

Layout : as in the main table a23 (canonical coordinates, gaps are '.')
The two platforms are reported separately (absolute levels differ between platforms and are not compared)
"""
import io
import os
import sys
import gzip
import collections
import numpy as np

# The data and output roots are environment variables with no default (CONVENTIONS, "the deposited copy wins").
#  Paths to files inside the repository are resolved relative to __file__: a_lib, bamlite and p02_lib are deposited in this folder.
#  The table of variables is in README.md, section 'Environment variables'.
_HERE = os.path.dirname(os.path.abspath(__file__))


def _need(k):
    try:
        return os.environ[k]
    except KeyError:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s. "
            "The deposited copy has no default; see README.md, section 'Environment variables'."
            % (k, os.path.basename(__file__)))

A = _need("RDNA_ACCESS_DIR")
#  This is written to: if it points at the published input folder, those tables are overwritten.
HF = _need("RDNA_FIBERSEQ_DIR")
OUT = os.path.join(HF, "01_latest_results", "REF2_20261010c_merged_profile.tsv")
NPZ = os.path.join(_need("RDNA_FHMM_DIR"), "zg1_profile.npz")   # written here
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

UNIT, _ = L.load_units()
DINUC, LARGE = 210, 150
E = 20          # the number of units counted by absolute index from each end
NB = 10         # the number of relative bins in the interior

NMAX = collections.Counter()
for (smp, slot, ap) in UNIT:
    NMAX[(smp, slot)] = max(NMAX[(smp, slot)], int(ap))

# -- beta, per position bin
BETA = collections.defaultdict(list)


def slot_cell(ap, n):
    """position -> bin number (0..E-1 distal, E..E+NB-1 middle, E+NB..2E+NB-1 proximal); None if absent"""
    if n < 2 * E + 4:
        return None
    if ap <= E:
        return ap - 1
    if ap > n - E:
        return E + NB + (E - (n - ap) - 1)
    r = (ap - E - 1) / float(n - 2 * E)
    b = min(NB - 1, int(r * NB))
    return E + b


with gzip.open(os.path.join(A, "a13_unit_state_gmm.tsv.gz"), "rt", encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(chr(10)).split(chr(9))
        try:
            ap = int(x[4]); b = float(x[7])
        except Exception:
            continue
        n = NMAX.get((x[1], x[3]), 0)
        c = slot_cell(ap, n)
        if c is not None:
            BETA[(x[0], c)].append(b)
print("beta bins %d" % len(BETA))

NCELL = 2 * E + NB
din_bp = collections.Counter(); lmsp_bp = collections.Counter(); cov_bp = collections.Counter()
n_din = collections.Counter(); n_patch = collections.Counter(); n_mol = collections.Counter()

for smp, src, rid, slot, arr, lo, us in L.iter_reads(UNIT, 0, 0):
    if src == "PACBIO_FIBERSEQ" and smp in L.PB_EXCL:
        continue
    n = NMAX.get((smp, slot), 0)
    if n < 2 * E + 4:
        continue
    # which unit each base belongs to
    own = np.full(len(arr), -1, dtype=np.int32)
    for ap, _t in us:
        s, e, _ = UNIT[(smp, slot, ap)]
        own[s - 1 - lo:e - lo] = int(ap)
    cod = np.zeros(len(arr), dtype=np.int8)
    cod[arr == ord("N")] = 1
    cod[arr == ord("M")] = 2
    ch = np.flatnonzero(np.diff(cod)) + 1
    b = np.concatenate(([0], ch, [len(cod)]))
    segs = [(b[i], b[i + 1], cod[b[i]]) for i in range(len(b) - 1)]
    seen = set()
    for i in range(1, len(segs) - 1):
        a0, e0, c = segs[i]
        if c == 0 or segs[i - 1][2] == 0 or segs[i + 1][2] == 0:
            continue
        ln = e0 - a0
        big = (c == 1 and ln > DINUC) or (c == 2 and ln > LARGE)
        if not big:
            continue
        for p in range(a0, e0):
            ap = own[p]
            if ap < 0:
                continue
            cc = slot_cell(int(ap), n)
            if cc is None:
                continue
            if c == 1:
                din_bp[(src, cc)] += 1
            else:
                lmsp_bp[(src, cc)] += 1
        cm = slot_cell(int(own[(a0 + e0) // 2]), n) if own[(a0 + e0) // 2] >= 0 else None
        if cm is not None:
            (n_din if c == 1 else n_patch)[(src, cm)] += 1
    ok = (cod > 0) & (own >= 0)
    for p in np.flatnonzero(ok):
        cc = slot_cell(int(own[p]), n)
        if cc is not None:
            cov_bp[(src, cc)] += 1
            seen.add(cc)
    for cc in seen:
        n_mol[(src, cc)] += 1

print("bins %d (distal %d, interior %d, proximal %d)" % (NCELL, E, NB, E))


def lab(c):
    if c < E:
        return "원위%d" % (c + 1)
    if c < E + NB:
        return "가운데 %d/%d" % (c - E + 1, NB)
    return "근위%d" % (2 * E + NB - c)


rows = []
for src in ("UCSC_ONT_FIBERSEQ", "PACBIO_FIBERSEQ"):
    plat = "ONT" if src.startswith("UCSC") else "PacBio"
    for c in range(NCELL):
        cv = cov_bp[(src, c)]
        if cv < 1000:
            continue
        bs = BETA.get((src, c), [])
        rows.append((plat, lab(c), c,
                     100.0 * din_bp[(src, c)] / cv,
                     100.0 * lmsp_bp[(src, c)] / cv,
                     1000.0 * n_din[(src, c)] / cv,
                     1000.0 * n_patch[(src, c)] / cv,
                     np.mean(bs) if len(bs) >= 20 else float("nan"),
                     n_mol[(src, c)], cv, len(bs)))

with io.open(OUT, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# profile by position over all reads combined (pooled): sum of numerators / sum of denominators" + chr(10))
    g.write("# positions: distal 1..%d absolute, interior %d bins relative, proximal 1..%d absolute "
            "(array_pos 1 = distal/DJ, N = proximal/centromere)" % (E, NB, E) + chr(10))
    g.write("# di-nuc = fraction of bases covered by >210bp nucleosomes, patch = fraction of bases covered by >150bp large MSP"
            + chr(10))
    g.write(chr(9).join(["플랫폼", "자리", "칸번호", "dinuc_덮은%", "patch_덮은%",
                         "dinuc_개수_per_kb", "patch_개수_per_kb", "beta",
                         "분자수", "호출염기", "유닛수"]) + chr(10))
    for r in rows:
        g.write(chr(9).join([r[0], r[1], str(r[2]), "%.3f" % r[3], "%.3f" % r[4],
                             "%.4f" % r[5], "%.4f" % r[6],
                             ("%.4f" % r[7]) if r[7] == r[7] else "NA",
                             str(r[8]), str(r[9]), str(r[10])]) + chr(10))

print()
print("=== pooled profile (bases covered %)")
for plat in ("ONT", "PacBio"):
    sel = [r for r in rows if r[0] == plat]
    if not sel:
        continue
    print("   [%s]  %-14s %8s %8s %8s %8s" % (plat, "자리", "di-nuc%", "patch%", "beta", "분자"))
    for r in sel:
        if r[2] < 5 or (E <= r[2] < E + NB and (r[2] - E) % 3 == 0) or r[2] >= 2 * E + NB - 5:
            print("          %-14s %8.3f %8.3f %8s %8d"
                  % (r[1], r[3], r[4], ("%.3f" % r[7]) if r[7] == r[7] else "NA", r[8]))
    print()

np.savez_compressed(NPZ, E=E, NB=NB,
                    **{"%s|%s" % (r[0], k): np.array([x[i] for x in rows if x[0] == r[0]])
                       for r in rows[:1] for k, i in (("cell", 2), ("din", 3), ("patch", 4),
                                                      ("dinkb", 5), ("patchkb", 6), ("beta", 7),
                                                      ("nmol", 8))})
print("written ->", OUT)
