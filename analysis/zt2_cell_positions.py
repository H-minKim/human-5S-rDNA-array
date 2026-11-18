# -*- coding: utf-8 -*-
"""Count the unit positions per bin, to give a mean per position.

The profile reports the mean per unit position, not the count per bin.

The molecule count of a bin grows with the width of the bin:
   a middle bin merges 4.9 unit positions, and a molecule that touches only one of them still counts as 1,
   so it is not on the same ruler as an end bin (1 position).

Mean depth per position = (read x unit observations in the bin) / (positions merged in the bin)
   - read x unit observations = the number of rows of a13 (the unit count of the pooled profile)
   - positions = the number of distinct (individual, slot, array_pos) falling in the bin
The a01 tracks need not be read again: one read of a13 gives both.
The bin rule is that of zg1 (E=20, NB=10).
"""
import io
import os
import sys
import gzip
import collections

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

#  This is written to: if it points at the published input folder, those tables are overwritten.
A = _need("RDNA_ACCESS_DIR")
HF = _need("RDNA_FIBERSEQ_DIR")
OUT = os.path.join(HF, "01_latest_results", "REF2_20261012_by_cell_positions_and_depth.tsv")
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

E, NB = 20, 10
NCELL = 2 * E + NB
UNIT, _ = L.load_units()
NMAX = collections.Counter()
for (smp, slot, ap) in UNIT:
    NMAX[(smp, slot)] = max(NMAX[(smp, slot)], int(ap))


def slot_cell(ap, n):
    """The rule of zg1."""
    if n < 2 * E + 4:
        return None
    if ap <= E:
        return ap - 1
    if ap > n - E:
        return E + NB + (E - (n - ap) - 1)
    return E + min(NB - 1, int((ap - E - 1) / float(n - 2 * E) * NB))


def lab(c):
    if c < E:
        return "원위%d" % (c + 1)
    if c < E + NB:
        return "가운데 %d/%d" % (c - E + 1, NB)
    return "근위%d" % (2 * E + NB - c)


OBS = collections.Counter()                       # (source, cell) -> the number of read x unit observations
POS = collections.defaultdict(set)                # (source, cell) -> the set of (individual, slot, position)
nline = 0
with gzip.open(os.path.join(A, "a13_unit_state_gmm.tsv.gz"), "rt", encoding="utf8") as f:
    head = f.readline().rstrip(chr(10)).split(chr(9))
    print("a13 columns :", " | ".join(head))
    for l in f:
        x = l.rstrip(chr(10)).split(chr(9))
        nline += 1
        try:
            ap = int(x[4])
        except Exception:
            continue
        smp, src, slot = x[1], x[0], x[3]
        n = NMAX.get((smp, slot), 0)
        c = slot_cell(ap, n)
        if c is None:
            continue
        OBS[(src, c)] += 1
        POS[(src, c)].add((smp, slot, ap))
print("a13 rows %d, bins %d" % (nline, len(OBS)))

rows = []
for src in ("UCSC_ONT_FIBERSEQ", "PACBIO_FIBERSEQ"):
    plat = "ONT" if src.startswith("UCSC") else "PacBio"
    for c in range(NCELL):
        o = OBS.get((src, c), 0)
        k = len(POS.get((src, c), ()))
        if k == 0:
            continue
        rows.append((plat, lab(c), c, o, k, o / float(k)))

with io.open(OUT, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# per bin: number of unit positions and mean depth per position" + chr(10))
    g.write("# observations = number of read x unit rows in a13, positions = distinct (individual,slot,array_pos)"
            + chr(10))
    g.write("# depth per position = observations / positions: a measure not affected by bin width" + chr(10))
    g.write("# bin rule same as zg1 (both ends %d absolute, interior %d relative)" % (E, NB) + chr(10))
    g.write(chr(9).join(["플랫폼", "자리", "칸번호", "관측수", "자리수", "자리당깊이"]) + chr(10))
    for r in rows:
        g.write(chr(9).join([r[0], r[1], str(r[2]), str(r[3]), str(r[4]),
                             "%.2f" % r[5]]) + chr(10))

print()
print("=== view: end bins hold 1 position, interior bins hold several")
print("%-8s %-13s %10s %9s %11s" % ("플랫폼", "자리", "관측수", "자리수", "자리당깊이"))
for r in rows:
    if r[2] in (0, 1, 19, 20, 24, 29, 30, 48, 49):
        print("%-8s %-13s %10d %9d %11.2f" % (r[0], r[1], r[3], r[4], r[5]))
print()
print("written ->", os.path.basename(OUT))
