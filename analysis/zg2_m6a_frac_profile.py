# -*- coding: utf-8 -*-
"""Measures (1) and (4): m6A with depth as the denominator.

The `M6A_SITE_*.tsv.gz` tables on the node hold this quantity.
   columns : sample | data_source | slot | array_pos | unit_pos | mod_strand | n | n_m6a | frac
        n     = number of molecules covering that base (depth)
        n_m6a = number of molecules with m6A at that base
        frac  = n_m6a / n  = "the fraction of molecules with m6A at that base"
   This is the reference definition : *"counting total occurrences across all molecules at a bp
      and dividing by the coverage at that bp"*

The `m6A concentration` has the number of A sites as its denominator. It is a different quantity.

Written
  (1) m6A fraction per position (inside the array): 15 distal absolute, 10 middle relative, 15 proximal absolute
  (4) profile per position inside the unit (unit_pos), per zone
Pooling : pooled (sum n_m6a / sum n), all molecules together
"""
import io
import os
import sys
import gzip
import glob
import collections
import numpy as np

# The data and output roots are environment variables with no default.
#  The table of variables is in README.md, section 'Environment variables'.

def _need(k):
    try:
        return os.environ[k]
    except KeyError:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s. "
            "The deposited copy has no default; see README.md, section 'Environment variables'."
            % (k, os.path.basename(__file__)))

W = _need("RDNA_FHMM_DIR")
SRCD = os.path.join(W, "m6a_site")
R = os.path.join(_need("RDNA_V12_CANON_DIR"), "02_canonical")
#  This is written to: if it points at the published input folder, those tables are overwritten.
HF = _need("RDNA_FIBERSEQ_DIR")
O1 = os.path.join(HF, "01_latest_results", "REF2_20261010d_m6a_depth_denominator_by_position.tsv")
O2 = os.path.join(HF, "01_latest_results", "REF2_20261010e_m6a_depth_denominator_in_unit.tsv")
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
E, NB = 20, 10      # the ends run to 20
NCELL = 2 * E + NB

NMAX = collections.Counter()
with io.open(os.path.join(R, "MAP_V11.tsv"), encoding="utf8") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip(chr(10)).split(chr(9)))}
    for l in f:
        x = l.rstrip(chr(10)).split(chr(9))
        if x[ix["seg_type"]] == "UNIT":
            k = (x[ix["sample"]], x[ix["slot"]])
            NMAX[k] = max(NMAX[k], int(x[ix["array_pos"]]))


def cell(ap, n):
    if n < 2 * E + 4:
        return None
    if ap <= E:
        return ap - 1
    if ap > n - E:
        return E + NB + (E - (n - ap) - 1)
    return E + min(NB - 1, int((ap - E - 1) / float(n - 2 * E) * NB))


def zone(c):
    return "원위" if c < E else ("가운데" if c < E + NB else "근위")


nm = collections.Counter(); nn = collections.Counter()           # ① (plat, cell)
um = collections.Counter(); un = collections.Counter()           # ④ (plat, zone, unit_pos)
rows_read = 0
for p in sorted(glob.glob(os.path.join(SRCD, "*.tsv.gz"))):
    plat = "ONT" if "fsmeth" in os.path.basename(p) else "PacBio"
    with gzip.open(p, "rt", encoding="utf8") as g:
        h = {k: i for i, k in enumerate(g.readline().rstrip(chr(10)).split(chr(9)))}
        for l in g:
            x = l.rstrip(chr(10)).split(chr(9))
            if len(x) != len(h):
                continue
            rows_read += 1
            try:
                ap = int(x[h["array_pos"]]); up = int(x[h["unit_pos"]])
                n = int(x[h["n"]]); k = int(x[h["n_m6a"]])
            except Exception:
                continue
            if n <= 0:
                continue
            N = NMAX.get((x[h["sample"]], x[h["slot"]]), 0)
            c = cell(ap, N)
            if c is None:
                continue
            nm[(plat, c)] += k; nn[(plat, c)] += n
            z = zone(c)
            um[(plat, z, up)] += k; un[(plat, z, up)] += n
print("rows read %d, position bins %d" % (rows_read, len(nn)))


def lab(c):
    if c < E:
        return "distal%d" % (c + 1)
    if c < E + NB:
        return "interior %d/%d" % (c - E + 1, NB)
    return "proximal%d" % (2 * E + NB - c)


with io.open(O1, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# metric (1): fraction of molecules with m6A at that base (denominator = depth), the reference definition" + chr(10))
    g.write("# pooled: sum(n_m6a) / sum(n), distal absolute %d, interior relative %d, proximal absolute %d"
            % (E, NB, E) + chr(10))
    g.write(chr(9).join(["platform", "position", "bin_index", "m6a_molecule_fraction_pct", "sum_n_m6a", "sum_n"]) + chr(10))
    for plat in ("ONT", "PacBio"):
        for c in range(NCELL):
            if nn[(plat, c)] < 1000:
                continue
            g.write(chr(9).join([plat, lab(c), str(c),
                                 "%.4f" % (100.0 * nm[(plat, c)] / nn[(plat, c)]),
                                 str(nm[(plat, c)]), str(nn[(plat, c)])]) + chr(10))

print()
print("=== metric (1): fraction of molecules with m6A (%), denominator is depth")
for plat in ("ONT", "PacBio"):
    v = [(lab(c), 100.0 * nm[(plat, c)] / nn[(plat, c)]) for c in range(NCELL)
         if nn[(plat, c)] >= 1000]
    if not v:
        continue
    print("   [%s]" % plat)
    for i, (L_, y) in enumerate(v):
        if i < 5 or (E <= i < E + NB and (i - E) % 3 == 0) or i >= NCELL - 5:
            print("      %-12s %6.3f" % (L_, y))
    d = np.mean([y for L_, y in v[:E]]); m = np.mean([y for L_, y in v[E:E + NB]])
    p_ = np.mean([y for L_, y in v[E + NB:]])
    print("      ---- distal mean %.3f, interior %.3f, proximal %.3f  (distal/interior %.2f-fold, proximal/interior %.2f-fold)"
          % (d, m, p_, d / m, p_ / m))

with io.open(O2, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# metric (4): fraction of molecules with m6A by position within the unit (unit_pos), per zone, denominator = depth" + chr(10))
    g.write(chr(9).join(["플랫폼", "구역", "unit_pos", "m6a_분자비율_%", "합_n"]) + chr(10))
    for plat in ("ONT", "PacBio"):
        for z in ("원위", "가운데", "근위"):
            for up in sorted({u for (p_, z_, u) in un if p_ == plat and z_ == z}):
                n = un[(plat, z, up)]
                if n < 200:
                    continue
                g.write(chr(9).join([plat, z, str(up),
                                     "%.4f" % (100.0 * um[(plat, z, up)] / n),
                                     str(n)]) + chr(10))
print()
print("written ->", os.path.basename(O1), "·", os.path.basename(O2))
