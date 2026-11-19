# -*- coding: utf-8 -*-
"""[Generated: `_make_cut508_variants.py` makes this file from zk1_m6a_by_state.py]
   Low and high methylation are split at `beta <= 0.508` (no U).
   m6A per position inside the unit, split into low and high methylation.

The split is by methylation state (low/high),
   not by region (distal/middle/proximal).
`M6A_SITE` is a table already pooled over reads, so no state can be attached to it.
`M6A_READUNIT` is used : it has read_id, and `calls` holds the ML score per position.
       columns : sample|data_source|cfg|kind|slot|array_pos|read_id|strand|mod_strand|n_a|n_m6a|frac|in_final|calls
       calls = list of "unit position:ML(0-255)"; m6A call = ML >= 128
       check : the total must agree with n_m6a of the table

The denominator is depth, as in the reference: the number of molecules covering that position.
The state is attached from the a13 GMM (per read x unit).
"""
import io
import os
import sys
import gzip
import glob
import collections
import numpy as np

# The data and output roots are environment variables with no default.
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

SRCD = os.path.join(_need("RDNA_FHMM_DIR"), "m6a_readunit")
#  This is written to: if it points at the published input folder, those tables are overwritten.
HF = _need("RDNA_FIBERSEQ_DIR")
OUT = os.path.join(HF, "01_latest_results", "REF2_20261011_m6a_in_unit_by_state_cut508.tsv")
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402


# ─────────────────────────────────────────────────────────────────
# Loader that splits low and high methylation at `beta <= 0.508`
#   - It returns the same key as a13: (sample, source, read_id, slot, array_pos) -> 'L' or 'H'.
#   - There is no undecided class: every unit is called one or the other.
#   - The robustness evidence is in the three-criteria robustness check, where all 18 combinations give the same sign.
CUT508 = 0.508


def _load_state_cut508():
    import gzip as _gz
    S = {}
    with _gz.open(os.path.join(L.D, "a13_unit_state_gmm.tsv.gz"), "rt",
                  encoding="utf8") as _f:
        _h = {k: i for i, k in enumerate(_f.readline().rstrip(chr(10)).split(chr(9)))}
        for _l in _f:
            _x = _l.rstrip(chr(10)).split(chr(9))
            try:
                _b = float(_x[_h["beta"]])
            except Exception:
                continue
            S[(_x[_h["sample"]], _x[_h["src"]], _x[_h["read_id"]],
               _x[_h["slot"]], _x[_h["array_pos"]])] = "L" if _b <= CUT508 else "H"
    return S
# ─────────────────────────────────────────────────────────────────

GMM = _load_state_cut508()
print("GMM %d" % len(GMM))
ML = 128

num = collections.Counter()      # (platform, state, unit position) -> the number of molecules carrying m6A
den = collections.Counter()      # (platform, state, unit position) -> the number of molecules covering that position
chk_m = chk_t = 0
nrow = 0
for p in sorted(glob.glob(os.path.join(SRCD, "*.tsv.gz"))):
    plat = "ONT" if "fsmeth" in os.path.basename(p) else "PacBio"
    src = "UCSC_ONT_FIBERSEQ" if plat == "ONT" else "PACBIO_FIBERSEQ"
    with gzip.open(p, "rt", encoding="utf8") as g:
        h = {k: i for i, k in enumerate(g.readline().rstrip(chr(10)).split(chr(9)))}
        for l in g:
            x = l.rstrip(chr(10)).split(chr(9))
            if len(x) != len(h):
                continue
            if "in_final" in h and x[h["in_final"]] != "Y":
                continue
            st = GMM.get((x[h["sample"]], src, x[h["read_id"]],
                          x[h["slot"]], x[h["array_pos"]]), "NA")
            if st not in ("L", "H"):
                continue
            nrow += 1
            c = x[h["calls"]]
            if not c:
                continue
            k_ = 0
            for it in c.split(","):
                if ":" not in it:
                    continue
                a_, b_ = it.split(":", 1)
                try:
                    up = int(a_); v = int(b_)
                except Exception:
                    continue
                den[(plat, st, up)] += 1
                if v >= ML:
                    num[(plat, st, up)] += 1
                    k_ += 1
            try:
                chk_m += abs(k_ - int(x[h["n_m6a"]])); chk_t += int(x[h["n_m6a"]])
            except Exception:
                pass
print("read x unit used %d" % nrow)
print("check: summed difference between m6A counted at ML>=%d and n_m6a of the table %d / %d (%.3f%%)"
      % (ML, chk_m, chk_t, 100.0 * chk_m / max(1, chk_t)))

rows = []
for (plat, st, up), d in sorted(den.items()):
    if d < 200:
        continue
    rows.append((plat, st, up, 100.0 * num[(plat, st, up)] / d, d))
with io.open(OUT, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# m6A by position within the unit: denominator is depth (number of molecules covering the position), by state" + chr(10))
    g.write("# m6A call = ML >= %d, check: difference from n_m6a of the table %.3f%%"
            % (ML, 100.0 * chk_m / max(1, chk_t)) + chr(10))
    g.write(chr(9).join(["플랫폼", "상태", "unit_pos", "m6a_분자비율_%", "깊이"]) + chr(10))
    for r in rows:
        g.write(chr(9).join([r[0], r[1], str(r[2]), "%.4f" % r[3], str(r[4])]) + chr(10))
print()
for plat in ("ONT", "PacBio"):
    for st in ("L", "H"):
        v = np.array([r[3] for r in rows if r[0] == plat and r[1] == st])
        if len(v):
            print("   %-7s %s: positions %4d, mean %6.2f%%"
                  % (plat, "저메틸" if st == "L" else "고메틸", len(v), v.mean()))
print()
print("written ->", os.path.basename(OUT))
