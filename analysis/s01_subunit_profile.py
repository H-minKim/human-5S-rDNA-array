# -*- coding: utf-8 -*-
"""S01: per-position (sub-unit) profile of m6A per A site inside the unit, separately for low and high methylation.

Why
  A bar per zone is a mean and does not show the shape inside the zone. The raw calls (unit_pos:ML) are counted as they are
  to give the profile in unit coordinates (real bp), separately for L and H.

Input (read only)
  <RDNA_V12_CANON_DIR>/08_methylation_s6/<node>/s6/{ont_meth,pb_meth}/M6A_READUNIT_*.tsv.gz
      calls column = "unit_pos:ML,..." (unit_pos = 1-based inside the unit, A sites of both strands), in_final=Y only
  unit coordinates and TSS : a_lib.load_units()  (MAP_V11 and the canonical sequence; TSS = the last GTCTACGGC in the unit)
  unit state : a13_unit_state_gmm.tsv.gz of the s6 version (L P>=0.9, H P<=0.1, U left out)

Method
  tx = unit_pos - tss_u  (tss_u = 1-based TSS position inside the unit) -> gene tx 0..118; AluY is about tx -2120..-1800
  per BIN bp bin, sum(m6a)/sum(A) per individual -> mean over individuals ± 95% CI (only bins with A>=30 within the individual)
Output
  s01_subunit_profile.tsv     (src, state, tx_bin, n_sample, mean, ci95, n_a)
  s01_element_bounds.tsv      (median zone boundaries, for the x-axis ticks of the figure)
usage : python s01_subunit_profile.py [BIN]      default BIN=20
"""
import io, os, sys, gzip, collections
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

sys.path.insert(0, _HERE)
import a_lib as L

BIN = int(sys.argv[1]) if len(sys.argv) > 1 else 20
#  This is written to: if it points at the published input folder, those tables are overwritten.
OUT = _need("RDNA_RESULTMAP_DIR")
S6 = os.path.join(L.R, "08_methylation_s6")
MDIR = {"UCSC_ONT_FIBERSEQ": "ont_meth", "PACBIO_FIBERSEQ": "pb_meth"}
# The PacBio cohort is the 48 jasmine samples, excluding two called with fibertools 0.4.2 (HG02630, hg002) and two from jasmine 2.2.0 with no m6A (HG02624, HG03672).
EXCL = {"PACBIO_FIBERSEQ": {"HG02630", "hg002", "HG02624", "HG03672"}, "UCSC_ONT_FIBERSEQ": set()}
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)

print("loading canonical units and TSS")
UNIT, _SEQ = L.load_units()
TSSU = {}                      # (sample, slot, array position) -> the TSS within the unit, 1-based, and the unit length
for k, (us, ue, tss) in UNIT.items():
    TSSU[k] = (tss - us + 1, ue - us + 1)
print("  units %d, median TSS within unit %d, median unit length %d"
      % (len(TSSU), int(np.median([v[0] for v in TSSU.values()])), int(np.median([v[1] for v in TSSU.values()]))))

print("loading s6 states")
ST = {}
with gzip.open(os.path.join(_need("RDNA_S6_DIR"), "a13_unit_state_gmm.tsv.gz"),
               "rt", encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip("\n").split("\t")
        ST[(x[1], x[0], x[2], x[3], x[4])] = x[10]     # (sample, src, read, slot, apos) -> L/H/U
print("  unit states %d" % len(ST))

acc = collections.defaultdict(lambda: [0, 0])          # (src, state, sample, bin) -> [A, m6A]
chk = {"row": 0, "no_unit": 0, "no_state": 0, "U": 0, "pos_gt_len": 0, "used": 0}
for src, sub in MDIR.items():
    for node in ("w3", "iB", "iA"):
        p = os.path.join(S6, node, "s6", sub, "M6A_READUNIT_%s.tsv.gz" % L.HOST[node])
        if not os.path.exists(p):
            print("  none", p); continue
        with gzip.open(p, "rt", encoding="utf8") as f:
            ix = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
            for l in f:
                x = l.rstrip("\n").split("\t")
                chk["row"] += 1
                if x[ix["data_source"]] != src or x[ix["in_final"]] != "Y":
                    continue
                smp = x[ix["sample"]]
                if smp in EXCL[src]:
                    continue
                key = (smp, x[ix["slot"]], x[ix["array_pos"]])
                if key not in TSSU:
                    chk["no_unit"] += 1; continue
                st = ST.get((smp, src, x[ix["read_id"]], x[ix["slot"]], x[ix["array_pos"]]))
                if st is None:
                    chk["no_state"] += 1; continue
                if st not in ("L", "H"):
                    chk["U"] += 1; continue
                tss_u, ulen = TSSU[key]
                calls = x[ix["calls"]]
                if not calls:
                    continue
                for c in calls.split(","):
                    q = c.split(":")
                    if len(q) != 2:
                        continue
                    pos = int(q[0])
                    if pos > ulen:
                        chk["pos_gt_len"] += 1
                    b = ((pos - tss_u) // BIN) * BIN
                    v = acc[(src, st, smp, b)]
                    v[0] += 1
                    if int(q[1]) >= 128:
                        v[1] += 1
                chk["used"] += 1
        print("  %s %s finished (cumulative read×unit %d)" % (src, node, chk["used"]))
print("check:", chk)

# the mean over individuals with a 95% confidence interval
per = collections.defaultdict(dict)
for (src, st, smp, b), (na, nm) in acc.items():
    if na >= 30:
        per[(src, st, b)][smp] = nm / float(na)
tot = collections.defaultdict(lambda: [0, 0])
for (src, st, smp, b), (na, nm) in acc.items():
    t = tot[(src, st, b)]
    t[0] += na; t[1] += nm

with io.open(os.path.join(OUT, "s01_subunit_profile.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("src\tstate\ttx_bin\tbin_bp\tn_sample\tmean\tci95\tpooled_frac\tn_a\n")
    for (src, st, b) in sorted(per, key=lambda k: (k[0], k[1], k[2])):
        v = np.array(list(per[(src, st, b)].values()), dtype=float)
        if len(v) < 5:
            continue
        ci = 1.96 * v.std(ddof=1) / np.sqrt(len(v))
        na, nm = tot[(src, st, b)]
        g.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\n"
                % (src, st, b, BIN, len(v), v.mean(), ci, nm / float(na), na))

# the segment boundaries, in within-unit coordinates converted to the transcription frame, as medians
tss_med = int(np.median([v[0] for v in TSSU.values()]))
bounds = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA", 1179, tss_med - 571),
          ("NTS2", tss_med - 570, tss_med - 91), ("upstream90", tss_med - 90, tss_med - 1),
          ("gene", tss_med, tss_med + 118)]
with io.open(os.path.join(OUT, "s01_element_bounds.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("element\tunit_start\tunit_end\ttx_start\ttx_end\ttss_u_median\n")
    for e, a, b in bounds:
        g.write("%s\t%d\t%d\t%d\t%d\t%d\n" % (e, a, b, a - tss_med, b - tss_med, tss_med))
print("saved: s01_subunit_profile.tsv, s01_element_bounds.tsv (BIN=%d, median TSS within unit %d)" % (BIN, tss_med))
