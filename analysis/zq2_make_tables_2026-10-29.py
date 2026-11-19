# -*- coding: utf-8 -*-
r"""Writes the three tables the figure reads, in one pass.

(1) `HMM_20261029_alu_by_zone.tsv`: the table for supplementary panel (c). Five windows.
     `far`, `mid`, `up` and `gene` must equal the published table (checked);
     `alu` is `alu_own` (bases of the unit's own only): the other four have `own_frac 1.0000`,
     so the five windows are on the same ruler.
(2) `HMM_20261029_paired_molecule_fraction.tsv`: the paired molecule fractions quoted in the legend of main panel (g),
     measured on the FiberHMM calls.
(3) `HMM_20261029_label_rate_ratio.tsv`: the labelling rate ratio quoted in the legend.
     It is measured per platform:
     from the published table `REF2_20261011_m6a_in_unit_by_state_cut508.tsv` (the table main panels (d)(e) draw),
     the depth-weighted mean is ONT 2.19 and PacBio 1.45,
     so it is written as a table that the legend reads, not as one number.

The rulers and gates are copied from `zy7_make_hmm_tables.py`; changing them breaks agreement with the four published tables.
"""
import io
import os
import sys
import gzip
import glob
import collections
import numpy as np
from scipy.stats import wilcoxon

# The data and output roots are environment variables with no default.
#  Paths to files inside the repository are resolved relative to __file__; a_lib is deposited in this folder.
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

MY = _need("RDNA_FHMM_DIR")
FH = os.path.join(MY, "data", "fhmm")
#  This is written to: if it points at the published input folder, those tables are overwritten.
HF = _need("RDNA_FIBERSEQ_DIR")
R = os.path.join(HF, "01_latest_results")
#  The output goes to that folder as well, which collects it in the same place as zy7.
OUTD = R
NL = chr(10)
TB = chr(9)
ST = chr(9733)
CUT, MINCOV, LARGE = 0.508, 1800, 150

REG = [("alu", -2122, -1805), ("far", -1500, -1000), ("mid", -900, -400),
       ("up", -300, -1), ("gene", 0, 118)]
OWNONLY = {"alu"}                     # only this window is counted in the bases of its own unit

sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

SHIFT = {("HG02559", "2"), ("HG03672", "1")}
UNIT, _ = L.load_units()

STATE = {}
with gzip.open(os.path.join(L.D, "a13_unit_state_gmm.tsv.gz"), "rt",
               encoding="utf8") as f:
    h = {k: i for i, k in enumerate(f.readline().rstrip(NL).split(TB))}
    for l in f:
        x = l.rstrip(NL).split(TB)
        try:
            b = float(x[h["beta"]])
        except Exception:                                # noqa: BLE001
            continue
        STATE[(x[h["sample"]], x[h["src"]], x[h["read_id"]],
               x[h["slot"]], x[h["array_pos"]])] = "L" if b <= CUT else "H"

SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FH, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        h = {k: i for i, k in enumerate(f.readline().rstrip(NL).split(TB))}
        for l in f:
            x = l.rstrip(NL).split(TB)
            if len(x) != len(h) or x[h["cls"]] == "tf":
                continue
            SEG[(x[h["sample"]], x[h["read_id"]][:30], x[h["slot"]])].append(
                (x[h["array_pos"]], x[h["cls"]],
                 int(x[h["rel_start"]]), int(x[h["rel_end"]])))
print("state bins %s, FHMM bundles %s" % (format(len(STATE), ","), format(len(SEG), ",")))

RG = collections.defaultdict(lambda: np.zeros(4))   # big, cov, own_cov, own_big
MOL = collections.defaultdict(set)                  # (plat,smp,rid30) -> {L,H}
nread = 0

for smp, src, rid, slot, arr, lo, us in L.iter_reads(UNIT, 0, 0):
    if src == "PACBIO_FIBERSEQ" and (smp in L.PB_EXCL or (smp, slot) in SHIFT):
        continue
    nread += 1
    plat = "ONT" if src.startswith("UCSC") else "PacBio"
    sg = SEG.get((smp, rid[:30], slot))
    if not sg:
        continue
    cod = np.zeros(len(arr), dtype=np.int8)
    for ap, c, a0, b0 in sg:
        u = UNIT.get((smp, slot, ap))
        if u is None:
            continue
        i0 = (u[2] - 1) + a0 - lo
        i1 = (u[2] - 1) + b0 + 1 - lo
        if i1 <= 0 or i0 >= len(arr):
            continue
        cod[max(0, i0):min(len(arr), i1)] = 1 if c == "nuc" else 2
    own = np.full(len(arr), -1, dtype=np.int32)
    tss = {}
    for ap, _t in us:
        s_, e_, t_ = UNIT[(smp, slot, ap)]
        own[s_ - 1 - lo:e_ - lo] = int(ap)
        tss[int(ap)] = t_ - 1 - lo
    ch = np.flatnonzero(np.diff(cod)) + 1
    bb = np.concatenate(([0], ch, [len(cod)]))
    segs = [(bb[i], bb[i + 1], int(cod[bb[i]])) for i in range(len(bb) - 1)]
    big = np.zeros(len(cod), dtype=np.int8)
    for i in range(1, len(segs) - 1):
        a0, e0, c = segs[i]
        if c == 0 or segs[i - 1][2] == 0 or segs[i + 1][2] == 0:
            continue
        if c == 2 and (e0 - a0) > LARGE:
            big[a0:e0] = 1
    for ap, _t in us:
        s_, e_, _t2 = UNIT[(smp, slot, str(ap))]
        a0, b0 = s_ - 1 - lo, e_ - lo
        sub = cod[a0:b0]
        if int((sub > 0).sum()) < MINCOV:
            continue
        st = STATE.get((smp, src, rid, slot, str(ap)))
        if st not in ("L", "H"):
            continue
        MOL[(plat, smp, rid[:30])].add(st)
        t0 = tss[int(ap)]
        for nm, ra, rb in REG:
            i0, i1 = t0 + ra, t0 + rb + 1
            if i0 < 0 or i1 > len(cod):
                continue
            s2 = cod[i0:i1]
            mm = s2 > 0
            if not mm.any():
                continue
            om = mm & (own[i0:i1] == int(ap))
            RG[(plat, smp, st, nm)] += (
                float((big[i0:i1][mm] == 1).sum()), int(mm.sum()),
                float(om.sum()), float((big[i0:i1][om] == 1).sum()))
print("read %s" % format(nread, ","))


def tally(nm):
    own = nm in OWNONLY
    ci, bi = (2, 3) if own else (1, 0)
    out = {}
    for plat in ("ONT", "PacBio"):
        aL, aH = [], []
        for s in sorted({k[1] for k in RG if k[0] == plat and k[3] == nm}):
            vL, vH = RG.get((plat, s, "L", nm)), RG.get((plat, s, "H", nm))
            if vL is None or vH is None or vL[ci] < 500 or vH[ci] < 500:
                continue
            aL.append(vL[bi] / vL[ci])
            aH.append(vH[bi] / vH[ci])
        if len(aL) < 5:
            out[plat] = None
            continue
        aL, aH = np.array(aL), np.array(aH)
        out[plat] = (len(aL), aL.mean(), aH.mean(), aL.mean() - aH.mean(),
                     wilcoxon(aL, aH)[1])
    return out


RES = collections.OrderedDict((nm, tally(nm)) for nm, _a, _b in REG)

# -- 1) per zone, over the five windows
O1 = os.path.join(OUTD, "HMM_20261029_alu_by_zone.tsv")
with io.open(O1, "w", encoding="utf8", newline=NL) as g:
    g.write("# FiberHMM version: bases covered by large MSP(>150bp) per region, individual-level Wilcoxon" + NL)
    g.write("# region section as in zy7_make_hmm_tables.py, with the alu window added" + NL)
    g.write("# windows (relative to TSS) alu -2122..-1805, far -1500..-1000, mid -900..-400 "
            "· up -300..-1 · gene 0..118" + NL)
    g.write("# alu counts only bases of its own unit (own); the other four have own_frac 1.0000, so "
            "the measure is the same. Counting alu plainly, the preceding unit contributes 1.1%" + NL)
    g.write(TB.join(["kind", "platform", "key", "n", "v_L", "v_H", "diff", "p"]) + NL)
    for nm, _a, _b in REG:
        for plat in ("ONT", "PacBio"):
            r = RES[nm][plat]
            if r is None:
                continue
            g.write(TB.join(["region", plat, nm, str(r[0]), "%.5f" % r[1],
                             "%.5f" % r[2], "%.5f" % r[3], "%.3g" % r[4]]) + NL)

# -- 2) the paired-molecule fraction, on the FiberHMM measure
O2 = os.path.join(OUTD, "HMM_20261029_paired_molecule_fraction.tsv")
PAIRPCT = {}
with io.open(O2, "w", encoding="utf8", newline=NL) as g:
    g.write("# fraction of molecules spanning both a low-methylation and a high-methylation unit, FiberHMM version" + NL)
    g.write("# denominator = molecules with at least one unit with a called state (after unit filter called bases>=1800)" + NL)
    g.write("# molecule key = (individual, read_id[:30], slot), same key as zy7" + NL)
    g.write(TB.join(["platform", "mol_all", "mol_pair", "pair_pct"]) + NL)
    for plat in ("ONT", "PacBio"):
        ks = [v for k, v in MOL.items() if k[0] == plat]
        a = len(ks)
        b = sum(1 for v in ks if len(v) == 2)
        PAIRPCT[plat] = 100.0 * b / a if a else float("nan")
        g.write(TB.join([plat, str(a), str(b), "%.4f" % PAIRPCT[plat]]) + NL)

# -- 3) the fold change in marker rate, from the published table
O3 = os.path.join(OUTD, "HMM_20261029_label_rate_ratio.tsv")
P3 = os.path.join(R, "REF2_20261011_m6a_in_unit_by_state_cut508.tsv")
rows3 = []
hdr3 = None
for l in io.open(P3, encoding="utf8"):
    if l.startswith("#"):
        continue
    x = l.rstrip(NL).split(TB)
    if hdr3 is None:
        hdr3 = {k: i for i, k in enumerate(x)}
        continue
    if len(x) == len(hdr3):
        rows3.append(x)
LAB = {}
with io.open(O3, "w", encoding="utf8", newline=NL) as g:
    g.write("# m6A marking rate: low-methylation units ÷ high-methylation units" + NL)
    g.write("# measure = fraction of molecules with m6A at that position (%), depth-weighted unit mean" + NL)
    g.write("# input table = REF2_20261011_m6a_in_unit_by_state_cut508.tsv (the table plotted in main-text panels (d)(e))" + NL)
    g.write(TB.join(["platform", "m6a_L_pct", "m6a_H_pct", "ratio_L_over_H",
                     "n_pos_L", "n_pos_H"]) + NL)
    for plat in ("ONT", "PacBio"):
        o = {}
        for st in ("L", "H"):
            v = [(float(r[hdr3["m6a_molecule_fraction_pct"]]), float(r[hdr3["depth"]]))
                 for r in rows3 if r[hdr3["platform"]] == plat and r[hdr3["state"]] == st]
            a = np.array([x for x, _ in v])
            w = np.array([y for _, y in v])
            o[st] = (float(np.average(a, weights=w)), len(a))
        LAB[plat] = o["L"][0] / o["H"][0]
        g.write(TB.join([plat, "%.4f" % o["L"][0], "%.4f" % o["H"][0],
                         "%.4f" % LAB[plat], str(o["L"][1]), str(o["H"][1])]) + NL)

# ───────────────────────────────── checks
print()
print(ST + ST + " check (1): are the four published windows reproduced as is (alu is an added window)")
pub = {}
for l in io.open(os.path.join(R, "HMM_20261014_by_zone.tsv"), encoding="utf8"):
    if l.startswith("#") or l.startswith("kind"):
        continue
    x = l.rstrip(NL).split(TB)
    if len(x) == 8:
        pub[(x[1], x[2])] = x
ok = bad = 0
for nm, _a, _b in REG:
    if nm == "alu":
        continue
    for plat in ("ONT", "PacBio"):
        r = RES[nm][plat]
        w = pub[(plat, nm)]
        got = "%s %s %s %s %s" % (r[0], "%.5f" % r[1], "%.5f" % r[2],
                                  "%.5f" % r[3], "%.3g" % r[4])
        want = "%s %s %s %s %s" % (w[3], w[4], w[5], w[6], w[7])
        ok += got == want
        bad += got != want
        print("   %-7s %-5s %-30s %s" % (plat, nm, want,
                                         "ok" if got == want else ST + "다름 " + got))
print("   ⇒ " + ST + "같음 %d · 어긋남 %d" % (ok, bad))
if bad:
    raise SystemExit(ST + ST + "the four published windows do not match: stopping.")

print()
print(ST + ST + " check (2): is alu equal to alu_own (compared with zq1)")
for plat, want in (("ONT", "0.02913"), ("PacBio", "0.03454")):
    got = "%.5f" % RES["alu"][plat][3]
    print("   %-7s zq1 alu_own %s vs this alu %s %s"
          % (plat, want, got, "ok" if got == want else ST + "다름"))

print()
print(ST + ST + " check (3): values of the added tables")
for plat in ("ONT", "PacBio"):
    print("   %-7s paired molecules %.2f%% (denominator %s, pairs %s), marking-rate fold %.2f"
          % (plat, PAIRPCT[plat],
             format(len([1 for k in MOL if k[0] == plat]), ","),
             format(len([1 for k, v in MOL.items() if k[0] == plat and len(v) == 2]), ","),
             LAB[plat]))
print()
for p in (O1, O2, O3):
    print("written : %s" % p)
