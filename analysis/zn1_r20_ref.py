# -*- coding: utf-8 -*-
"""r20 with the reference threshold; F2 (a)(b) read this table.

The definition of `actuation`
   r20  : bases covered by MSPs / called bases               (no threshold)
   here : bases covered by large MSPs (>150 bp) / called bases
        *"The cutoff for identifying di-nucleosomes and large MSPs was 210 and 150"*

Everything else is as in r20, so the two can be compared at the same positions
   - zones : far -1500..-1000, mid -900..-400, up -300..-1, gene 0..+118  (from the TSS)
   - per-position profile : positions 1..2241 of the unit, each unit stacked with its own TSS at OFF=2123
   - test : Wilcoxon signed-rank per individual (low against high methylation)
   - excluded : 2 PacBio fibertools individuals, 2 PacBio slots with shifted positions
The layout is that of the main table a23 (canonical coordinates).
"""
import io
import os
import sys
import collections
import numpy as np
from scipy.stats import wilcoxon

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
W = _need("RDNA_RESULTMAP_DIR")
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

ULEN, OFF, LARGE = 2241, 2123, 150
REG = [("far", -1500, -1000), ("mid", -900, -400), ("up", -300, -1), ("gene", 0, 118)]
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
UNIT, _ = L.load_units()
GMM = L.load_state_gmm()

# (platform, sample, state, region) -> [large bp, called bp]
RG = collections.defaultdict(lambda: np.zeros(2))
# (platform, state, unit position) -> [large bp, called bp, molecules]
PR = collections.defaultdict(lambda: np.zeros(3))
nread = 0
for smp, src, rid, slot, arr, lo, us in L.iter_reads(UNIT, 0, 0):
    plat = "ONT" if src.startswith("UCSC") else "PacBio"   # the names the figure looks for
    if src == "PACBIO_FIBERSEQ" and (smp in L.PB_EXCL or (smp, slot) in SHIFT):
        continue
    nread += 1
    cod = np.zeros(len(arr), dtype=np.int8)
    cod[arr == ord("N")] = 1
    cod[arr == ord("M")] = 2
    ch = np.flatnonzero(np.diff(cod)) + 1
    bb = np.concatenate(([0], ch, [len(cod)]))
    segs = [(bb[i], bb[i + 1], cod[bb[i]]) for i in range(len(bb) - 1)]
    big = np.zeros(len(arr), dtype=bool)
    for i in range(1, len(segs) - 1):
        a0, e0, c = segs[i]
        if c == 2 and segs[i - 1][2] != 0 and segs[i + 1][2] != 0 and (e0 - a0) > LARGE:
            big[a0:e0] = True
    for ap, _t in us:
        s_, e_, tss = UNIT[(smp, slot, str(ap))]
        st = GMM.get((smp, src, rid, slot, str(ap)), "NA")
        if st not in ("L", "H"):
            continue
        a0, b0 = s_ - 1 - lo, e_ - lo
        sc, sb = cod[a0:b0], big[a0:b0]
        if len(sc) == 0:
            continue
        tss_u = tss - s_ + 1
        up = np.arange(1, len(sc) + 1)
        rel = up - tss_u
        ok = sc > 0
        for nm, r0, r1 in REG:
            m = ok & (rel >= r0) & (rel <= r1)
            k = int(m.sum())
            if k:
                v = RG[(plat, smp, st, nm)]
                v[0] += int(sb[m].sum()); v[1] += k
        # the profile by position, with the TSS aligned to OFF
        tgt = OFF + rel
        inr = ok & (tgt >= 1) & (tgt <= ULEN)
        if inr.any():
            tt = tgt[inr]; bbv = sb[inr]
            for t_, b_ in zip(tt, bbv):
                p = PR[(plat, st, int(t_))]
                p[0] += int(b_); p[1] += 1
            for t_ in np.unique(tt):
                PR[(plat, st, int(t_))][2] += 1
print("reads %d, zone bins %d, position bins %d" % (nread, len(RG), len(PR)))

# ───────── unified: per zone, with a Wilcoxon over individuals
OUT1 = os.path.join(W, "r20_actuation_unified_REF.tsv")
rows = []
for plat in ("ONT", "PacBio"):
    for nm, _a, _b in REG:
        aL, aH = [], []
        smps = {s for (p_, s, t_, r_) in RG if p_ == plat and r_ == nm}
        for s in sorted(smps):
            vL, vH = RG.get((plat, s, "L", nm)), RG.get((plat, s, "H", nm))
            if vL is None or vH is None or vL[1] < 500 or vH[1] < 500:
                continue
            aL.append(vL[0] / vL[1]); aH.append(vH[0] / vH[1])
        if len(aL) < 5:
            continue
        aL, aH = np.array(aL), np.array(aH)
        _s, p = wilcoxon(aL, aH)
        rows.append(("region", plat, nm, len(aL), aL.mean(), aH.mean(),
                     aL.mean() - aH.mean(), p))
with io.open(OUT1, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# r20 with the reference threshold: actuation = bases covered by large MSP(>150bp) / called bases"
            + chr(10))
    g.write("# zones, test and exclusions are the same as r20 (individual-level Wilcoxon)" + chr(10))
    g.write(chr(9).join(["kind", "platform", "key", "n", "v_L", "v_H", "diff", "p"]) + chr(10))
    for r in rows:
        g.write(chr(9).join([r[0], r[1], r[2], str(r[3]), "%.5f" % r[4], "%.5f" % r[5],
                             "%.5f" % r[6], "%.3g" % r[7]]) + chr(10))

print()
print("=== by zone (large MSP >150 bp)")
print("   %-5s %-6s %4s %9s %9s %9s %10s" % ("platform", "구역", "n", "저메틸", "고메틸", "차", "p"))
for r in rows:
    print("   %-5s %-6s %4d %9.5f %9.5f %9.5f %10.2e%s"
          % (r[1], r[2], r[3], r[4], r[5], r[6], r[7], "" if r[7] < 0.01 else "  n.s."))

# ───────── the profile, by position
OUT2 = os.path.join(W, "r20_actuation_profile_REF.tsv")
with io.open(OUT2, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# actuation by position = bases covered by large MSP(>150bp) / called bases "
            "(each unit's own TSS aligned to %d)" % OFF + chr(10))
    g.write(chr(9).join(["platform", "state", "unit_pos", "n_mol", "actuation"]) + chr(10))
    n2 = 0
    for (plat, st, upx), v in sorted(PR.items()):
        if v[1] < 200:
            continue
        g.write(chr(9).join([plat, st, str(upx), str(int(v[2])),
                             "%.5f" % (v[0] / v[1])]) + chr(10))
        n2 += 1
print()
print("position rows %d -> %s" % (n2, os.path.basename(OUT2)))
print("zone rows %d -> %s" % (len(rows), os.path.basename(OUT1)))
