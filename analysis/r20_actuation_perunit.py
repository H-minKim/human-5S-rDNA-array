# -*- coding: utf-8 -*-
"""R20: actuation per zone, with the zones cut at each unit's own TSS.

Tracks are indexed by position inside the unit, and the TSS varies between units
      (measured 2113 to 2126), so a zone cut at a fixed offset would misplace the gene zone.
Stacking : position i of a unit is stacked at OFF+(i-tss_u), with OFF=2123,
           so OFF points at the TSS of every unit and the zones below are cut at OFF.
The output is written as r20_*.


Unit of quantification, as in the reference (FIRE/fibertools)
  percent actuation = per position : (molecules with an MSP at that position) / (molecules covering that position)
  Note: nucfrac (= N bp/(N+M) bp) is a base-pair fraction and is not the reference unit.

Main platform = ONT; PacBio (HiFi) is reported alongside. The PacBio-only tests (mono/kb, nucA) are not in the main result.

Produced in one run (all in the same unit)
  (1) actuation per zone (per individual, Wilcoxon)
  (2) window dependence : zone window against the whole unit
  (3) depth dependence : strata of molecules per position
  (4) quality confounding : molecules reweighted to the H distribution in 1-point QUAL bins
  (5) untreated background : the same calculation in the untreated control (raw, and after quality matching)
  (6) per-position profile (for the figure)
Output r20_actuation_unified.tsv, r20_actuation_profile.tsv
"""
import io, os, sys, gzip, glob, struct, bisect, collections
import numpy as np
from scipy.stats import wilcoxon

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
from p02_lib import seq2ref

#  This is written to: if it points at the published input folder, those tables are overwritten.
W = _need("RDNA_RESULTMAP_DIR")
R = _need("RDNA_V12_CANON_DIR")
S6 = _need("RDNA_S6_DIR")
#  File-name labels (RFC1918 private range; these are not reachable addresses). NODE_LABELS overrides them.
HOST = {"w3": "ip-172-31-58-58", "iB": "ip-172-31-32-12", "iA": "ip-172-31-57-143"}
if os.environ.get("NODE_LABELS"):
    HOST = dict(kv.split("=", 1) for kv in os.environ["NODE_LABELS"].split(",") if "=" in kv)
SUB = {"UCSC_ONT_FIBERSEQ": "ont_meth", "PACBIO_FIBERSEQ": "pb_meth"}
EXCL = {"PACBIO_FIBERSEQ": {"HG02630", "hg002", "HG02624", "HG03672"}, "UCSC_ONT_FIBERSEQ": set()}
OFF, ULEN = 2123, 2241
REG = [("far", -1500, -1000), ("mid", -900, -400), ("up", -300, -1), ("gene", 0, 118)]
MINCOV = 20
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
_CIG = "MIDNSHP=X"


def read_bam_q(path):
    with gzip.open(path, "rb") as f:
        data = f.read()
    assert data[:4] == b"BAM\x01"
    i = 4
    lt = struct.unpack_from("<i", data, i)[0]; i += 4 + lt
    nr = struct.unpack_from("<i", data, i)[0]; i += 4
    refs = []
    for _ in range(nr):
        ln = struct.unpack_from("<i", data, i)[0]; i += 4
        refs.append(data[i:i + ln - 1].decode()); i += ln + 4
    while i < len(data):
        bs = struct.unpack_from("<i", data, i)[0]
        rec = data[i + 4:i + 4 + bs]; i += 4 + bs
        refid, pos, lrn, mapq, _b, ncig, flag, lseq = struct.unpack_from("<iiBBHHHi", rec, 0)
        o = 32
        name = rec[o:o + lrn - 1].decode(); o += lrn
        cig = []
        for k in range(ncig):
            v = struct.unpack_from("<I", rec, o + 4 * k)[0]
            cig.append((_CIG[v & 15], v >> 4))
        o += 4 * ncig + (lseq + 1) // 2
        qual = np.frombuffer(rec[o:o + lseq], dtype=np.uint8).astype(np.int16)
        yield dict(name=name, flag=flag, ref=refs[refid] if refid >= 0 else None,
                   pos=pos, cigar=cig, qual=qual)


UNITS = collections.defaultdict(list)
with io.open(os.path.join(R, "02_canonical", "MAP_V11.tsv"), encoding="utf8") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for l in f:
        x = l.rstrip("\n").split("\t")
        if x[ix["seg_type"]] == "UNIT":
            UNITS[(x[ix["sample"]], x[ix["slot"]])].append(
                (int(x[ix["unit_start"]]), int(x[ix["unit_end"]]), x[ix["array_pos"]]))

ST = {}
with gzip.open(os.path.join(S6, "a13_unit_state_gmm.tsv.gz"), "rt", encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip("\n").split("\t")
        ST[(x[1], x[0], x[2], x[3], x[4])] = x[10]

# ---------- first pass: the mean QUAL per unit
import sys as _sys
_sys.path.insert(0, _HERE)
import a_lib as _alib
print("pass 0: load TSS per unit", flush=True)
_UN, _ = _alib.load_units()
TSSU = {k: (v[2] - v[0] + 1) for k, v in _UN.items()}
del _UN
print("  units %d, TSS median %d" % (len(TSSU), int(np.median(list(TSSU.values())))), flush=True)

print("pass 1: per-unit QUAL from BAM", flush=True)
QU = {}
BAMS = []
for src, sub in SUB.items():
    for node in HOST:
        BAMS += [(src, q) for q in sorted(glob.glob(os.path.join(R, "08_methylation_s6", node, "s6", sub, "bam", "*.G8v3.bam")))]
for src, bam in BAMS:
    smp = os.path.basename(bam).split(".")[0]
    if smp in EXCL.get(src, set()):
        continue
    for r in read_bam_q(bam):
        if r["flag"] & 0x900 or r["ref"] is None:
            continue
        try:
            sm, slot = r["ref"].rsplit("_h", 1)
        except ValueError:
            continue
        us = sorted(UNITS.get((sm, slot), []))
        if not us:
            continue
        ustart = [u[0] for u in us]
        acc = collections.defaultdict(lambda: [0, 0])
        for sp, g in seq2ref(r).items():
            j = bisect.bisect_right(ustart, g) - 1
            if j < 0:
                continue
            a, b, ap = us[j]
            if a <= g <= b and sp < len(r["qual"]):
                v = acc[ap]; v[0] += int(r["qual"][sp]); v[1] += 1
        for ap, v in acc.items():
            if v[1] >= 500:
                QU[(smp, src, r["name"], slot, ap)] = v[0] / v[1]
print("  QUAL units %d" % len(QU), flush=True)

# ---------- second pass: stream the track and accumulate per position
NB = ULEN + 2
P_SMP = collections.defaultdict(lambda: np.zeros((NB, 2)))      # (src,smp,state)
P_Q = collections.defaultdict(lambda: np.zeros((NB, 2)))        # (src,state,qbin)
QCNT = collections.Counter()                                     # (source, state, quality bin) -> number of molecules
NOTSS = collections.Counter() # units with no TSS found
print("pass 2: accumulate tracks", flush=True)
for p in sorted(glob.glob(os.path.join(S6, "a01_tracks", "*.tsv.gz"))):
    with gzip.open(p, "rt", encoding="utf8") as f:
        ix = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
        for l in f:
            x = l.rstrip("\n").split("\t")
            smp, src = x[ix["sample"]], x[ix["src"]]
            if smp in EXCL.get(src, set()):
                continue
            key = (smp, src, x[ix["read_id"]], x[ix["slot"]], x[ix["array_pos"]])
            st = ST.get(key, "U")
            if st not in ("L", "H"):
                continue
            tss_u = TSSU.get((smp, x[ix["slot"]], x[ix["array_pos"]]))
            if tss_u is None:
                NOTSS[src] += 1
                continue
            a = np.frombuffer(x[ix["track"]].encode(), dtype=np.uint8)
            cov = (a != ord(".")).astype(np.float64)
            msp = (a == ord("M")).astype(np.float64)
            # position i, 1-based, maps to OFF + (i - tss_u); anything outside the array is clipped
            sh = OFF - tss_u
            j0 = 1 + sh                      # the bin position 1 falls in
            k0 = max(0, 1 - j0)   # lower bound 1: nothing is written to bin 0, the same boundary as r18
            k1 = min(len(a), ULEN + 1 - j0)   # upper bound ULEN (2241), the same boundary as r18
            if k1 <= k0:
                continue
            d0, d1 = j0 + k0, j0 + k1
            A = P_SMP[(src, smp, st)]
            A[d0:d1, 0] += cov[k0:k1]; A[d0:d1, 1] += msp[k0:k1]
            q = QU.get(key)
            if q is not None:
                qb = int(round(q))
                B = P_Q[(src, st, qb)]
                B[d0:d1, 0] += cov[k0:k1]; B[d0:d1, 1] += msp[k0:k1]
                QCNT[(src, st, qb)] += 1
print("  individual×state %d, QUAL bins %d, no TSS %s" % (len(P_SMP), len(QCNT), dict(NOTSS)), flush=True)


def act(A, a0, b0, mincov=MINCOV):
    m = A[a0:b0, 0] >= mincov
    if m.sum() < 30:
        return float("nan")
    return float((A[a0:b0, 1][m] / A[a0:b0, 0][m]).mean())


rows = []
SRCS = (("UCSC_ONT_FIBERSEQ", "ONT"), ("PACBIO_FIBERSEQ", "PacBio"))

print("\n" + "=" * 100)
print("(1) percent actuation by zone (individual level)")
print("%-7s %-5s %6s %9s %9s %11s %11s" % ("src", "구역", "개체", "act_L", "act_H", "L-H", "p"))
for src, s6 in SRCS:
    smps = sorted({k[1] for k in P_SMP if k[0] == src})
    for rn, ra, rb in REG:
        a0, b0 = OFF + ra, OFF + rb + 1
        dL, dH = [], []
        for smp in smps:
            AL = P_SMP.get((src, smp, "L")); AH = P_SMP.get((src, smp, "H"))
            if AL is None or AH is None:
                continue
            vl, vh = act(AL, a0, b0), act(AH, a0, b0)
            if np.isfinite(vl) and np.isfinite(vh):
                dL.append(vl); dH.append(vh)
        if len(dL) < 8:
            continue
        dL, dH = np.array(dL), np.array(dH)
        _, pv = wilcoxon(dL, dH)
        print("%-7s %-5s %6d %9.4f %9.4f %+11.4f %11.2e" %
              (s6, rn, len(dL), dL.mean(), dH.mean(), (dL - dH).mean(), pv))
        rows.append(("region", s6, rn, len(dL), round(float(dL.mean()), 5),
                     round(float(dH.mean()), 5), round(float((dL - dH).mean()), 5), "%.3g" % pv))

print("\n" + "=" * 100)
print("(2) window dependence: zone window vs whole unit")
print("%-7s %-12s %6s %11s %11s" % ("src", "창", "개체", "L-H", "p"))
for src, s6 in SRCS:
    smps = sorted({k[1] for k in P_SMP if k[0] == src})
    for wn, a0, b0 in (("zone window", OFF - 1500, OFF + 119), ("whole unit", 1, ULEN + 1)):
        dL, dH = [], []
        for smp in smps:
            AL = P_SMP.get((src, smp, "L")); AH = P_SMP.get((src, smp, "H"))
            if AL is None or AH is None:
                continue
            vl, vh = act(AL, a0, b0), act(AH, a0, b0)
            if np.isfinite(vl) and np.isfinite(vh):
                dL.append(vl); dH.append(vh)
        if len(dL) < 8:
            continue
        dL, dH = np.array(dL), np.array(dH)
        _, pv = wilcoxon(dL, dH)
        print("%-7s %-12s %6d %+11.4f %11.2e" % (s6, wn, len(dL), (dL - dH).mean(), pv))
        rows.append(("window", s6, wn, len(dL), round(float(dL.mean()), 5),
                     round(float(dH.mean()), 5), round(float((dL - dH).mean()), 5), "%.3g" % pv))

print("\n" + "=" * 100)
print("(3) stratum dependence: strata of molecules per position (whole-unit window)")
print("%-7s %-12s %10s %10s %11s" % ("src", "molecule count stratum", "act_L", "act_H", "L-H"))
CB = [(20, 50), (50, 100), (100, 200), (200, 500), (500, 10 ** 9)]
for src, s6 in SRCS:
    TL = np.zeros((NB, 2)); TH = np.zeros((NB, 2))
    for k, A in P_SMP.items():
        if k[0] != src:
            continue
        (TL if k[2] == "L" else TH)[:] += A
    for lo, hi in CB:
        m = (TL[:, 0] >= lo) & (TL[:, 0] < hi) & (TH[:, 0] >= lo) & (TH[:, 0] < hi)
        if m.sum() < 30:
            continue
        vl = float((TL[m, 1] / TL[m, 0]).mean()); vh = float((TH[m, 1] / TH[m, 0]).mean())
        print("%-7s %-12s %10.4f %10.4f %+11.4f   (positions %d)" %
              (s6, "%d-%s" % (lo, "" if hi > 10 ** 8 else str(hi)), vl, vh, vl - vh, int(m.sum())))
        rows.append(("by_cov", s6, "%d-%s" % (lo, "" if hi > 10 ** 8 else str(hi)), int(m.sum()),
                     round(vl, 5), round(vh, 5), round(vl - vh, 5), ""))

print("\n" + "=" * 100)
print("(4) quality confounding: L reweighted to the H distribution in 1-point QUAL bins (whole-unit window)")
print("%-7s %11s %11s %9s" % ("src", "before matching", "맞춤 후", "retention"))
for src, s6 in SRCS:
    qbs = sorted({k[2] for k in QCNT if k[0] == src})
    TL = np.zeros((NB, 2)); TH = np.zeros((NB, 2))
    WL = np.zeros((NB, 2)); WH = np.zeros((NB, 2))
    for qb in qbs:
        nl = QCNT.get((src, "L", qb), 0); nh = QCNT.get((src, "H", qb), 0)
        AL = P_Q.get((src, "L", qb)); AH = P_Q.get((src, "H", qb))
        if AL is None or AH is None:
            continue
        TL += AL; TH += AH
        if nl >= 10 and nh >= 10:
            WL += AL * (nh / nl)           # L reweighted to H's QUAL distribution
            WH += AH
    def _a(Z):
        m = Z[:, 0] >= MINCOV
        return float((Z[m, 1] / Z[m, 0]).mean()) if m.sum() > 30 else float("nan")
    raw = _a(TL) - _a(TH); adj = _a(WL) - _a(WH)
    print("%-7s %11.4f %11.4f %8.0f%%" % (s6, raw, adj, 100 * adj / raw))
    rows.append(("qual", s6, "matched", 0, round(raw, 5), round(adj, 5),
                 round(100 * adj / raw, 1), ""))

with io.open(os.path.join(W, "r20_actuation_unified.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("kind\tplatform\tkey\tn\tv_L\tv_H\tdiff\tp\n")
    for x in rows:
        g.write("\t".join(str(q) for q in x) + "\n")

# ---------- the profile by position
prof = []
for src, s6 in SRCS:
    for st in ("L", "H"):
        T = np.zeros((NB, 2))
        for k, A in P_SMP.items():
            if k[0] == src and k[2] == st:
                T += A
        for i in range(1, ULEN + 1):
            if T[i, 0] >= 200:
                prof.append((s6, st, i, int(T[i, 0]), T[i, 1] / T[i, 0]))
with io.open(os.path.join(W, "r20_actuation_profile.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("platform\tstate\tunit_pos\tn_mol\tactuation\n")
    for r in prof:
        g.write("%s\t%s\t%d\t%d\t%.5f\n" % r)
print("\nsaved: r20_actuation_unified.tsv, r20_actuation_profile.tsv")
