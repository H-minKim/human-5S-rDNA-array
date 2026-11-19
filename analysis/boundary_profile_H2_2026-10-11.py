# -*- coding: utf-8 -*-
# Boundary profile - H2 gate and cut 0.508
#
#   Why the H2 gate and the 0.508 cut
#     (a) the `RU.tsv.gz` and `CALLS.tsv.gz` on disk have no `data_source` column
#        (`load_readunit_v2.py` writes it, but the files are the v1 generation), so the data source cannot be filtered from them.
#        Of 303 samples, 289 are HPRC only, 9 are 1kGP and HPRC name collisions,
#        1 is GIAB_SOMATIC (the cancer sample HG008-T), 3 are GIAB_ONT and 1 is 1kGP only.
#        On clean windows, only 90.1% / 91.0% fall inside H2 of the registry.
#     (b) the cut in `boundary_profile_wide.py` is 0.5, while the project constant is 0.508.
#
#   Everything is produced in one pass - 788 MB is read only once
#     (1) BOUNDARY_WIDE on the new ruler (H2 and 0.508)
#     (2) the same on the old ruler (no gate, 0.5) as well, for the reproduction test
#     (3) the full table of clean 3+3 windows (new ruler) - for reselecting the numerator of Supplementary Figure 7
#     (4) the mean over the contrast-passing windows alone - to tell whether the grey curve is shallow because of (i) overlap or (ii) weak windows
#     (5) the transition position of each window - turning "AluY 4/4" into a distribution
#
#   Everything is written under new names, so the outputs of `boundary_profile_wide.py` are not overwritten.
import collections
import csv
import gzip
import io
import os
import sys

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    The roots of the data this script reads are listed in the
    "Running this code" section of README.md."""
    import os as _o, sys as _s
    _x = _o.environ.get(_v)
    if not _x:
        _s.exit(
            "environment variable %s is not set; it is required by %s"
            % (_v, _o.path.basename(__file__))
            + "\n(a module stops the script that imports it)."
            + "\nThe deposited copy has no default; see README.md"
            + " section \"Running this code\"."
        )
    return _x


W = _need("RDNA_WORK_ROOT")
V = os.path.join(_need("RDNA_VERSION_ROOT"), "08_versions", "v11")
T = _need("RDNA_V11_ROOT")
RU = os.path.join(W, "out", "RU.tsv.gz")
CA = os.path.join(W, "out", "CALLS.tsv.gz")
SUM = os.path.join(V, "03_maps", "g8canonv10", "CONSENSUS_SUM_G8_full.tsv")
REG = os.path.join(T, "02_tables", "TIERS_H0H3_SLOTS_2026-10-05.tsv")
OUT = os.path.join(W, "out_H2_2026-10-11")
os.makedirs(OUT, exist_ok=True)

UL, MLTHR, MINCPG = 2240, 128, 20
NBEF, NAFT = 3, 3
NU = NBEF + NAFT
NL = chr(10)
#  Two rulers side by side: (name, unit hard-call cut, whether the H2 gate is applied)
RULERS = [("H2_508", 0.508, True), ("old_500", 0.500, False)]
#  Contrast: the minimum of the three high-side units >= HI and the maximum of the three low-side units <= LO
LO, HI = 0.25, 0.70
#  Transition-position ruler: at least RF of the RW CpG on either side are in that state
RW, RF = 20, 0.90
EL = [("AluY", 1, 318), ("spacer-1", 319, 1178), ("GA-rich", 1179, 1551),
      ("spacer-2", 1552, 2031), ("5'-flank", 2032, 2121), ("5S gene", 2122, 2240)]


def elem(p):
    for nm, a, b in EL:
        if a <= p <= b:
            return nm
    return "?"


H2 = {(r["sample"], r["slot"]) for r in
      csv.DictReader(io.open(REG, encoding="utf8"), delimiter=chr(9)) if r["tier"] == "H2"}
sys.stderr.write("H2 슬롯 %d%s" % (len(H2), NL))

BRK = set()
for i, L in enumerate(io.open(SUM, encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) >= 8 and int(f[7]) > 0:
        BRK.add((f[0], f[1]))

#  Unit hard call and beta - separately for each ruler
CLS = {k: {} for k, _, _ in RULERS}
BETA = {}
for i, L in enumerate(gzip.open(RU, "rt", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 10 or int(f[7]) < MINCPG or (f[0], f[1]) in BRK:
        continue
    k3 = (f[0], f[1], f[2], int(f[3]))
    b = float(f[9])
    BETA[k3] = b
    for nm, thr, gate in RULERS:
        if gate and (f[0], f[1]) not in H2:
            continue
        CLS[nm][k3] = 1 if b >= thr else 0
for nm, _, _ in RULERS:
    sys.stderr.write("%s 유닛 %d%s" % (nm, len(CLS[nm]), NL))

WIN = {nm: {"H2L": collections.defaultdict(lambda: [0, 0]),
            "L2H": collections.defaultdict(lambda: [0, 0])} for nm, _, _ in RULERS}
#  (4) contrast-passing windows only (new ruler only)
WINC = {"H2L": collections.defaultdict(lambda: [0, 0]),
        "L2H": collections.defaultdict(lambda: [0, 0])}
CTRL = {nm: {0: collections.defaultdict(lambda: [0, 0]),
             1: collections.defaultdict(lambda: [0, 0])} for nm, _, _ in RULERS}
N = collections.Counter()
ROWS = []        # (3) the full set of clean windows (new ruler)
STEPS = []       # (5) transition positions
NEW = RULERS[0][0]


def step_pos(seq_pos, seq_m, A):
    """Transition positions from the raw CpGs. seq_pos and seq_m keep the order inside the window (walking distal to proximal)."""
    n = len(seq_m)
    if n < 2 * RW:
        return None
    B = 1 - A
    back = [0] * n
    fwd = [0] * n
    c = 0
    for i in range(n):
        c += 1 if seq_m[i] == A else 0
        if i >= RW:
            c -= 1 if seq_m[i - RW] == A else 0
        back[i] = c / float(min(i + 1, RW))
    c = 0
    for i in range(n - 1, -1, -1):
        c += 1 if seq_m[i] == B else 0
        if i + RW < n:
            c -= 1 if seq_m[i + RW] == B else 0
        fwd[i] = c / float(min(n - i, RW))
    la = [i for i in range(n) if back[i] >= RF]
    if not la:
        return None
    ia = max(la)
    fb = [i for i in range(ia + 1, n) if fwd[i] >= RF]
    if not fb:
        return None
    ib = min(fb)
    return ia, ib


def process(read, key3):
    for nm, thr, gate in RULERS:
        cl = CLS[nm]
        idxs = sorted(read, reverse=True)
        if len(idxs) < 3:
            continue
        c = [cl.get(key3 + (i,)) for i in idxs]
        #  Control: units inside a block whose two neighbours are of the same class
        for a in range(1, len(idxs) - 1):
            if idxs[a - 1] - idxs[a] != 1 or idxs[a] - idxs[a + 1] != 1:
                continue
            if c[a] is None or c[a - 1] != c[a] or c[a + 1] != c[a]:
                continue
            d = CTRL[nm][c[a]]
            for p, m in read[idxs[a]]:
                v = d[p]
                v[0] += 1
                v[1] += m
        #  Transition windows
        for a in range(len(idxs) - 1):
            if c[a] is None or c[a + 1] is None or c[a] == c[a + 1]:
                continue
            lo, hi = a - NBEF + 1, a + NAFT
            if lo < 0 or hi >= len(idxs):
                continue
            if any(idxs[t] - idxs[t + 1] != 1 for t in range(lo, hi)):
                continue
            if any(c[t] != c[a] for t in range(lo, a + 1)):
                continue
            if any(c[t] != c[a + 1] for t in range(a + 1, hi + 1)):
                continue
            kind = "H2L" if c[a] == 1 else "L2H"
            N["%s|%s" % (nm, kind)] += 1
            units = [idxs[t] for t in range(lo, hi + 1)]
            bts = [BETA.get(key3 + (u,)) for u in units]
            for j, t in enumerate(range(lo, hi + 1)):
                base = j * UL
                for p, m in read[idxs[t]]:
                    v = WIN[nm][kind][base + p]
                    v[0] += 1
                    v[1] += m
            if nm != NEW:
                continue
            #  ---- new ruler only: window table, contrast, transition position
            ncs = [len(read[u]) for u in units]
            if None in bts:
                continue
            hi_s = min(bts[:3]) if kind == "H2L" else min(bts[3:])
            lo_s = max(bts[3:]) if kind == "H2L" else max(bts[:3])
            okc = (hi_s >= HI) and (lo_s <= LO)
            ROWS.append((kind,) + key3 + (units[0], units[-1],
                                          ",".join(str(x) for x in ncs),
                                          ",".join("%.4f" % x for x in bts),
                                          "%.4f" % (hi_s - lo_s), 1 if okc else 0))
            if okc:
                for j, t in enumerate(range(lo, hi + 1)):
                    base = j * UL
                    for p, m in read[idxs[t]]:
                        v = WINC[kind][base + p]
                        v[0] += 1
                        v[1] += m
            #  (5) transition position: the window is unrolled into a CpG sequence and split on the raw ruler
            sp, sm, su = [], [], []
            for j, t in enumerate(range(lo, hi + 1)):
                for p, m in sorted(read[idxs[t]]):
                    sp.append(p)
                    sm.append(m)
                    su.append(j - (NBEF - 1))
            r = step_pos(sp, sm, c[a])
            if r:
                ia, ib = r
                STEPS.append((kind, key3[0], key3[1], su[ia], sp[ia], elem(sp[ia]),
                              su[ib], sp[ib], elem(sp[ib]), ib - ia - 1,
                              1 if okc else 0))


cur, buf, seen, skip, n = None, collections.defaultdict(list), set(), False, 0
for L in gzip.open(CA, "rt", encoding="utf8"):
    if L[0] == "s":
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 6:
        continue
    rid = (f[0], f[1], f[2])
    if rid != cur:
        if cur is not None and not skip:
            process(buf, cur)
        buf = collections.defaultdict(list)
        cur = rid
        skip = rid in seen
        seen.add(rid)
    if not skip and ((f[0], f[1], f[2], int(f[3])) in BETA):
        buf[int(f[3])].append((int(f[4]), 1 if int(f[5]) >= MLTHR else 0))
    n += 1
    if n % 50000000 == 0:
        sys.stderr.write("  %d 백만%s" % (n // 1000000, NL))
        sys.stderr.flush()
if cur is not None and not skip:
    process(buf, cur)
print("number of windows: %s" % dict(N))


def dump(name, win, ctrl):
    o = io.open(os.path.join(OUT, name), "w", encoding="utf8", newline=NL)
    o.write("kind\tg\tunit_off\tunit_pos\tn\tbeta\tn_LL\tb_LL\tn_HH\tb_HH\tscore" + NL)
    k = 0
    for kind in ("H2L", "L2H"):
        for g in range(1, NU * UL + 1):
            v = win[kind].get(g)
            if not v or v[0] < 20:
                continue
            up = (g - 1) % UL + 1
            off = (g - 1) // UL - (NBEF - 1)
            cl, ch = ctrl[0].get(up), ctrl[1].get(up)
            if not cl or not ch or cl[0] < 50 or ch[0] < 50:
                continue
            bl, bh = cl[1] / cl[0], ch[1] / ch[0]
            b = v[1] / v[0]
            s = (b - bl) / (bh - bl) if bh - bl >= 0.15 else ""
            o.write("%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%.4f\t%d\t%.4f\t%s%s"
                    % (kind, g, off, up, v[0], b, cl[0], bl, ch[0], bh,
                       ("%.4f" % s) if s != "" else "", NL))
            k += 1
    o.close()
    print("  %-34s %s lines" % (name, format(k, ",")))


for nm, _, _ in RULERS:
    dump("BOUNDARY_WIDE_%s.tsv" % nm, WIN[nm], CTRL[nm])
dump("BOUNDARY_WIDE_%s_대비통과만.tsv" % NEW, WINC, CTRL[NEW])

with io.open(os.path.join(OUT, "CLEAN_WINDOWS_%s.tsv" % NEW), "w", encoding="utf8",
             newline=NL) as o:
    o.write("kind\tsample\tslot\tread_id\tunit_first\tunit_last\tncpg_per_unit"
            "\tbeta_per_unit\tgap\tcontrast_ok" + NL)
    for r in ROWS:
        o.write("\t".join(str(x) for x in r) + NL)
print("  clean windows %s" % format(len(ROWS), ","))

with io.open(os.path.join(OUT, "STEP_POS_%s.tsv" % NEW), "w", encoding="utf8",
             newline=NL) as o:
    o.write("kind\tsample\tslot\tunit_off_before\tunit_pos_before\telem_before"
            "\tunit_off_after\tunit_pos_after\telem_after\tcpg_between\tcontrast_ok" + NL)
    for r in STEPS:
        o.write("\t".join(str(x) for x in r) + NL)
print("  transition positions %s" % format(len(STEPS), ","))
print("output: %s" % OUT)
