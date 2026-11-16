# -*- coding: utf-8 -*-
# The v11 array reference: by the same rule and the same bins as the dispersed loci
#   State = the CpG on NTS1 (319-1178), NTS2 and the promoter (1552-2121) and the gene (2122-2240), at least 10 of them
#          < 0.3 low, > 0.7 high, otherwise mid   (same as struct_aln3.py, ML >= 128)
#   Bin   = (unit_pos - 1) // 50   (0..44)   - coordinate frame canon_ctx == unit_canon (0 mismatches over 2,240 bp)
#   Position = distal end 1-10, proximal end 1-10, interior (SUBUNIT prox_idx / dist_idx)
# Input out/CALLS.tsv.gz (CpG calls taken from the deduplicated v11 READMAP_METH) and out/SUBUNIT.tsv.gz
import io, os, gzip, collections, sys

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
OUT = os.path.join(W, "dispersed_vs_array")
SUM = os.path.join(V, "03_maps", "g8canonv10", "CONSENSUS_SUM_G8_full.tsv")
BRK = set()
for i, L in enumerate(io.open(SUM, encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    if len(f) >= 8 and int(f[7]) > 0:
        BRK.add((f[0], f[1]))
POS = {}
for i, L in enumerate(gzip.open(os.path.join(W, "out", "SUBUNIT.tsv.gz"), "rt", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    pi, di = int(f[3]), int(f[4])
    POS[(f[0], f[1], f[2], pi)] = "distal1_10" if di <= 10 else ("proximal1_10" if pi <= 10 else "interior")
sys.stderr.write("SUBUNIT 위치 %d\n" % len(POS))


def spec(u):
    return 319 <= u <= 1178 or 1552 <= u <= 2240


BINS = collections.defaultdict(lambda: [0, 0])     # (pos, state, bin) -> n, meth
NU = collections.Counter()
cur, seen, skip = None, set(), False
buf = collections.defaultdict(list)


def flush(key, b):
    for pi, calls in b.items():
        cs = sum(1 for u, m in calls if spec(u))
        if cs < 10:
            st = "NA"
        else:
            beta = sum(m for u, m in calls if spec(u)) / float(cs)
            st = "low" if beta < 0.3 else ("high" if beta > 0.7 else "mid")
        pos = POS.get(key + (pi,), "unknown")
        NU[(pos, st)] += 1
        for u, m in calls:
            for pp in (pos, "ALL"):
                v = BINS[(pp, st, (u - 1) // 50)]
                v[0] += 1
                v[1] += m


n = 0
for L in gzip.open(os.path.join(W, "out", "CALLS.tsv.gz"), "rt", encoding="utf8"):
    if L[0] == "s":
        continue
    f = L.rstrip("\n").split("\t")
    if len(f) < 6:
        continue
    rid = (f[0], f[1], f[2])
    if rid != cur:
        if cur is not None and not skip:
            flush(cur, buf)
        buf = collections.defaultdict(list)
        cur = rid
        skip = (rid in seen) or ((f[0], f[1]) in BRK)
        seen.add(rid)
    if not skip:
        u = int(f[4])
        if 1 <= u <= 2240:
            buf[int(f[3])].append((u, 1 if int(f[5]) >= 128 else 0))
    n += 1
    if n % 50000000 == 0:
        sys.stderr.write("  %d 백만\n" % (n // 1000000))
        sys.stderr.flush()
if cur is not None and not skip:
    flush(cur, buf)
with io.open(os.path.join(OUT, "ARRAY_V11_BINS.tsv"), "w", encoding="utf8", newline="\n") as O:
    O.write("pos\tstate\tbin\tn\tmeth\n")
    for (p, s, b), (nn, k) in sorted(BINS.items()):
        O.write("%s\t%s\t%d\t%d\t%d\n" % (p, s, b, nn, k))
with io.open(os.path.join(OUT, "ARRAY_V11_UNITS.tsv"), "w", encoding="utf8", newline="\n") as O:
    O.write("pos\tstate\tn_units\n")
    for (p, s), c in sorted(NU.items()):
        O.write("%s\t%s\t%d\n" % (p, s, c))
print("유닛 상태 :", dict(collections.Counter({s: sum(c for (p, s2), c in NU.items() if s2 == s) for s in ("low", "mid", "high", "NA")})))
print("saved ARRAY_V11_BINS.tsv")
