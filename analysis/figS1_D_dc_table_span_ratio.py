# -*- coding: utf-8 -*-
# d is counted as the number of distinct array_pos values a read covers.
#   Counting rows instead double-counts a read across the two placements and configurations.
import io, os, gzip, glob, collections, statistics

# -- (repo) paths come from the environment --------------------------------
#   The roots are required; the output and table directories fall back to this file's folder.
V = os.environ["V_ROOT"]

MAP = os.path.join(V, r"03_maps/g8canonv10/CONSENSUS_MAP_G8_full.tsv")
HSN = os.path.join(V, r"02_canonical/hap_seq_node_v8.tsv")
RAW = os.path.join(V, r"08_methylation/raw")
#   DC_OUT sets the location, defaulting to this file's folder.
OUT = os.environ.get("DC_OUT",
                     os.path.join(os.path.dirname(os.path.abspath(__file__)), "DC_TABLE.tsv"))

c = collections.Counter()
f = io.open(MAP, encoding="utf8")
h = f.readline().rstrip("\n").split("\t")
ix = {k: i for i, k in enumerate(h)}
for L in f:
    r = L.rstrip("\n").split("\t")
    if r[ix["seg_type"]] == "UNIT":
        c[(r[0], r[1])] += 1

BODY = {"arr": "HPRC", "hiarr": "hifiasm", "t1": "span", "c2seq": "cons2", "pat": "T2T"}
src = {}
for L in io.open(HSN, encoding="utf8"):
    r = L.rstrip("\n").split("\t")
    if len(r) >= 3:
        src[(r[0], r[1])] = BODY.get(r[2].split("__")[-1].split("_")[0], "?")

pos = collections.defaultdict(lambda: collections.defaultdict(set))   # (smp,slot) -> read -> {array_pos}
kinds = collections.defaultdict(lambda: collections.defaultdict(set))  # read -> {AO/JX}
for path in sorted(glob.glob(os.path.join(RAW, "METH_READUNIT_*.tsv.gz"))):
    tag = "JX" if "_jx" in os.path.basename(path) else "AO"
    with gzip.open(path, "rt", encoding="utf8") as g:
        hh = g.readline().rstrip("\n").split("\t")
        j = {k: i for i, k in enumerate(hh)}
        for L in g:
            r = L.rstrip("\n").split("\t")
            k = (r[0], r[j["slot"]])
            pos[k][r[j["read_id"]]].add(int(r[j["array_pos"]]))
            kinds[k][r[j["read_id"]]].add(tag)
    print("  read %s" % os.path.basename(path), flush=True)

rows = []
for k in sorted(c):
    if k not in pos or not pos[k]:
        continue
    best = max(pos[k].items(), key=lambda kv: len(kv[1]))
    d = len(best[1])
    rows.append((k[0], k[1], src.get(k, "?"), c[k], d, d / float(c[k]), len(pos[k]),
                 "+".join(sorted(kinds[k][best[0]]))))

with io.open(OUT, "w", encoding="utf8", newline="\n") as g:
    g.write("sample\tslot\tsource\tc_units\td_longest_read_units\td_over_c\tn_reads\tlongest_read_kind\n")
    for r in rows:
        g.write("%s\t%s\t%s\t%d\t%d\t%.4f\t%d\t%s\n" % r)
print("table %d rows -> %s" % (len(rows), OUT))

print()
print("## d/c by source  (d = number of distinct units covered by the longest read)")
print("| source | haplotypes | c median | d median | d/c median | d/c min | d/c max | d/c>1 |")
print("|---|---:|---:|---:|---:|---:|---:|---:|")
by = collections.defaultdict(list)
for smp, sl, s, cc, dd, q, nr, kd in rows:
    by[s].append((cc, dd, q, smp, sl))
for s in sorted(by, key=lambda x: -len(by[x])):
    v = by[s]
    qs = sorted(x[2] for x in v)
    print("| %s | %d | %d | %d | %.3f | **%.3f** | %.3f | %d |"
          % (s, len(v), statistics.median([x[0] for x in v]), statistics.median([x[1] for x in v]),
             statistics.median(qs), qs[0], qs[-1], sum(1 for x in qs if x > 1.0)))

hi = by.get("hifiasm", [])
mn = min(hi, key=lambda x: x[2])
thr = mn[2]
print()
print("hifiasm haplotypes %d, min d/c = %.4f  (%s slot %s, c=%d, d=%d)"
      % (len(hi), thr, mn[3], mn[4], mn[0], mn[1]))
qs = sorted(x[2] for x in hi)
print("   hifiasm d/c distribution: q10 %.3f, q25 %.3f, median %.3f, q75 %.3f, max %.3f"
      % (qs[int(len(qs) * .1)], qs[int(len(qs) * .25)], statistics.median(qs), qs[int(len(qs) * .75)], qs[-1]))

print()
print("## reads to remove by threshold (all haplotypes)")
tot = sum(len(pos[k]) for k in pos)
print("| threshold d/c | reads to remove | fraction | haplotypes affected |")
print("|---|---:|---:|---:|")
for t2 in (thr, 0.20, 0.30, 0.50, 0.80, 1.00):
    rem = 0; aff = 0
    for k in pos:
        if k not in c or not c[k]:
            continue
        r2 = sum(1 for u in pos[k].values() if len(u) / float(c[k]) > t2)
        rem += r2
        if r2:
            aff += 1
    print("| %.4f%s | %s | %.2f%% | %d / %d |"
          % (t2, " 최소" if abs(t2 - thr) < 1e-9 else "", format(rem, ","), 100.0 * rem / tot, aff, len(pos)))
print("read total %s" % format(tot, ","))
