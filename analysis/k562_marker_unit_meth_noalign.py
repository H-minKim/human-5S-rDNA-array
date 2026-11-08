# -*- coding: utf-8 -*-
# . K562: are the marker variants in units whose methylation state is the expected one? Measured inside the reads, with no assembly and no placement (read only, new folder).
#   The question this script answers: looking at the K562 methylation reads without alignment, do all the marker variants fall in units with the methylation we expect? Look unit by unit, and among those take the ones that carry a marker variant. Is it close to what we expect?
#   ① Reads (K562r.hit.fq with 5mC MM/ML) are aligned to the canonical unit only (unit_canon, 2240 bp) and cut into units (minimap2, query = read, target = canonical unit, aligned canonical unit length >= 1800).
#   ② For each unit, the base at every marker position (canonical unit coordinates, 70 K562 markers) gives the markers that unit carries.
#   ③ For each unit, gene body (2122 to 2240) CpG methylation: ML >= 128 counts as methylated, and only units with at least 3 calls are used.
#   ④ Expectation: proximal and distal markers sit in low methylation units, interior markers in high methylation units. The marker zones were defined on the assembly, but this count does not use the assembly.
import os, re, sys, csv, subprocess, collections, array
import pysam
P = os.path.expanduser("~/p11_v11"); W = os.path.join(P, "K562MK_NOALIGN"); os.makedirs(W, exist_ok=True)
FQ = os.environ.get("FQ", os.path.join(P, "AO_K562r/fq/K562r.hit.fq")); UNIT = os.path.expanduser("~/s5s_tuc/ref/unit_canon.fa")
MK = os.path.join(W, "P11_K562_MARKERS_NEWZONE_2026-10-13.tsv")       # the marker table uploaded from the local machine
G0, G1, MINCALL, MLT = 2122, 2240, 3, 128
COMP = str.maketrans("ACGTN", "TGCAN")
markers = [r for r in csv.DictReader(open(MK), delimiter="\t")]
ZN = {"근위(합)": "proximal", "가운데": "interior", "원위(합)": "distal"}
MKP = collections.defaultdict(list)
for r in markers: MKP[int(r["pos"])].append((r["variant"], r["alt"], ZN[r["zone"]]))
print("표지", len(markers), collections.Counter(ZN[r["zone"]] for r in markers))
R = {}
with open(FQ) as g:
    while True:
        h = g.readline()
        if not h: break
        s = g.readline().rstrip(); g.readline(); q = g.readline().rstrip()
        f = h[1:].rstrip("\n").split(); rid = f[0]
        tags = {x[:2]: x[5:] for x in f[1:] if len(x) > 5 and x[2] == ":"}
        R[rid] = (s, q, tags)
with open(os.path.join(W, "reads.fa"), "w") as o:
    for k, v in R.items(): o.write(">%s\n%s\n" % (k, v[0]))
out = subprocess.run("minimap2 -x map-hifi -c -N 200 -p 0.05 --secondary=yes %s %s 2>/dev/null" % (UNIT, os.path.join(W, "reads.fa")), shell=True, capture_output=True, text=True).stdout
HITS = collections.defaultdict(list)
for l in out.splitlines():
    f = l.split("\t")
    qn, L, qs, qe, st, ts, te = f[0], int(f[1]), int(f[2]), int(f[3]), f[4], int(f[7]), int(f[8])
    if te - ts < 1800: continue
    cg = [x[5:] for x in f[12:] if x.startswith("cg:Z:")][0]
    HITS[qn].append((qs, qe, st, ts, te, L, cg))
print("read", len(R), "· 유닛이 맞은 read", len(HITS), "· 유닛", sum(len(v) for v in HITS.values()))
rows = []
for rid, hits in HITS.items():
    s, q, tags = R[rid]
    mm, ml = tags.get("MM") or tags.get("Mm"), tags.get("ML") or tags.get("Ml")
    cm = {}
    if mm and ml:
        a = pysam.AlignedSegment(); a.query_name = rid; a.query_sequence = s; a.flag = 4
        a.set_tag("MM", mm, value_type="Z"); a.set_tag("ML", array.array("B", [int(x) for x in (ml.split(",")[1:] if ml.startswith("C") else ml.split(","))]))
        for k, v in (a.modified_bases or {}).items():
            if k[0] == "C" and k[2] == "m": cm.update({p: x for p, x in v})
    for qs, qe, st, ts, te, L, cg in sorted(hits):
        t = ts; r = qs if st == "+" else L - qe                     # r = the read coordinate in the aligned orientation, that is, the canonical unit orientation
        t2i = {}
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg):
            n = int(n)
            if op in "M=X":
                for k in range(n): t2i[t + k] = (r + k) if st == "+" else (L - 1 - (r + k))
                t += n; r += n
            elif op == "I": r += n
            elif op in "DN": t += n
        def base(tp):                                                # the read base at a canonical unit coordinate (0-based), reported in the canonical unit orientation
            i = t2i.get(tp)
            if i is None: return "-"
            b = s[i]; return b if st == "+" else b.translate(COMP)
        carried = []
        for p, lst in MKP.items():
            if ts <= p - 1 < te:
                b = base(p - 1)
                for v, alt, z in lst:
                    if b == alt: carried.append((v, z))
        win = sorted({t2i[tp] for tp in range(G0 - 1, G1) if tp in t2i})
        lo, hi = (win[0], win[-1]) if win else (0, -1)
        cpg = [i for i in range(lo, hi + 1) if s[i:i + 2] == "CG"]   # the C on the read forward strand, because MM is given on the read forward strand
        calls = [cm[i] for i in cpg if i in cm]
        n = len(calls); m = sum(1 for x in calls if x >= MLT)
        rows.append((rid, qs, qe, st, n, m, ";".join("%s:%s" % c for c in carried)))
with open(os.path.join(W, "unit_marker_meth.tsv"), "w", newline="\n") as o:
    o.write("read\tq_start\tq_end\tstrand\tn_call\tn_meth\tmarkers\n")
    for x in rows: o.write("\t".join(map(str, x)) + "\n")
print("-> unit_marker_meth.tsv · 유닛", len(rows), "· call ≥3", sum(1 for x in rows if x[4] >= MINCALL))
