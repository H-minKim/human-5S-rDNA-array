# -*- coding: utf-8 -*-
# Arrays are grouped into types by comparing spacers only; the 119 bp gene is shared across all species and must not be used for grouping.
#   A unit ends as spacer then gene, so the spacer is unit[:-gene_length].
#   Same type means blastn identity >= ID and the alignment covering at least COV of the longer spacer (defaults 80% and 0.50).
#   The distribution of gaps between gene ends (the period) is also reported per array, to show whether one unit or several were taken as one.
#   Outputs: SPACER_TYPES.tsv, PERIOD_DIST.tsv.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import glob
import subprocess
import collections
import statistics

NL = chr(10)
TAB = chr(9)
X = os.path.expanduser("~/xsp")
ID = float(sys.argv[1]) if len(sys.argv) > 1 else 80.0
COV = float(sys.argv[2]) if len(sys.argv) > 2 else 0.50
OUT = os.path.join(X, "unit_types")
os.makedirs(OUT, exist_ok=True)
ST, ix = {}, {}
f = io.open(os.path.join(X, "struct", "ALL.struct.tsv"), encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ix = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    ST["%s__%s_%s_%s" % (r[ix["species"]], r[ix["chrom"]], r[ix["array_start"]], r[ix["array_end"]])] = r
UNIT, SPACER = {}, {}
for p in sorted(glob.glob(os.path.join(X, "smap2", "canon_unit", "*.unit.fa"))):
    aid = os.path.basename(p).replace(".unit.fa", "")
    s = "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()
    UNIT[aid] = s
    gl = int(ST[aid][ix["gene_len_med"]]) if aid in ST else 119
    SPACER[aid] = s[:-gl] if len(s) > gl + 20 else ""
BYSP = collections.defaultdict(list)
for a in UNIT:
    BYSP[a.split("__")[0]].append(a)
rep = io.open(os.path.join(OUT, "SPACER_TYPES.tsv"), "w", encoding="utf8", newline=NL)
rep.write(TAB.join(["species", "array_id", "copies", "period_med", "unit_len", "spacer_len", "gc_spacer",
                    "spacer_type", "type_arrays", "type_copies", "pident_to_type_rep", "aln_frac_of_longer"]) + NL)
for sp in sorted(BYSP):
    aids = sorted([a for a in BYSP[sp] if SPACER.get(a)], key=lambda a: -int(ST[a][ix["copies"]]))
    if not aids:
        continue
    fa = os.path.join(OUT, sp + ".spacers.fa")
    with io.open(fa, "w", encoding="utf8", newline=NL) as g:
        for a in aids:
            g.write(">" + a + NL + SPACER[a] + NL)
    best = {}
    if len(aids) > 1:
        r = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1e-3 -query %s -subject %s "
                           "-outfmt '6 qseqid sseqid pident length bitscore' 2>/dev/null" % (fa, fa),
                           shell=True, capture_output=True, text=True)
        for L in r.stdout.splitlines():
            q, s, pid, aln, bs = L.split(TAB)
            if q == s:
                continue
            v = (float(bs), float(pid), int(aln))
            if (q, s) not in best or v[0] > best[(q, s)][0]:
                best[(q, s)] = v
    adj = collections.defaultdict(set)
    for (q, s), (bs, pid, aln) in best.items():
        mx = max(len(SPACER[q]), len(SPACER[s]))
        if pid >= ID and aln >= COV * mx:
            adj[q].add(s)
            adj[s].add(q)
    seen, types = set(), []
    for a in aids:
        if a in seen:
            continue
        st, comp = [a], []
        seen.add(a)
        while st:
            x = st.pop()
            comp.append(x)
            for y in adj.get(x, ()):
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        types.append(sorted(comp, key=lambda z: -int(ST[z][ix["copies"]])))
    types.sort(key=lambda t: -sum(int(ST[z][ix["copies"]]) for z in t))
    print("== %s : array %d, spacer types %d" % (sp, len(aids), len(types)))
    for i, t in enumerate(types, 1):
        cps = sum(int(ST[z][ix["copies"]]) for z in t)
        sl = sorted(len(SPACER[z]) for z in t)
        gc = [float(ST[z][ix["gc_spacer"]]) for z in t if ST[z][ix["gc_spacer"]] != "NA"]
        print("   type %d : array %2d, copy %5d, spacer %s bp, GC %s" % (
            i, len(t), cps, ("%d-%d" % (sl[0], sl[-1])) if sl[0] != sl[-1] else str(sl[0]),
            ("%.1f-%.1f" % (min(gc), max(gc))) if gc else "NA"))
        for z in t[:3]:
            print("        %s copy %s" % (z.split("__")[1], ST[z][ix["copies"]]))
        if len(t) > 3:
            print("        ... %d more" % (len(t) - 3))
        r0 = t[0]
        for z in t:
            b = best.get((z, r0))
            pid = "100.00" if z == r0 else ("%.2f" % b[1] if b else "NA")
            af = "1.000" if z == r0 else ("%.3f" % (b[2] / float(max(len(SPACER[z]), len(SPACER[r0])))) if b else "NA")
            rr = ST[z]
            rep.write(TAB.join([sp, z, rr[ix["copies"]], rr[ix["period_med"]], str(len(UNIT[z])),
                                str(len(SPACER[z])), rr[ix["gc_spacer"]], str(i), str(len(t)), str(cps), pid, af]) + NL)
rep.close()
# The period distribution: was one unit taken, or several as one?
pd = io.open(os.path.join(OUT, "PERIOD_DIST.tsv"), "w", encoding="utf8", newline=NL)
pd.write(TAB.join(["array_id", "n_gap", "p10", "median", "p90", "min", "max", "modal_bin", "bins"]) + NL)
HIT = collections.defaultdict(list)
for sp in sorted(BYSP):
    tbl = os.path.join(X, "cms", sp + ".tbl")
    if not os.path.exists(tbl):
        continue
    H = collections.defaultdict(list)
    for L in io.open(tbl, encoding="utf8"):
        if L.startswith("#"):
            continue
        ff = L.split()
        a, b = int(ff[7]), int(ff[8])
        H[ff[0]].append((min(a, b), max(a, b)))
    for c in H:
        H[c].sort()
    for aid in sorted(BYSP[sp], key=lambda a: -int(ST[a][ix["copies"]])):
        r = ST[aid]
        c, s0, e1 = r[ix["chrom"]], int(r[ix["array_start"]]), int(r[ix["array_end"]])
        hs = [h for h in H[c] if s0 <= h[0] and h[1] <= e1]
        ends = sorted(h[1] for h in hs)
        gaps = [ends[i + 1] - ends[i] for i in range(len(ends) - 1)]
        if not gaps:
            continue
        gs = sorted(gaps)
        bins = collections.Counter(int(g / 100) * 100 for g in gaps)
        top = bins.most_common(4)
        pd.write(TAB.join([aid, str(len(gaps)), str(gs[int(0.1 * len(gs))]), str(int(statistics.median(gs))),
                           str(gs[int(0.9 * (len(gs) - 1))]), str(gs[0]), str(gs[-1]), str(top[0][0]),
                           ",".join("%d:%d" % (k, v) for k, v in top)]) + NL)
pd.close()
print("tables : SPACER_TYPES.tsv, PERIOD_DIST.tsv")
