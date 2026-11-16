# -*- coding: utf-8 -*-
# 5S in other species: build a reference (representative) unit per array, with the human G8v4 method carried across per species.
#   Frame: a unit runs from gene end + 1 to the next gene end, so it ends as [spacer][gene], the same phase as the human canonical unit.
#   Gene boundaries: blastn with that array's gene consensus, then extend to the full length using the query coordinates.
#                 (plus strand: start = sstart - (qstart - 1), end = send + (L - qend); the minus strand is the reverse)
#                 Measured in human: if the last base of the gene carries a variant, the HSP stops one or two bases short.
#   Per species: PAD = 5% (at least 20), CTX = 9% (at least 40), LO/HI = -6% / +7% of the period (at least +-30 bp).
#   Representative unit: a per-position majority over the units of modal length, then a check from the sequence that the gene lies at the end of the unit.
#   Outputs: refunit/<array_id>.unit.fa, refunit/<array_id>.gene.fa, REFUNIT_PARAMS.tsv.
#   Arguments: species key, fasta, tblout, structure table, output folder.
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
SP, FA, TBL, STRUCT, OUT = sys.argv[1:6]
os.makedirs(os.path.join(OUT, "refunit"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(COMP.get(x, "N") for x in reversed(s))


def faidx(reg):
    r = subprocess.run(["samtools", "faidx", FA, reg], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


LEN = {}
for L in io.open(FA + ".fai", encoding="utf8"):
    f = L.split(TAB)
    LEN[f[0]] = int(f[1])
HIT = collections.defaultdict(list)
for L in io.open(TBL, encoding="utf8"):
    if L.startswith("#"):
        continue
    f = L.split()
    a, b = int(f[7]), int(f[8])
    HIT[f[0]].append((min(a, b), max(a, b), f[9]))
for c in HIT:
    HIT[c].sort()
ARR = []
f = io.open(STRUCT, encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ix = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["species"]] == SP:
        ARR.append((r[ix["chrom"]], int(r[ix["array_start"]]), int(r[ix["array_end"]]), r[ix["strand_major"]][0]))
P = io.open(os.path.join(OUT, "REFUNIT_PARAMS_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
P.write(TAB.join(["species", "array_id", "chrom", "array_start", "array_end", "strand", "n_gene_hsp", "period_med",
                  "LO", "HI", "PAD", "CTX", "n_units_in_range", "n_break", "refunit_len", "gene_len",
                  "gene_start_in_refunit", "gene_at_end", "gene_identity_pct", "units_ge99_pct_blastn",
                  "units_median_identity", "refunit_usable", "note"]) + NL)
for (c, s0, e1, strand) in ARR:
    aid = "%s__%s_%d_%d" % (SP, c, s0, e1)
    hits = [h for h in HIT[c] if s0 <= h[0] and h[1] <= e1]
    if len(hits) < 3:
        continue
    mar = 3000
    L0, R1 = max(1, s0 - mar), min(LEN[c], e1 + mar)
    seq = faidx("%s:%d-%d" % (c, L0, R1))
    if not seq:
        continue
    # 1) the gene consensus of that array: a per-position majority over the cmsearch hit sequences, using the modal length only
    gs = []
    for (a, b, st) in hits:
        g = seq[a - L0:b - L0 + 1]
        gs.append(rc(g) if st == "-" else g)
    gl = collections.Counter(len(x) for x in gs).most_common(1)[0][0]
    same = [x for x in gs if len(x) == gl]
    gene = "".join(collections.Counter(x[i] for x in same).most_common(1)[0][0] for i in range(gl))
    gp = os.path.join(OUT, "refunit", aid + ".gene.fa")
    io.open(gp, "w", encoding="utf8", newline=NL).write(">%s_gene len=%d%s%s%s" % (aid, gl, NL, gene, NL))
    # 2) gene HSPs by blastn, extended to the full length using the query coordinates
    rp = os.path.join(OUT, "refunit", aid + ".region.fa")
    io.open(rp, "w", encoding="utf8", newline=NL).write(">region%s%s%s" % (NL, seq, NL))
    rb = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1e-5 -query %s -subject %s "
                        "-outfmt '6 qstart qend sstart send pident bitscore' 2>/dev/null" % (gp, rp),
                        shell=True, capture_output=True, text=True)
    G = []
    for L in rb.stdout.splitlines():
        q1, q2, s1, s2, pid, bs = L.split(TAB)
        q1, q2, s1, s2 = int(q1), int(q2), int(s1), int(s2)
        if s1 <= s2:
            st, a, b = "+", s1 - (q1 - 1), s2 + (gl - q2)
        else:
            st, a, b = "-", s2 - (gl - q2), s1 + (q1 - 1)
        G.append((a, b, st, float(bs)))
    G.sort()
    keep, last = [], -1
    for (a, b, st, bs) in G:
        if a > last:
            keep.append((a, b, st))
            last = b
    major = collections.Counter(x[2] for x in keep).most_common(1)[0][0] if keep else strand
    if major == "-":
        n = len(seq)
        wseq = rc(seq)
        gg = sorted([(n - b + 1, n - a + 1) for (a, b, st) in keep if st == major])
    else:
        wseq = seq
        gg = sorted([(a, b) for (a, b, st) in keep if st == major])
    if len(gg) < 3:
        P.write(TAB.join([SP, aid, c, str(s0), str(e1), strand, str(len(keep)), "NA", "NA", "NA", "NA", "NA",
                          "0", "NA", "NA", str(gl), "NA", "NA", "NA", "NA", "NA", "N", "유전자 HSP 부족"]) + NL)
        continue
    per = [gg[i + 1][1] - gg[i][1] for i in range(len(gg) - 1)]
    per.sort()
    pm = per[len(per) // 2]
    # The human values (LO 2100 and HI 2400 for a period of 2240, that is -6.25% / +7.1%) are carried across per species, taking whichever bound is the more permissive.
    LO = min(pm - 30, int(round(pm * 0.94)))
    HI = max(pm + 30, int(round(pm * 1.07)))
    PAD = max(20, int(round(pm * 0.05)))
    CTX = max(40, int(round(pm * 0.09)))
    units, brk = [], 0
    for i in range(len(gg) - 1):
        a, b = gg[i][1] + 1, gg[i + 1][1]
        n = b - a + 1
        if LO <= n <= HI:
            units.append(wseq[a - 1:b])
        else:
            brk += 1
    if not units:
        P.write(TAB.join([SP, aid, c, str(s0), str(e1), major, str(len(gg)), str(pm), str(LO), str(HI), str(PAD),
                          str(CTX), "0", str(brk), "NA", str(gl), "NA", "NA", "NA", "NA", "NA", "N", "범위 안 유닛 없음"]) + NL)
        continue
    ml = collections.Counter(len(u) for u in units).most_common(1)[0][0]
    pool = [u for u in units if len(u) == ml]
    ref = "".join(collections.Counter(u[i] for u in pool).most_common(1)[0][0] for i in range(ml))
    up = os.path.join(OUT, "refunit", aid + ".unit.fa")
    io.open(up, "w", encoding="utf8", newline=NL).write(">%s_refunit len=%d gene_len=%d%s%s%s" % (aid, ml, gl, NL, ref, NL))
    # 3) checks: does the gene lie at the end of the unit (from the sequence), and how closely do the units match the representative
    tail = ref[-gl:]
    idp = 100.0 * sum(1 for i in range(gl) if tail[i] == gene[i]) / gl
    rb2 = subprocess.run("blastn -task blastn -word_size 11 -dust no -query %s -subject %s "
                         "-outfmt '6 sstart send pident bitscore' 2>/dev/null" % (gp, up),
                         shell=True, capture_output=True, text=True)
    gstart = "NA"
    bestbs = -1
    for L in rb2.stdout.splitlines():
        s1, s2, pid, bs = L.split(TAB)
        if float(bs) > bestbs:
            bestbs = float(bs)
            gstart = min(int(s1), int(s2))
    atend = "Y" if gstart != "NA" and abs((ml - gl + 1) - gstart) <= 2 else "N"
    # Homogeneity is measured by alignment (blastn); a direct position-by-position comparison breaks on a single indel.
    ufp = os.path.join(OUT, "refunit", aid + ".units.fa")
    with io.open(ufp, "w", encoding="utf8", newline=NL) as g2:
        for k, u in enumerate(units):
            g2.write(">u%d%s%s%s" % (k + 1, NL, u, NL))
    rb3 = subprocess.run("blastn -task blastn -word_size 11 -dust no -query %s -subject %s "
                         "-outfmt '6 qseqid pident length qlen bitscore' 2>/dev/null" % (ufp, up),
                         shell=True, capture_output=True, text=True)
    bb = {}
    for L in rb3.stdout.splitlines():
        q, pid, aln, qlen, bs = L.split(TAB)
        if int(aln) < 0.5 * int(qlen):
            continue
        if q not in bb or float(bs) > bb[q][0]:
            bb[q] = (float(bs), float(pid))
    pids = sorted(v[1] for v in bb.values())
    medid = "%.2f" % pids[len(pids) // 2] if pids else "NA"
    ge = sum(1 for v in bb.values() if v[1] >= 99.0)
    gepct = "%.1f" % (100.0 * ge / len(units)) if units else "NA"
    try:
        os.remove(ufp)
    except OSError:
        pass
    # Is the representative usable: at least 5 units, the gene at the end, and at least half the units aligning to the representative.
    why = []
    if len(units) < 5:
        why.append("유닛<5")
    if atend != "Y":
        why.append("유전자가 끝 아님")
    if len(bb) < 0.5 * len(units):
        why.append("정렬된 유닛<50%")
    usable = "Y" if not why else "N"
    P.write(TAB.join([SP, aid, c, str(s0), str(e1), major, str(len(gg)), str(pm), str(LO), str(HI), str(PAD), str(CTX),
                      str(len(units)), str(brk), str(ml), str(gl), str(gstart), atend, "%.1f" % idp,
                      gepct, medid, usable, ";".join(why)]) + NL)
    print("%s : 유전자HSP %d · 주기 %d (LO %d HI %d) · 유닛 %d (끊김 %d) · 대표 %d bp · 유전자 끝 %s (자리 %s · 일치 %.1f%%)"
          % (aid, len(gg), pm, LO, HI, len(units), brk, ml, atend, gstart, idp))
    for x in (rp,):
        try:
            os.remove(x)
        except OSError:
            pass
P.close()
