# -*- coding: utf-8 -*-
# 5S in other species: the canonical sequence (a flank-array-flank representative) plus a consensus map (a per-unit variant table), in the same form as for human.
#   The canonical sequence carries FL bp of flank on each side and is oriented along the array's gene direction (reverse-complemented if that is the minus strand).
#   Unit frame: spacer then gene, as in the human canonical unit, so unit i runs from the previous gene end + 1 to this gene end.
#            (The human canonical unit runs AluY .. NTS2 and ends with the gene at 2122-2240.)
#   Units are numbered from 1 along the array's gene direction. The human convention that position 1 is distal has no basis in other species and is not used.
#            Instead the chromosome coordinate is recorded with every unit, so which end it lies at can always be recovered.
#   Variants: the same method as for human, taken with minimap2 -cx asm20 --cs against each species' representative unit, asm20 so that a highly diverged array is still found, and written as position plus base, positionD n, or positionI n.
#   Arguments: species key, fasta, tblout, structure table, output folder, [flank bp = 20000].
#   Outputs: canon/<species>__<chrom>_<start>_<end>.fa, the canonical sequence
#          map/<species>__<chrom>_<start>_<end>.map.tsv (consensus map)
#          canon_unit/<species>__....unit.fa, that array's representative unit
#          MAP_SUMMARY.tsv, one row per array
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
SP, FA, TBL, STRUCT, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
FL = int(sys.argv[6]) if len(sys.argv) > 6 else 20000
for d in ("canon", "map", "canon_unit"):
    os.makedirs(os.path.join(OUT, d), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(COMP.get(x, "N") for x in reversed(s))


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
    if r[ix["species"]] != SP:
        continue
    ARR.append((r[ix["chrom"]], int(r[ix["array_start"]]), int(r[ix["array_end"]]), r[ix["strand_major"]][0], int(r[ix["copies"]])))


def faidx(reg):
    r = subprocess.run(["samtools", "faidx", FA, reg], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


def cs_to_vars(cs, rs):
    # minimap2 cs (short) -> the same tokens as the human map: 123A for a substitution, 123D2 for a deletion, 123I3 for an insertion, in 1-based representative-unit coordinates.
    out, p, i = [], rs, 0
    while i < len(cs):
        c = cs[i]
        if c == ":":
            j = i + 1
            while j < len(cs) and cs[j].isdigit():
                j += 1
            p += int(cs[i + 1:j])
            i = j
        elif c == "*":
            out.append("%d%s" % (p + 1, cs[i + 2].upper()))
            p += 1
            i += 3
        elif c == "-":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            n = j - i - 1
            out.append("%dD%d" % (p + 1, n))
            p += n
            i = j
        elif c == "+":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            out.append("%dI%d" % (p + 1, j - i - 1))
            i = j
        elif c == "~":
            j = i + 1
            while j < len(cs) and not (cs[j] in ":*-+~"):
                j += 1
            i = j
        else:
            i += 1
    return out


SUM = io.open(os.path.join(OUT, "MAP_SUMMARY_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
SUM.write(TAB.join(["species", "array_id", "chrom", "array_start", "array_end", "strand", "copies", "units_total",
                    "units_unaligned", "canon_unit_len", "units_identical_pct_of_aligned", "var_sites", "var_per_unit_med",
                    "canon_fa", "map_tsv"]) + NL)
for (c, s0, e1, strand, cop) in ARR:
    aid = "%s__%s_%d_%d" % (SP, c, s0, e1)
    L0, R1 = max(1, s0 - FL), min(LEN[c], e1 + FL)
    seq = faidx("%s:%d-%d" % (c, L0, R1))
    if not seq:
        continue
    hits = [h for h in HIT[c] if s0 <= h[0] and h[1] <= e1]
    # The working orientation is the array's gene direction; on the minus strand everything is flipped.
    if strand == "-":
        wseq = rc(seq)
        gen = sorted([(len(seq) - (h[1] - L0) - 1, len(seq) - (h[0] - L0) - 1) for h in hits])
    else:
        wseq = seq
        gen = sorted([(h[0] - L0, h[1] - L0) for h in hits])
    io.open(os.path.join(OUT, "canon", aid + ".fa"), "w", encoding="utf8", newline=NL).write(
        ">%s flank%d_array_flank%d strand=%s chrom=%s array=%d-%d%s%s%s" % (aid, FL, FL, strand, c, s0, e1, NL, wseq, NL))
    units, qf = [], io.open(os.path.join(OUT, "map", aid + ".units.fa"), "w", encoding="utf8", newline=NL)
    for i in range(1, len(gen)):
        us, ue = gen[i - 1][1] + 1, gen[i][1]
        if ue <= us:
            continue
        u = wseq[us:ue + 1]
        units.append((i, us, ue, u, gen[i][0] - us))
        qf.write(">u%d%s%s%s" % (i, NL, u, NL))
    qf.close()
    if not units:
        continue
    ln = collections.Counter(len(u[3]) for u in units).most_common(1)[0][0]
    same = [u[3] for u in units if len(u[3]) == ln]
    canon = "".join(collections.Counter(x[i] for x in same).most_common(1)[0][0] for i in range(ln)) if same else units[0][3]
    cp = os.path.join(OUT, "canon_unit", aid + ".unit.fa")
    io.open(cp, "w", encoding="utf8", newline=NL).write(">%s_canon_unit len=%d%s%s%s" % (aid, ln, NL, canon, NL))
    r = subprocess.run("minimap2 -cx asm20 --cs -t 4 %s %s 2>/dev/null" % (cp, os.path.join(OUT, "map", aid + ".units.fa")),
                       shell=True, capture_output=True, text=True)
    best = {}
    for L in r.stdout.splitlines():
        f2 = L.split(TAB)
        if len(f2) < 12:
            continue
        q, rs, bs = f2[0], int(f2[7]), int(f2[13].split(":")[-1]) if f2[13].startswith("AS") else 0
        cs = ""
        for x in f2[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        if q not in best or bs > best[q][0]:
            best[q] = (bs, rs, cs, int(f2[9]), int(f2[10]))
    M = io.open(os.path.join(OUT, "map", aid + ".map.tsv"), "w", encoding="utf8", newline=NL)
    M.write(TAB.join(["species", "array_id", "chrom", "strand", "array_pos", "chrom_start", "chrom_end", "unit_len",
                      "gene_start_in_unit", "gene_len", "aln_len", "identity_pct", "n_var", "varlist"]) + NL)
    nv, ident, allv, noaln = [], 0, collections.Counter(), 0
    for (i, us, ue, u, gs) in units:
        b = best.get("u%d" % i)
        # a unit that failed to align is marked NOALN, not counted as carrying no variants
        if b:
            v = cs_to_vars(b[2], b[1])
            vstr = ",".join(v) if v else "-"
            nvstr = str(len(v))
            idps = "%.2f" % (100.0 * b[3] / b[4]) if b[4] else "NA"
            alns = str(b[4])
            nv.append(len(v))
            ident += (len(v) == 0)
        else:
            v, vstr, nvstr, idps, alns = [], "NOALN", "NA", "NA", "0"
            noaln += 1
        if strand == "-":
            cs_, ce_ = R1 - ue, R1 - us
        else:
            cs_, ce_ = L0 + us, L0 + ue
        M.write(TAB.join([SP, aid, c, strand, str(i), str(cs_), str(ce_), str(len(u)), str(gs + 1),
                          str(gen[i][1] - gen[i][0] + 1), alns, idps, nvstr, vstr]) + NL)
        for x in v:
            allv[x] += 1
    M.close()
    nv.sort()
    na = len(units) - noaln
    SUM.write(TAB.join([SP, aid, c, str(s0), str(e1), strand, str(cop), str(len(units)), str(noaln), str(ln),
                        "%.1f" % (100.0 * ident / na) if na else "NA", str(len(allv)), str(nv[len(nv) // 2]) if nv else "NA",
                        os.path.join("canon", aid + ".fa"), os.path.join("map", aid + ".map.tsv")]) + NL)
    print("%s %s : units %d (alignment failed %d), representative unit %d bp, among aligned variants none %s, variant positions %d"
          % (SP, aid, len(units), noaln, ln, ("%.1f%%" % (100.0 * ident / na)) if na else "NA", len(allv)))
SUM.close()
