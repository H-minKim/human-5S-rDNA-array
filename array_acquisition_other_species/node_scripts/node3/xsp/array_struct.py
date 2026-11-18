# -*- coding: utf-8 -*-
# Measuring the 5S array structure per species: whether it is flank-array-flank, and the internal composition of the unit.
#   An array groups cmsearch hits within 10 kb and needs at least 3 copies, the same rule as the 5S array catalogue.
#   A unit is the interval between neighbouring hit starts, anchored on the gene.
#   Gene and spacer are the hit interval and the remainder; GC and CpG density are measured separately for each.
#   The flank is FL bp on each side of the array; the 5S hit count, the N count and the distance to the end of the contig are recorded.
#   Outputs: <species>.struct.tsv (one row per array), <species>.units.fa (up to 30 interior units) and <species>.flanks.fa (both sides).
#   Arguments: species name, fasta, tblout, output folder, [flank bp = 20000].
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import subprocess
import collections
import statistics

NL = chr(10)
TAB = chr(9)
SP, FA, TBL, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
FL = int(sys.argv[5]) if len(sys.argv) > 5 else 20000
GAP, MINC = 10000, 3
os.makedirs(OUT, exist_ok=True)
LEN = {}
for L in io.open(FA + ".fai", encoding="utf8"):
    f = L.split(TAB)
    LEN[f[0]] = int(f[1])
H = collections.defaultdict(list)
for L in io.open(TBL, encoding="utf8"):
    if L.startswith("#"):
        continue
    f = L.split()
    a, b = int(f[7]), int(f[8])
    H[f[0]].append((min(a, b), max(a, b), f[9], float(f[14])))
for c in H:
    H[c].sort()
ARR = []
for c, v in H.items():
    cur = [v[0]]
    for h in v[1:]:
        if h[0] - cur[-1][1] <= GAP:
            cur.append(h)
        else:
            if len(cur) >= MINC:
                ARR.append((c, cur))
            cur = [h]
    if len(cur) >= MINC:
        ARR.append((c, cur))
ARR.sort(key=lambda x: -len(x[1]))


def faidx(regs):
    if not regs:
        return {}
    p = os.path.join(OUT, "_r.txt")
    io.open(p, "w", newline=NL).write(NL.join(regs) + NL)
    r = subprocess.run(["samtools", "faidx", "-r", p, FA], capture_output=True, text=True)
    d, k = {}, None
    for L in r.stdout.splitlines():
        if L.startswith(">"):
            k = L[1:].strip()
            d[k] = []
        elif k:
            d[k].append(L.strip())
    return dict((x, "".join(y).upper()) for x, y in d.items())


def gc_cpg(s):
    s = s.replace("N", "")
    if not s:
        return "NA", "NA"
    gc = 100.0 * (s.count("G") + s.count("C")) / len(s)
    cpg = 1000.0 * s.count("CG") / len(s)
    return "%.1f" % gc, "%.1f" % cpg


O = io.open(os.path.join(OUT, SP + ".struct.tsv"), "w", encoding="utf8", newline=NL)
O.write(TAB.join(["species", "chrom", "array_start", "array_end", "span_bp", "copies", "strand_major",
                  "period_med", "period_sd", "gene_len_med", "spacer_len_med", "gc_gene", "gc_spacer",
                  "cpg_per_kb_gene", "cpg_per_kb_spacer", "unit_identity_pct",
                  "flank_L_hits", "flank_R_hits", "flank_L_N", "flank_R_N", "flank_L_gc", "flank_R_gc",
                  "dist_contig_start", "dist_contig_end", "contig_len", "verdict"]) + NL)
UF = io.open(os.path.join(OUT, SP + ".units.fa"), "w", encoding="utf8", newline=NL)
FF = io.open(os.path.join(OUT, SP + ".flanks.fa"), "w", encoding="utf8", newline=NL)
for ai, (c, hits) in enumerate(ARR):
    s0, e1 = hits[0][0], hits[-1][1]
    st = [h[0] for h in hits]
    strand = collections.Counter(h[2] for h in hits).most_common(1)[0][0]
    glen = [h[1] - h[0] + 1 for h in hits]
    per = [st[i + 1] - st[i] for i in range(len(st) - 1)]
    if not per:
        continue
    pm = statistics.median(per)
    psd = statistics.pstdev(per) if len(per) > 1 else 0.0
    mid = len(per) // 2
    lo = max(0, mid - 15)
    idx = list(range(lo, min(len(per), lo + 30)))
    regs = ["%s:%d-%d" % (c, st[i], st[i + 1] - 1) for i in idx]
    seqs = faidx(regs)
    units = [(regs[j], seqs.get(regs[j], "")) for j in range(len(regs))]
    gseq, sseq = [], []
    for j, i in enumerate(idx):
        u = units[j][1]
        if not u:
            continue
        UF.write(">%s_arr%d_u%d%s%s%s" % (SP, ai, i, NL, u, NL))
        gl = hits[i][1] - hits[i][0] + 1
        gseq.append(u[:gl])
        sseq.append(u[gl:])
    gg, gc_ = gc_cpg("".join(gseq))
    sg, sc_ = gc_cpg("".join(sseq))
    # How similar the units are: one interior unit is taken as the reference and the others compared to it position by position.
    ident = "NA"
    us = [u for _n, u in units if u]
    if len(us) >= 3:
        ref = us[len(us) // 2]
        sc = []
        for u in us:
            if u is ref or not u:
                continue
            n = min(len(u), len(ref))
            if n < 50:
                continue
            m = sum(1 for k in range(n) if u[k] == ref[k])
            sc.append(100.0 * m / n)
        if sc:
            ident = "%.1f" % statistics.median(sc)
    Ls, Le = max(1, s0 - FL), s0 - 1
    Rs, Re = e1 + 1, min(LEN[c], e1 + FL)
    fr = []
    if Le > Ls:
        fr.append("%s:%d-%d" % (c, Ls, Le))
    if Re > Rs:
        fr.append("%s:%d-%d" % (c, Rs, Re))
    fs = faidx(fr)
    lh = sum(1 for h in H[c] if Ls <= h[0] <= Le)
    rh = sum(1 for h in H[c] if Rs <= h[0] <= Re)
    lseq = fs.get("%s:%d-%d" % (c, Ls, Le), "")
    rseq = fs.get("%s:%d-%d" % (c, Rs, Re), "")
    if lseq:
        FF.write(">%s_arr%d_L%s%s%s" % (SP, ai, NL, lseq, NL))
    if rseq:
        FF.write(">%s_arr%d_R%s%s%s" % (SP, ai, NL, rseq, NL))
    lgc = gc_cpg(lseq)[0]
    rgc = gc_cpg(rseq)[0]
    ln, rn = lseq.count("N"), rseq.count("N")
    dstart, dend = s0 - 1, LEN[c] - e1
    ok = (len(lseq) >= FL * 0.9 and len(rseq) >= FL * 0.9 and ln == 0 and rn == 0 and lh == 0 and rh == 0)
    verdict = "flank-array-flank" if ok else "확인필요"
    O.write(TAB.join([SP, c, str(s0), str(e1), str(e1 - s0 + 1), str(len(hits)), strand,
                      "%.0f" % pm, "%.1f" % psd, "%.0f" % statistics.median(glen),
                      "%.0f" % (pm - statistics.median(glen)), gg, sg, gc_, sc_, ident,
                      str(lh), str(rh), str(ln), str(rn), lgc, rgc,
                      str(dstart), str(dend), str(LEN[c]), verdict]) + NL)
O.close()
UF.close()
FF.close()
try:
    os.remove(os.path.join(OUT, "_r.txt"))
except OSError:
    pass
print(SP, "array", len(ARR))
