# -*- coding: utf-8 -*-
# In each assembly contig, find the 'distal end' by flank alignment,
# cut at unit boundary = 5S gene END, and extract the 20 distal-side units.
# Uses the same minimap2 parameters as asm_audit.py (numeric reproducibility).
import os, re, subprocess, sys, json

W = os.path.expanduser("~/p11_encode")
GQ = W + "/v5/5S_gene_query.fa"
FD = W + "/ref/FLANK_distal_50000.fa"
FP = W + "/ref/FLANK_proximal_50000.fa"
OUT = "/tmp/distal20"
os.makedirs(OUT, exist_ok=True)

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda x: x.translate(COMP)[::-1]

def rd(p):
    n, s = None, []
    for L in open(p):
        if L[0] == ">":
            if n: yield n, "".join(s)
            n = L[1:].split()[0]; s = []
        else: s.append(L.strip())
    if n: yield n, "".join(s)

def run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout

TARGETS = [
    ("K562_hapA",    W + "/k562ul2/K562ul2.bp.hap1.p_ctg.fa", "h1tg000001l"),
    ("K562_hapB",    W + "/k562ul2/K562ul2.bp.hap2.p_ctg.fa", "h2tg000001l"),
    ("GM12878_hap1", W + "/gm12878/units/hap1.fa",            "hap1"),
    ("GM12878_hap2", W + "/gm12878/units/hap2.fa",            "hap2"),
    ("HepG2_piece4", W + "/hepg2_pieces.fa",                  "h1tg000004l"),
    ("HepG2_piece5", W + "/hepg2_pieces.fa",                  "h1tg000005l"),
]

tsv = open(OUT + "/DISTAL20_UNITS.tsv", "w")
tsv.write("assembly\tcontig\tarray_pos\tcontig_start_1b\tcontig_end_1b\tcontig_strand\tlength\tgene_end_boundary_distal\tgene_end_boundary_prox\n")
log = open(OUT + "/DISTAL20_LOG.txt", "w")
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s); log.write(s + "\n")

TMP = OUT + "/_one.fa"
for label, path, ctg in TARGETS:
    seqs = dict(rd(path))
    if ctg not in seqs:
        P("!! %s : contig %s 없음 (%s)" % (label, ctg, path)); continue
    s = seqs[ctg]
    open(TMP, "w").write(">%s\n%s\n" % (ctg, s))
    P("\n" + "=" * 90)
    P("### %s  %s  %s  len=%d  N=%d" % (label, os.path.basename(path), ctg, len(s), s.count("N")))

    # --- flank alignment (same parameters as asm_audit)
    fl = {}
    for side, ref in (("prox", FP), ("dist", FD)):
        rows = []
        for L in run("minimap2 -c -x asm20 -t 8 %s %s 2>/dev/null" % (ref, TMP)).splitlines():
            f = L.split("\t")
            rows.append((int(f[10]), int(f[9]) / int(f[10]), int(f[2]), int(f[3])))
        rows.sort(reverse=True)
        fl[side] = rows[0] if rows and rows[0][0] >= 5000 else None
        r = fl[side]
        P("  %s flank : %s" % (side, ("block %d id %.4f contig %d-%d" % r if r else "없음")))

    # --- gene hits (same parameters as asm_audit) : start,end,strand
    hits = []
    for L in run("minimap2 -c -x map-ont -N 900 -p 0.2 -t 8 %s %s 2>/dev/null" % (GQ, TMP)).splitlines():
        f = L.split("\t")
        if int(f[10]) >= 100 and int(f[9]) / int(f[10]) >= 0.85:
            hits.append((int(f[2]), int(f[3]), f[4]))
    hits.sort()
    ded = []
    for h in hits:
        if not ded or h[0] - ded[-1][0] > 800: ded.append(h)
    P("  유전자 히트 %d 개 · 구간 %d-%d · strand %s" %
      (len(ded), ded[0][0], ded[-1][0],
       dict((x, sum(1 for h in ded if h[2] == x)) for x in set(h[2] for h in ded))))
    gaps = [b[0] - a[0] for a, b in zip(ded, ded[1:])]
    if gaps:
        import statistics
        inr = [g for g in gaps if 1500 < g < 3500]
        P("  간격 중앙 %s · 1500-3500 밖 %d %s" % (int(statistics.median(inr)) if inr else "-",
                                                len(gaps) - len(inr), sorted(g for g in gaps if not (1500 < g < 3500))[:8]))

    # --- orientation : DESC if the distal flank has higher coordinates than the gene cluster
    d = fl["dist"]
    if d is None:
        P("  !! 원위 flank 없음 → 건너뜀"); continue
    gmid = (ded[0][0] + ded[-1][0]) / 2
    desc = d[2] > gmid          # True: distal = higher coordinates → units run in decreasing coordinate direction
    P("  방향 : 원위 flank %s → %s" % ("HIGHER" if desc else "LOWER",
                                      "contig 좌표 감소" if desc else "contig 좌표 증가"))

    # --- unit boundary = gene END (contig coordinate of canonical base 2240)
    # desc (= stored in reverse): the canonical end is the left end of the hit, f[2]
    # asc (= stored forward): the canonical end is the right end of the hit, f[3]
    if desc:
        B = sorted([h[0] for h in ded], reverse=True)   # 0-based, coordinate just right of the boundary
        flank_edge = d[2]                                # array-side end of the distal flank
    else:
        B = sorted([h[1] for h in ded])                  # 0-based exclusive end
        flank_edge = d[3]
    P("  원위 flank 끝 → 첫 유전자 END 거리 = %d bp" % abs(B[0] - flank_edge))

    if len(B) < 21:
        P("  !! 경계 %d 개뿐 → 20 유닛 불가" % len(B))
    nmax = min(20, len(B) - 1)
    fa = open(OUT + "/%s_distal20.fa" % label, "w")
    for k in range(nmax):
        b1, b2 = B[k], B[k + 1]        # b1 = distal-side boundary, b2 = proximal-side boundary
        if desc:
            sub = s[b2:b1]; sub = rc(sub); cs, ce, st = b2 + 1, b1, "-"
        else:
            sub = s[b1:b2];              cs, ce, st = b1 + 1, b2, "+"
        ap = k + 1
        fa.write(">%s_pos%03d %s:%d-%d(%s) len=%d\n%s\n" % (label, ap, ctg, cs, ce, st, len(sub), sub))
        tsv.write("%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n" % (label, ctg, ap, cs, ce, st, len(sub), b1, b2))
    fa.close()
    P("  → %s_distal20.fa : %d 유닛 · 길이 %s" %
      (label, nmax, [ (B[k]-B[k+1] if desc else B[k+1]-B[k]) for k in range(min(nmax,20)) ]))

tsv.close(); log.close()
print("\nDONE")
