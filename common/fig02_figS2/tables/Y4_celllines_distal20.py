# -*- coding: utf-8 -*-
r"""
Y4 - is the '2-3 position spacing' rule reproduced in the distal 20 units of cell lines (K562, GM12878, HepG2)

Data (read only): `$RDNA_DATA_ROOT/P11_ENCODE_unit_placement/04_assembly_handover/*_distal20.fa` (, handed over 2026-09-24)
  - array_pos 1 = distal (DJ) end, sequences already in canonical orientation, unit boundary = 5S gene END
Method (no aligner)
  - unit similarity = Jaccard of 31-mer sets (a flat comparison without variant calling)
  - terminal-type marker = context of the 1 bp deletion at position 1096 + context of 635 C>T (context ±8 bp, same as for HiFi)
    At 1096 the canonical base is T, but 92% of real units are T>G -> the common allele is set to G
  - measured: terminal type by depth, spacing between terminal types, similarity by spacing (spacing 1 vs spacing 2), number of groups within a haplotype
Output: $OUTDIR\Y4_cellline_units.tsv, Y4_cellline_summary.tsv
"""
import io, os, sys, glob, collections
import numpy as np
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of the 12 diagnostic variants)
from Bio import Align
# cell line units are also called with the 12 diagnostic variants. Align globally to CAN, call v4 tokens and count with TD.n_diag
#   (test Y4ALN_TEST_1031: 12 units with the context-1096 deletion = 12 units with alignment D1096:1, those 12 have 11-12 diagnostics, the other 103 units have 0-4)
_AL = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3, open_gap_score=-5, extend_gap_score=-1)
def aln_tokens(s):
    a = _AL.align(G.CAN, s)[0]
    t, q = a.aligned
    out = []
    for (t0, t1), (q0, q1) in zip(t, q):
        for i in range(t1 - t0):
            if G.CAN[t0 + i] != s[q0 + i] and s[q0 + i] in "ACGT": out.append("%d%s" % (t0 + i + 1, s[q0 + i]))
    for j in range(len(t) - 1):
        tg = t[j + 1][0] - t[j][1]; qg = q[j + 1][0] - q[j][1]
        if tg > 0: out.append("%dD%d" % (t[j][1] + 1, tg))
        if qg > 0: out.append("%dI%d" % (t[j][1] + 1, qg))
    return {y for y in (G.to_v4(x) for x in out) if y}
OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
SRC = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\P11_ENCODE_unit_placement\04_assembly_handover")
CAN = G.CAN; K = 31; W = 8
P1096, P635 = 1095, 634
ctx_G = CAN[P1096 - W:P1096] + "G" + CAN[P1096 + 1:P1096 + W + 1]
ctx_T = CAN[P1096 - W:P1096 + W + 1]
ctx_del = CAN[P1096 - W:P1096] + CAN[P1096 + 1:P1096 + W + 1]
ctx635_n = CAN[P635 - W:P635 + W + 1]
ctx635_T = CAN[P635 - W:P635] + "T" + CAN[P635 + 1:P635 + W + 1]
def read_fa(p):
    seqs = {}; name = None; buf = []
    for L in io.open(p, encoding="utf8"):
        if L.startswith(">"):
            if name: seqs[name] = "".join(buf).upper()
            name = L[1:].split()[0]; buf = []
        else: buf.append(L.strip())
    if name: seqs[name] = "".join(buf).upper()
    return seqs
def kset(s): return set(s[i:i + K] for i in range(len(s) - K + 1))
def jac(a, b): return len(a & b) / max(len(a | b), 1)
files = sorted(glob.glob(os.path.join(SRC, "*_distal20.fa")))
print("파일 %d : %s" % (len(files), ", ".join(os.path.basename(f) for f in files)))
g = io.open(os.path.join(OUT, "Y4_cellline_units.tsv"), "w", encoding="utf8", newline="\n")
g.write("assembly\tdepth\tunit_name\tlen\tcall_1096\tcall_635\ttype\tjac_to_prev\tn_diag12\n")
summ = []
for p in files:
    asm = os.path.basename(p).replace("_distal20.fa", "")
    seqs = read_fa(p)
    order = sorted(seqs, key=lambda n: int("".join(c for c in n.split("_pos")[-1][:3] if c.isdigit()) or 0))
    ks = {n: kset(seqs[n]) for n in order}
    types = []; prev = None
    for d, n in enumerate(order, 1):
        s = seqs[n]
        c1 = "del" if ctx_del in s else ("G" if ctx_G in s else ("T" if ctx_T in s else "?"))
        c2 = "T" if ctx635_T in s else ("norm" if ctx635_n in s else "?")
        kd = TD.n_diag(aln_tokens(s))
        if TD.MODE == "D1096": t = "말단형" if c1 == "del" else ("보통" if c1 in ("G", "T") else "판정불가")
        else: t = "말단형" if kd >= TD.TMIN else "보통"
        jp = jac(ks[n], ks[prev]) if prev else float("nan")
        g.write("%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%d\n" % (asm, d, n, len(s), c1, c2, t, ("%.4f" % jp) if prev else "—", kd))
        types.append(t); prev = n
    mk = [i + 1 for i, t in enumerate(types) if t == "말단형"]
    gaps = [mk[i + 1] - mk[i] for i in range(len(mk) - 1)]
    lag = {}
    for L in (1, 2, 3):
        v = [jac(ks[order[i]], ks[order[i + L]]) for i in range(len(order) - L)]
        lag[L] = float(np.mean(v)) if v else float("nan")
    # number of groups within a haplotype (Jaccard 0.85 single linkage, same threshold as for the human data)
    lab = list(range(len(order)))
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            if jac(ks[order[i]], ks[order[j]]) >= 0.85:
                a, b = lab[i], lab[j]
                if a != b: lab = [a if x == b else x for x in lab]
    ng = len(set(lab))
    summ.append((asm, len(order), ",".join(map(str, mk)) or "없음", ",".join(map(str, gaps)) or "—", lag[1], lag[2], lag[3], ng,
                 sum(1 for t in types if t == "판정불가")))
g.close()
g = io.open(os.path.join(OUT, "Y4_cellline_summary.tsv"), "w", encoding="utf8", newline="\n")
g.write("assembly\tn_units\tterminal_depths\tgaps\tjac_lag1\tjac_lag2\tjac_lag3\tn_groups\tn_uncalled\n")
print("\n=== 세포주 원위 20 유닛 ===")
print("  조립본                유닛  말단형 자리        간격      닮음 간격1  간격2  간격3   그룹수")
for r in summ:
    g.write("%s\t%d\t%s\t%s\t%.4f\t%.4f\t%.4f\t%d\t%d\n" % r)
    print("  %-20s %4d  %-16s %-8s %8.3f %6.3f %6.3f %7d" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]))
g.close()
allg = [x for r in summ for x in ([int(y) for y in r[3].split(",")] if r[3] != "—" else [])]
print("\n  모든 세포주 간격 모음 : " + (" · ".join("%d 칸 %d" % (k, v) for k, v in sorted(collections.Counter(allg).items())) if allg else "말단형 2 개 이상인 조립본 없음"))
d = [r[5] - r[4] for r in summ]
print("  간격2 − 간격1 (조립본별) : " + " · ".join("%+.3f" % x for x in d) + "  · 평균 %+.4f" % np.mean(d))
print("\n저장 :", OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('Y4_cellline_summary.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
