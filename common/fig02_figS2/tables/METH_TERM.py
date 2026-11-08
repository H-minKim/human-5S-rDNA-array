# -*- coding: utf-8 -*-
r"""Methylation of terminal-type vs non-terminal-type units (2026-10-13, user question)
   After normalising at the haplotype level (unit beta - median unit beta of that haplotype)
   (1) at the same distal depth, terminal-type units vs non-terminal-type units (between haplotypes)
   (2) within the same haplotype, a terminal-type unit vs its immediate neighbour (depth +-1) non-terminal-type units (pairs)
   Ruler: hard-call beta (beta_call, dedup=pool), unit reads >= 5, cohort 541 haplotypes of 299 individuals (V11_SUBSET, repro_tables SUBSET541),
        terminal type = at least 10 of the 12 diagnostic variants (termdef.is_term, 2026-10-31)
   2026-11-07: the '# ' header lines written into the tables below still say 'H1 557' and 'D1096:1' (old wording), kept to preserve the canonical table bytes
   Read only, outputs go to scratchpad"""
import os, sys, io, collections
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of the 12 diagnostic variants)
OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
MINR = 5
U, _, META = G.load("h1")
M = G.load_meth(min_reads=MINR)
tok = collections.Counter(t for v in U.values() for r in v for t in r[4] if "1096" in t)
print("1096 이 든 토큰 :", tok.most_common(5))
TERM = "D1096:1"   # unused (old definition). Terminal type is decided by TD.is_term below
rows = []   # (slot, depth, N, is_term, beta)
for s, recs in U.items():
    for gp, gd, N, _, vs, _, info in recs:
        m = M.get((s[0], s[1], info["array_pos"]))
        if m is None or np.isnan(m["beta_call"]): continue
        rows.append((s, gd, N, TD.is_term(vs), m["beta_call"]))
print("메틸 있는 유닛 %d · 하플로타입 %d · 말단형 유닛 %d" % (len(rows), len({r[0] for r in rows}), sum(r[3] for r in rows)))
by = collections.defaultdict(list)
for r in rows: by[r[0]].append(r)
med = {s: np.median([r[4] for r in v]) for s, v in by.items() if len(v) >= 20}
ok = [r for r in rows if r[0] in med]
print("유닛 20 개 이상 하플로타입 %d · 유닛 %d" % (len(med), len(ok)))
# unit table for the figure (2026-10-14), depth 1-15, including haplotype-normalised values
with io.open(os.path.join(OUT, "METH_TERM_units.tsv"), "w", encoding="utf8", newline="\n") as gu:
    gu.write("# 말단형 vs 비말단형 유닛 메틸 · H1 557 · 하드콜 β(beta_call · dedup pool) · read>=%d · "
             "유닛 20 개 이상 하플로타입 · beta_norm = β - 그 하플로타입 유닛 β 중앙값 · 말단형 = D1096:1 보유\n" % MINR)
    gu.write("sample\tslot\tdepth\tis_term\tbeta\tbeta_norm\n")
    for s, d, N, t, b in ok:
        if d <= 15:
            gu.write("%s\t%s\t%d\t%d\t%.5f\t%.5f\n" % (s[0], s[1], d, int(t), b, b - med[s]))
D = collections.defaultdict(lambda: ([], [], [], []))   # depth -> (term_raw, non_raw, term_norm, non_norm)
for s, d, N, t, b in ok:
    if d > 15: continue
    x = D[d]; i = 0 if t else 1
    x[i].append(b); x[i + 2].append(b - med[s])
g = io.open(os.path.join(OUT, "METH_TERM_by_depth.tsv"), "w", encoding="utf8", newline="\n")
g.write("# 말단형 vs 비말단형 · H1 557 · 하드콜 β · read>=%d · 정규화 = β - 하플로타입 중앙값 · %s\n" % (MINR, META["preset"]))
g.write("depth\tn_term\tn_non\tterm_raw_med\tnon_raw_med\tterm_norm_med\tnon_norm_med\tdiff_norm_med\tMWU_p\n")
print("\n① 같은 깊이 (하플로타입 사이)")
print("깊이  n말단  n비말단   원값중앙(말/비)   정규화중앙(말/비)   차   p")
for d in sorted(D):
    tr, nr, tn, nn = D[d]
    if len(tr) < 5: continue
    p = stats.mannwhitneyu(tn, nn).pvalue
    dm = np.median(tn) - np.median(nn)
    print("%3d %6d %8d    %.3f / %.3f      %+.3f / %+.3f   %+.3f  %.2g" % (d, len(tr), len(nr), np.median(tr), np.median(nr), np.median(tn), np.median(nn), dm, p))
    g.write("%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4g\n" % (d, len(tr), len(nr), np.median(tr), np.median(nr), np.median(tn), np.median(nn), dm, p))
g.close()
# (2) pairs within the same haplotype: terminal-type unit beta - mean beta of neighbouring (+-1) non-terminal-type units
pair = []
for s, v in by.items():
    if s not in med: continue
    dd = {r[1]: r for r in v}
    for d, r in dd.items():
        if not r[3] or d > 15: continue
        nb = [dd[e][4] for e in (d - 1, d + 1) if e in dd and not dd[e][3]]
        if nb: pair.append((d, r[4] - np.mean(nb)))
a = np.array([x[1] for x in pair])
w = stats.wilcoxon(a)
print("\n② 같은 하플로타입 안 : 말단형 - 이웃 비말단형  짝 %d · 중앙 %+.4f · 평균 %+.4f · 양수 %.3f · Wilcoxon p %.3g"
      % (len(a), np.median(a), a.mean(), (a > 0).mean(), w.pvalue))
for d in (2, 5, 7):
    b = np.array([x[1] for x in pair if x[0] == d])
    if len(b) >= 5: print("   깊이 %d : 짝 %d · 중앙 %+.4f · 양수 %.3f · p %.3g" % (d, len(b), np.median(b), (b > 0).mean(), stats.wilcoxon(b).pvalue))
io.open(os.path.join(OUT, "METH_TERM_pairs.tsv"), "w", encoding="utf8", newline="\n").write(
    "depth\tterm_minus_neighbour\n" + "\n".join("%d\t%.5f" % x for x in pair) + "\n")

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('METH_TERM_units.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
