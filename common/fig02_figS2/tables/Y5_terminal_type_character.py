# -*- coding: utf-8 -*-
r"""
Y5: what exactly is different in terminal-type units (v11-fix2, primary, distal depth <=15)

What is measured
  1. The list of variants that almost always come together in terminal-type units (carrying D1096:1), carrier rate >=0.8, the 'defining set'
  2. Other variants skewed to the terminal type (top of terminal-type carrier rate - ordinary carrier rate)
  3. Variant counts per element: AluY(1-318), NTS1(319-1178), GA-rich(1179-1551), NTS2(1552-2031), promoter(2032-2121), gene(2122-2240)
     -> are the gene and promoter intact
  4. Unit length (sum of indels), methylation beta (compared with ordinary units at the same depth)
Output: $OUTDIR\Y5_terminal_defining_variants.tsv, Y5_terminal_by_element.tsv
"""
import io, os, sys, collections
import numpy as np
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of the 12 diagnostic variants)
OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
MARK = "D1096:1"; DMAX = 15
U, A_, M_ = G.load("primary", verbose=False)
MT = G.load_meth(min_reads=1, verbose=False)
BLOCKS = G.BLOCKS
ELEM = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA-rich", 1179, 1551), ("NTS2", 1552, 2031),
        ("promoter", 2032, 2121), ("gene", 2122, 2240)]
def pos_of(t):
    if t.startswith("DB") or t.startswith("IB"):
        b = int(t[2:].split(":")[0]); return BLOCKS[b - 1][0] if 1 <= b <= len(BLOCKS) else None
    if t[0] == "D": return int(t[1:].split(":")[0])
    if t[0] == "I": return int(t[1:].split(":")[0])
    if t[0] == "R": return int(t[1:-2])
    return int(t[:-2])
def elem_of(p):
    for nm, a, b in ELEM:
        if a <= p <= b: return nm
    return "?"
def indel_bp(t):
    if t.startswith("DB") or t.startswith("IB"): return (-1 if t[0] == "D" else 1) * int(t.split(":")[1])
    if t[0] == "D": return -int(t.split(":")[1])
    if t[0] == "I":
        q = t[1:].split(":")
        return int(q[1]) if len(q) > 1 and q[1].isdigit() else 1
    return 0
term = []; ordu = []
for k in sorted(U):
    for r in U[k]:
        if r[1] > DMAX: continue
        m = MT.get((k[0], k[1], r[6]["array_pos"]))
        b = m["beta_call"] if m and m["n_reads"] >= 3 else None
        rec = (r[4], r[1], b)
        (term if TD.is_term(r[4]) else ordu).append(rec)
print("말단형 유닛 %d · 보통 유닛 %d (원위 깊이 ≤%d)" % (len(term), len(ordu), DMAX))
ct = collections.Counter(); co = collections.Counter()
for v, d, b in term: ct.update(v)
for v, d, b in ordu: co.update(v)
NT, NO = len(term), len(ordu)
g = io.open(os.path.join(OUT, "Y5_terminal_defining_variants.tsv"), "w", encoding="utf8", newline="\n")
g.write("variant\tpos\telement\tfrac_terminal\tfrac_ordinary\tdiff\n")
rows = []
for t, n in ct.items():
    ft, fo = n / NT, co.get(t, 0) / NO
    p = pos_of(t)
    rows.append((t, p, elem_of(p) if p else "?", ft, fo, ft - fo))
rows.sort(key=lambda r: -r[5])
for r in rows[:200]: g.write("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\n" % r)
g.close()
print("\n=== 1·2. 말단형을 정의하는 변이 (보유율 차 상위 20) ===")
print("  변이          자리   구획        말단형   보통    차이")
for r in rows[:20]:
    print("  %-12s %5s %-10s %6.3f %6.3f %+6.3f" % (r[0], r[1], r[2], r[3], r[4], r[5]))
core = [r for r in rows if r[3] >= 0.8]
print("\n  말단형의 ≥80%% 가 가진 변이 %d 개 : %s" % (len(core), ", ".join(r[0] for r in core[:25])))
print("  그중 보통 유닛에서도 흔한 것(≥0.5) : %s" % (", ".join(r[0] for r in core if r[4] >= 0.5) or "없음"))
# 3. per element
def elem_count(recs):
    c = collections.Counter(); n = len(recs)
    for v, d, b in recs:
        for t in v:
            p = pos_of(t)
            if p: c[elem_of(p)] += 1
    return {k2: c[k2] / n for k2 in [e[0] for e in ELEM]}
et, eo = elem_count(term), elem_count(ordu)
g = io.open(os.path.join(OUT, "Y5_terminal_by_element.tsv"), "w", encoding="utf8", newline="\n")
g.write("element\tbp\tterminal_per_unit\tordinary_per_unit\tratio\tterminal_per_kb\tordinary_per_kb\n")
print("\n=== 3. 구획별 유닛당 변이 수 ===")
print("  구획        bp     말단형   보통    비(말단/보통)   말단형 /kb   보통 /kb")
for nm, a, b in ELEM:
    bp = b - a + 1
    g.write("%s\t%d\t%.4f\t%.4f\t%.3f\t%.3f\t%.3f\n" % (nm, bp, et[nm], eo[nm], et[nm] / max(eo[nm], 1e-9), 1000 * et[nm] / bp, 1000 * eo[nm] / bp))
    print("  %-10s %5d %8.3f %7.3f %12.2f %11.2f %9.2f" % (nm, bp, et[nm], eo[nm], et[nm] / max(eo[nm], 1e-9), 1000 * et[nm] / bp, 1000 * eo[nm] / bp))
# 4. length and methylation
lt = [2240 + sum(indel_bp(t) for t in v) for v, d, b in term]
lo = [2240 + sum(indel_bp(t) for t in v) for v, d, b in ordu]
print("\n=== 4. 유닛 길이 · 메틸 ===")
print("  길이 중앙 : 말단형 %.0f bp · 보통 %.0f bp" % (np.median(lt), np.median(lo)))
bt = [b for v, d, b in term if b is not None]; bo = [b for v, d, b in ordu if b is not None]
print("  β 평균(깊이 보정 전) : 말단형 %.3f (n %d) · 보통 %.3f (n %d)" % (np.mean(bt), len(bt), np.mean(bo), len(bo)))
for dep in (2, 5, 7, 9, 12):
    a = [b for v, d, b in term if d == dep and b is not None]
    c = [b for v, d, b in ordu if d == dep and b is not None]
    if len(a) >= 5 and len(c) >= 5:
        print("   깊이 %2d : 말단형 %.3f (n %3d) · 보통 %.3f (n %3d)" % (dep, np.mean(a), len(a), np.mean(c), len(c)))
print("\n저장 :", OUT)
