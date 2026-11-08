# -*- coding: utf-8 -*-
r"""
T1 - parent-to-child transmission of unit variants, tested with trios.
     (data = v11-addition2 trio version, 2026-09-20 combined handover note section 3)

What is examined
  (A) Find the transmitted haplotype - match the **per-unit variant sets** of child haplotypes against parent haplotypes.
      Null = all unrelated (child, parent) pairs.
  (B) Do units change within one generation - in the matched units, count variants found only in the child or only in the parent and look at their **position (depth from either end)**.
      If the S27 claim that 'distal 3-22 has frequent unit replacement and homogenisation' is right, more mismatches should appear distally.
  (C) Check of the 'rare single copy (R1)' proxy - the proxy used here is "at most 3 cohort haplotypes, 1 unit in the array" as
      a proxy for a recent event. If a child R1 variant **is also present in the transmitting parent, it was inherited** and is not a new event.
      Only those absent from the transmitting parent are de novo candidates.
  (D) Are the unit count and inherited-type blocks transmitted unchanged.

Data and frame
  Map: `G8MAP_trio.tsv` in the map folder of the v11-addition2 trio version (v11 cohort + trio combined, UNIT rows, verify OK)
  Child values come from the v11 cohort rows, parent values from the trio rows, **the same G8 frame**, so they can be matched by array_pos.
  Matching: shift the units of two haplotypes by -5..+5 anchored at the distal end (array_pos 1) and at the proximal end,
        and choose the offset with the largest mean per-unit variant Jaccard (unit counts may differ).
Limits (handover note section 5-2)
  - hg003 (father) -> hg002 is **excluded** because its canonical version is inferred. Output only as a control.
  - The child uses the v11-revision2 placement and the parent the chain AO placement. Only **variants of the canonical assembly** are used here, so the placement method does not matter.
  - Assembly errors may be present on both sides (mismatch = upper bound on de novo events).
Output: T1_pair_summary.tsv, T1_unit_mismatch.tsv, T1_R1_check.tsv, log
"""
import io, os, sys, collections, math
import numpy as np
V11CODE = __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir))
sys.path.insert(0, V11CODE)
import v11_io as G
import unitdef_v11 as FC

OUT = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\5S_position_variant_work\2026-09-20_trio")
MAPF = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11-add2_trio\03_maps\G8MAP_trio.tsv")
print("지도 :", MAPF, "md5", G.md5(MAPF))

# ---------- read the map ----------
SLOT = collections.defaultdict(list)     # (sample,slot) -> [(array_pos, N, set(tokens))]
for f, H in G._rows(MAPF):
    if f[H["seg_type"]] != "UNIT" or f[H["verify"]] != "OK": continue
    raw = [] if f[H["varlist"]] == "-" else f[H["varlist"]].split(",")
    if any(t.endswith("N") for t in raw): continue
    vs = set(x for x in (G.to_v4(t) for t in raw) if x is not None)
    SLOT[(f[0], f[H["slot"]])].append((int(f[H["array_pos"]]), int(f[H["N"]]), vs))
for k in SLOT: SLOT[k].sort()
print("슬롯 %d · 유닛 %d" % (len(SLOT), sum(len(v) for v in SLOT.values())))

FAM = [("HG01252", "HG01250", "부모1"), ("HG01252", "HG01251", "부모2"),
       ("HG00558", "HG00556", "부모"), ("HG03816", "HG03815", "부모"),
       ("HG01192", "HG01190", "부모"), ("HG01975", "HG01973", "부모"),
       ("hg002", "hg004", "모"), ("hg002", "hg003", "부(추론정본·제외)")]
CHILDREN = sorted(set(c for c, p, r in FAM))
PARENTS = sorted(set(p for c, p, r in FAM))

def jac(a, b):
    if not a and not b: return 1.0
    u = len(a | b)
    return len(a & b) / u if u else 1.0
def align(cs, ps):
    """child unit list cs, parent ps -> (best_mean_jac, mode, shift, pairs[(ci,pi)])"""
    best = None
    for mode in ("dist", "prox"):
        nc, npar = len(cs), len(ps)
        for sh in range(-5, 6):
            pairs = []
            for i in range(nc):
                j = i + sh if mode == "dist" else (npar - nc) + i + sh
                if 0 <= j < npar: pairs.append((i, j))
            if len(pairs) < 10: continue
            m = float(np.mean([jac(cs[i][2], ps[j][2]) for i, j in pairs]))
            ex = float(np.mean([cs[i][2] == ps[j][2] for i, j in pairs]))
            if best is None or m > best[0]: best = (m, ex, mode, sh, pairs)
    return best

g1 = io.open(os.path.join(OUT, "T1_pair_summary.tsv"), "w", encoding="utf8", newline="\n")
g1.write("child\tchild_slot\tchild_N\tparent\tparent_slot\tparent_N\trelation\tmode\tshift\tn_pairs\tmean_jaccard\tfrac_exact\t"
         "n_child_only\tn_parent_only\tper_unit_child_only\tper_unit_parent_only\n")
BEST = {}
print("\n=== (A) 자식 하플로타입 × 부모 하플로타입 — 유닛 변이로 맞추기 ===")
print("  자식        부모        관계        맞춤  자카드 평균  완전일치 비율  자식만/부모만 변이")
for ch in CHILDREN:
    for cs_slot in ("1", "2"):
        cs = SLOT.get((ch, cs_slot))
        if not cs: continue
        for pa in PARENTS:
            rel = next((r for c, p, r in FAM if c == ch and p == pa), "무관")
            for ps_slot in ("1", "2"):
                ps = SLOT.get((pa, ps_slot))
                if not ps: continue
                a = align(cs, ps)
                if a is None: continue
                m, ex, mode, sh, pairs = a
                co = sum(len(cs[i][2] - ps[j][2]) for i, j in pairs)
                po = sum(len(ps[j][2] - cs[i][2]) for i, j in pairs)
                g1.write("%s\t%s\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%d\t%d\t%.4f\t%.4f\n" % (
                    ch, cs_slot, len(cs), pa, ps_slot, len(ps), rel, mode, sh, len(pairs), m, ex, co, po,
                    co / len(pairs), po / len(pairs)))
                key = (ch, cs_slot)
                if key not in BEST or m > BEST[key][0]: BEST[key] = (m, ex, pa, ps_slot, rel, mode, sh, pairs, cs, ps)
for (ch, cs_slot), (m, ex, pa, ps_slot, rel, mode, sh, pairs, cs, ps) in sorted(BEST.items()):
    co = sum(len(cs[i][2] - ps[j][2]) for i, j in pairs); po = sum(len(ps[j][2] - cs[i][2]) for i, j in pairs)
    print("  %-9s s%s  %-9s s%s  %-10s %s%+d  %.4f      %.4f        %d / %d" % (ch, cs_slot, pa, ps_slot, rel, mode, sh, m, ex, co, po))
g1.close()

# related vs unrelated control
rows = [l.rstrip("\n").split("\t") for l in io.open(os.path.join(OUT, "T1_pair_summary.tsv"), encoding="utf8")][1:]
kin = [float(r[10]) for r in rows if r[6] != "무관"]; unk = [float(r[10]) for r in rows if r[6] == "무관"]
print("\n  자카드 평균 : 친족 짝 중앙 %.4f (최소 %.4f) · 무관 짝 중앙 %.4f (최대 %.4f) · 짝 %d/%d" % (
    np.median(kin), min(kin), np.median(unk), max(unk), len(kin), len(unk)))

# ---------- (B)(C)(D) confirmed transmission pairs only ----------
# Transmission pairs are not written by hand but chosen from the (A) result: the best (child slot -> parent slot) pair that is related and has Jaccard >= 0.95
CONF = []; CONF_EXTRA = []
for (ch_, cs_), (m_, ex_, pa_, ps_, rel_, mode_, sh_, pairs_, cs__, ps__) in sorted(BEST.items()):
    if rel_ == "무관" or m_ < 0.95: continue
    (CONF_EXTRA if pa_ == "hg003" else CONF).append((ch_, cs_, pa_, ps_))   # hg003 (inferred canonical version) is to be excluded
print("\n  확정 전달 짝 (친족 · 자카드>=0.95) : " + " · ".join("%s s%s<-%s s%s" % t for t in CONF) +
      "   | 참고(제외) : " + " · ".join("%s s%s<-%s s%s" % t for t in CONF_EXTRA))
g2 = io.open(os.path.join(OUT, "T1_unit_mismatch.tsv"), "w", encoding="utf8", newline="\n")
g2.write("child\tchild_slot\tparent\tparent_slot\tarray_pos_child\tap_prox\tap_dist\tside_depth\tchild_only\tparent_only\n")
g3 = io.open(os.path.join(OUT, "T1_R1_check.tsv"), "w", encoding="utf8", newline="\n")
g3.write("child\tchild_slot\tparent\tvariant\tn_hap_cohort\tin_parent_same_unit\tin_parent_array\tap_prox\tap_dist\torigin_call\n")

# cohort frequency and origin call (v11-revision2 stable version), used for the child R1 definition
U11, _, _ = G.load("primary", verbose=False)
hapcnt = collections.Counter()
for k, rs in U11.items():
    s = set()
    for r in rs: s |= r[4]
    hapcnt.update(s)
ORIGIN = {}
for i, Ln in enumerate(io.open(os.path.join(G.OUT, "export", "V11_ORIGIN_CALL.tsv"), encoding="utf8")):
    f = Ln.rstrip("\n").split("\t")
    if i == 0: H = {c: j for j, c in enumerate(f)}; continue
    ORIGIN[f[0]] = f[H["origin_call"]]
INH = {"계승형", "조상형"}

print("\n=== (B) 전달 짝에서 한 세대 뒤 유닛 변이 불일치 · (C) R1 대리지표 검증 ===")
tot_units = 0; tot_co = 0; tot_po = 0; depth_co = collections.Counter(); depth_n = collections.Counter()
R1sum = collections.Counter()
for ch, cs_slot, pa, ps_slot in CONF + CONF_EXTRA:
    key = (ch, cs_slot)
    cs = SLOT[(ch, cs_slot)]; ps = SLOT[(pa, ps_slot)]
    a = align(cs, ps)
    m, ex, mode, sh, pairs = a
    N = len(cs)
    arr = collections.Counter(t for _, _, vs in cs for t in vs)
    pall = set().union(*[vs for _, _, vs in ps]) if ps else set()
    co_n = po_n = 0
    for i, j in pairs:
        apos = cs[i][0]; ad = apos; ap = N - apos + 1
        d = ("근%d" % ap) if ap <= ad else ("원%d" % ad)
        conly = cs[i][2] - ps[j][2]; ponly = ps[j][2] - cs[i][2]
        co_n += len(conly); po_n += len(ponly)
        if ch == "hg002" and pa == "hg003": pass
        else:
            k2 = min(ap, ad); depth_n[k2 if k2 <= 30 else 99] += 1
            depth_co[k2 if k2 <= 30 else 99] += len(conly)
        if conly or ponly:
            g2.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\n" % (ch, cs_slot, pa, ps_slot, apos, ap, ad, d,
                     ",".join(sorted(conly)) or "-", ",".join(sorted(ponly)) or "-"))
    # (C) R1 check
    n_r1 = 0; in_same = 0; in_arr = 0
    for i, j in pairs:
        apos = cs[i][0]; ad = apos; ap = N - apos + 1
        for t in cs[i][2]:
            if ORIGIN.get(t, "") in INH: continue
            if arr[t] != 1 or hapcnt.get(t, 0) > 3: continue
            n_r1 += 1
            s1 = t in ps[j][2]; s2 = t in pall
            in_same += s1; in_arr += s2
            g3.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n" % (ch, cs_slot, pa, t, hapcnt.get(t, 0), int(s1), int(s2), ap, ad, ORIGIN.get(t, "미검정")))
    tag = " (제외 대상)" if (ch, pa) == ("hg002", "hg003") else ""
    print("  %-8s s%s ← %-8s s%s : 유닛 %3d (맞춤 %s%+d · 자카드 %.4f) · 자식만 %3d · 부모만 %3d · R1 %2d 중 같은 유닛에 %2d · 배열 어딘가 %2d%s" % (
        ch, cs_slot, pa, ps_slot, len(pairs), mode, sh, m, co_n, po_n, n_r1, in_same, in_arr, tag))
    if not tag:
        tot_units += len(pairs); tot_co += co_n; tot_po += po_n
        R1sum["n"] += n_r1; R1sum["same"] += in_same; R1sum["arr"] += in_arr
g2.close(); g3.close()
print("\n  확정 전달 5 짝 합계 : 유닛 %d · 자식에만 있는 변이 %d (유닛당 %.4f) · 부모에만 %d (유닛당 %.4f)" % (
    tot_units, tot_co, tot_co / tot_units, tot_po, tot_po / tot_units))
print("  R1(희귀 단일 사본) 합계 %d 개 중 전달 부모의 같은 유닛에 %d (%.0f%%) · 부모 배열 어딘가에 %d (%.0f%%)" % (
    R1sum["n"], R1sum["same"], 100.0 * R1sum["same"] / max(R1sum["n"], 1), R1sum["arr"], 100.0 * R1sum["arr"] / max(R1sum["n"], 1)))
print("\n  불일치의 위치 (가까운 끝에서의 깊이 · 유닛 1000 개당 자식-only 변이)")
for d in list(range(1, 11)) + [99]:
    if depth_n[d] == 0: continue
    print("    %-6s 유닛 %4d · 자식만 %3d · 1000 유닛당 %.1f" % ("깊이%d" % d if d != 99 else "깊이>30", depth_n[d], depth_co[d], 1000.0 * depth_co[d] / depth_n[d]))
print("\n저장 : %s\\T1_*.tsv" % OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('T1_pair_summary.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
