# -*- coding: utf-8 -*-
"""
EXPORT_origin_call - calls each variant **inherited / recurrent** and attaches the call to the group table.

Basis of the call (a direct extension of the lineage test of S3_zone_specific_lineage.py)
  Background: haplotypes with 15 or more MID units (more than 30 units from both junctions).
         Background markers = variants whose MID carrier fraction is 20-80% of haplotypes.
         Haplotype x background-marker matrix -> 8 lineages with KMeans(8, n_init=20, random_state=0).
  Lineage concentration: the fraction of carrier haplotypes that fall in one lineage (top_share)
    is compared with a null of the same number of randomly drawn haplotypes (500 permutations) -> p_lineage
    If the variant itself is a background marker, it is removed from the background and the lineages are rebuilt (no self-reference).

Called on three axes
  (1) Lineage concentration: p_lineage < 0.05 means the carriers are confined to particular lineages = descended from one ancestor
  (2) Copies per array: number of units per carrier array. Large means it arose once and its copies expanded
  (3) Recurrence-prone context: C>T/G>A at CpG (deamination) or an indel inside a repeat block (slippage)

Call rules (conservative, 'undetermined' when ambiguous)
  ancestral            : carrier haplotype fraction >= 0.90        (present in almost every lineage, so the lineage test is meaningless)
  inherited            : p_lineage < 0.05                     (confined to lineages)
  recurrence-suggestive: p_lineage >= 0.20 and copies per array <= 1.5 and recurrence-prone context
  dispersed            : p_lineage >= 0.20 but the above conditions not met (not confined to lineages, but weak evidence of recurrence)
  undetermined         : fewer than 8 carrier haplotypes, or 0.05 <= p_lineage < 0.20

Limit: this is a call on **traces of transmission**, not a count of mutation events.
  Counting de novo events directly needs trio (parent-child) data, which this edition does not have.
Output: columns added to export\V11_VARIANT_GROUPS.tsv (original kept as .noorigin), and export\V11_ORIGIN_CALL.tsv
"""
import sys, io, os, collections, shutil
import numpy as np
from sklearn.cluster import KMeans
import unitdef_v11 as FC
import v11_io as G

PRESET = sys.argv[1] if len(sys.argv) > 1 else "primary"
OUTD = os.path.join(G.OUT, "export")
U, A, M = G.load(PRESET)
rng = np.random.default_rng(13)          # same seed 13 as S3

# ---------- background: haplotypes with 15 or more MID units ----------
HAPS = []
MIDSET = {}
for k, rs in U.items():
    mid = [r for r in rs if r[0] > 30 and r[1] > 30]
    if len(mid) < 15: continue
    HAPS.append(k)
    c = collections.Counter()
    for r in mid: c.update(r[4])
    MIDSET[k] = (c, len(mid))
NH = len(HAPS)
print("배경 하플로타입 %d (MID >= 15 유닛)" % NH)

# background markers: MID carrier fraction is 20-80% of haplotypes
have = collections.Counter()
for k in HAPS:
    c, n = MIDSET[k]
    for x in c: have[x] += 1
MK = sorted([x for x in have if 0.20 * NH <= have[x] <= 0.80 * NH])
print("배경 마커 %d" % len(MK))
MKI = {x: i for i, x in enumerate(MK)}
X = np.zeros((NH, len(MK)), dtype=np.float32)
for i, k in enumerate(HAPS):
    c, n = MIDSET[k]
    for x in c:
        j = MKI.get(x)
        if j is not None: X[i, j] = 1.0

LAB_FULL = KMeans(n_clusters=8, n_init=20, random_state=0).fit(X).labels_
def lineages_excluding(x):
    """If the variant is a background marker, drop it and rebuild the lineages (no self-reference)"""
    if x not in MKI: return LAB_FULL
    keep = [j for j in range(len(MK)) if j != MKI[x]]
    return KMeans(n_clusters=8, n_init=20, random_state=0).fit(X[:, keep]).labels_

# ---------- carrier haplotypes per variant ----------
carry = collections.defaultdict(set)
nu = collections.Counter(); nh_all = collections.defaultdict(set)
for k, rs in U.items():
    for r in rs:
        for x in r[4]:
            nu[x] += 1; nh_all[x].add(k)
for k in HAPS:
    c, n = MIDSET[k]
    for x in c: carry[x].add(k)
# variants seen only outside MID are also included (carrier status is judged on the whole array)
carry_all = collections.defaultdict(set)
for k in HAPS:
    s = set()
    for r in U[k]: s |= r[4]
    for x in s: carry_all[x].add(k)

NHAPTOT = len(U)
def conc(carmask, B=500, lab=None):
    lab = LAB_FULL if lab is None else lab
    n = int(carmask.sum())
    c = np.array([carmask[lab == q].sum() for q in range(8)], float)
    top = c.max() / n
    nul = np.array([max(np.sum(lab[rng.permutation(NH)[:n]] == q) for q in range(8)) / float(n) for _ in range(B)])
    return float(top), float(np.median(nul)), float((np.sum(nul >= top) + 1) / (B + 1.0)), int(np.argmax(c))

MINHAP = 8
print("판정 대상 : 배경 하플로타입 중 보유 %d 개 이상" % MINHAP)
rows = []
cand = [x for x in carry_all if len(carry_all[x]) >= MINHAP]
print("대상 변이 %d / 전체 %d" % (len(cand), len(nu)))
for i, x in enumerate(sorted(cand, key=lambda z: -nu[z])):
    cm = np.array([1 if k in carry_all[x] else 0 for k in HAPS])
    nc = int(cm.sum())
    hapfrac = len(nh_all[x]) / float(NHAPTOT)
    per = nu[x] / float(len(nh_all[x]))
    cl, pos, ref, alt, ln, blk, rep = FC.parse(x)
    inrep = bool(rep or blk)
    recur_ctx = ((cl == "snv" and (ref, alt) in (("C", "T"), ("G", "A")) and FC.iscpg(pos))
                 or (cl in ("del", "ins") and inrep))
    if hapfrac >= 0.90:
        call = "조상형"; top = nl = p = float("nan"); which = -1
    else:
        lab = lineages_excluding(x) if nc >= MINHAP else LAB_FULL
        top, nl, p, which = conc(cm, 500, lab)
        if p < 0.05: call = "계승형"
        elif p >= 0.20 and per <= 1.5 and recur_ctx: call = "재발시사형"
        elif p >= 0.20: call = "분산형"
        else: call = "판정불가"
    rows.append((x, nc, len(nh_all[x]), nu[x], hapfrac, per, top, nl, p, which,
                 "Y" if recur_ctx else "N", call))
    if (i + 1) % 200 == 0: print("   ... %d / %d" % (i + 1, len(cand)))

p1 = os.path.join(OUTD, "V11_ORIGIN_CALL.tsv")
g = io.open(p1, "w", encoding="utf8", newline="\n")
g.write("variant\tn_carrier_bg_hap\tn_haplotypes\tn_units\thap_freq\tunits_per_carrier_array\t"
        "top_lineage_share\ttop_lineage_null\tp_lineage\ttop_lineage\trecurrence_context\torigin_call\n")
for r in rows:
    g.write("%s\t%d\t%d\t%d\t%.4f\t%.2f\t%s\t%s\t%s\t%d\t%s\t%s\n" % (
        r[0], r[1], r[2], r[3], r[4], r[5],
        ("%.4f" % r[6]) if r[6] == r[6] else "", ("%.4f" % r[7]) if r[7] == r[7] else "",
        ("%.4g" % r[8]) if r[8] == r[8] else "", r[9], r[10], r[11]))
g.close()
print("저장:", p1)

cnt = collections.Counter(r[11] for r in rows)
print("")
print("=== 판정 분포 (대상 %d) ===" % len(rows))
for c, n in cnt.most_common():
    sub = [r for r in rows if r[11] == c]
    print("  %-10s %5d  · 배열당 사본 중앙 %5.2f · 보유 하플로 중앙 %4.0f · 재발문맥 %4.1f%%" % (
        c, n, np.median([r[5] for r in sub]), np.median([r[2] for r in sub]),
        100 * np.mean([r[10] == "Y" for r in sub])))
print("")
print("=== 판정별 대표 변이 ===")
for c in ("계승형", "재발시사형", "분산형", "조상형"):
    sub = sorted([r for r in rows if r[11] == c], key=lambda r: -r[3])[:4]
    for r in sub:
        print("  %-10s %-12s 유닛 %6d · 배열 %3d · 배열당 %5.1f · 최다계통 %s (귀무 %s) p=%s" % (
            c, r[0], r[3], r[2], r[5],
            ("%.3f" % r[6]) if r[6] == r[6] else "-", ("%.3f" % r[7]) if r[7] == r[7] else "-",
            ("%.3g" % r[8]) if r[8] == r[8] else "-"))

# ---------- attach columns to the group table ----------
src = os.path.join(OUTD, "V11_VARIANT_GROUPS.tsv")
bak = src + ".noorigin"
if not os.path.exists(bak): shutil.copy2(src, bak)
D = {r[0]: r for r in rows}
hd = None; out = []
for L in io.open(bak, encoding="utf8"):
    f = L.rstrip("\n").split("\t")
    if hd is None:
        hd = f
        out.append(f + ["origin_call", "p_lineage", "top_lineage_share", "top_lineage_null",
                        "units_per_carrier_array", "recurrence_context"])
        continue
    r = D.get(f[0])
    if r is None:
        out.append(f + ["판정불가(보유 하플로 부족)", "", "", "", "%.2f" % (float(f[hd.index("n_units")]) / max(float(f[hd.index("n_haplotypes")]), 1)), ""])
    else:
        out.append(f + [r[11], ("%.4g" % r[8]) if r[8] == r[8] else "",
                        ("%.4f" % r[6]) if r[6] == r[6] else "", ("%.4f" % r[7]) if r[7] == r[7] else "",
                        "%.2f" % r[5], r[10]])
g = io.open(src, "w", encoding="utf8", newline="\n")
for f in out: g.write("\t".join(f) + "\n")
g.close()
print("")
print("그룹표에 열 6 개 추가:", src, " (원본 보존:", os.path.basename(bak) + ")")
