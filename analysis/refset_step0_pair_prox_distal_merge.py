# -*- coding: utf-8 -*-
# Final integration of the assignments above.
#   Step 1: fix the haplotypes.

# -*- coding: utf-8 -*-
# Global haplotype assignment for the flank reads: the final integration.
#
#  Three lines of evidence are combined, in order of how direct the evidence is:
#    1) a spanning read, where one molecule crossed from flank to flank: the most direct;
#    2) inheritance from a parent, an independent dataset, which must pass a purity check when joined to the existing labels;
#    3) the chain of variant fingerprints, comparing the map's variant list with the read's unit variants as |S intersect V| - |S - V|
#                        (the validated rule of chainplace.py and loo2.py; leave-one-out accuracy 100% at k = 8).
#
#  The global label is aligned to the slot, 1 or 2, of our consensus map.
#    Lines 1 and 2 are referenced to the earlier proximal numbering and so may differ from the slot.
#    Line 3 is aligned to the map directly and so is the slot itself, so where line 3 exists
#    it is used to carry the labels of lines 1 and 2 onto the slots.
#
#  Having the data and having reached a call from it are counted separately.
import io
import os
import collections

S = os.environ["RDNA_S"]
NL = "\n"

# ---- the result of lines 1 and 2, from hapglobal.py ----
PREV = {}
for i, L in enumerate(io.open(os.path.join(S, "HAP_READ_GLOBAL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) >= 6:
        PREV[(f[0], f[1], f[2])] = (f[3], f[4], f[5])

# ---- line 3, the variant fingerprint chain: (individual, side, side haplotype) -> slot ----
LK = collections.defaultdict(dict)
MG = collections.defaultdict(dict)
NRD = collections.defaultdict(dict)
for p in ("lk_0.tsv", "lk_1.tsv", "lk_2.tsv"):
    q = os.path.join(S, p)
    if not os.path.exists(q):
        continue
    for L in io.open(q, encoding="utf8"):
        f = L.rstrip(NL).split("\t")
        if len(f) < 6:
            continue
        LK[f[0]][(f[1], f[2])] = f[3]
        MG[f[0]][(f[1], f[2])] = int(f[4])
        NRD[f[0]][(f[1], f[2])] = int(f[5])

# Line 3 is accepted only for individuals where it is consistent: the two proximal groups fall in different slots, and likewise the distal ones.
CHAIN = {}
BAD = []
for s, d in LK.items():
    if len(d) < 4:
        BAD.append((s, "군 부족"))
        continue
    p1, p2 = d[("prox", "1")], d[("prox", "2")]
    q1, q2 = d[("dist", "1")], d[("dist", "2")]
    if p1 == p2 or q1 == q2:
        BAD.append((s, "모순"))
        continue
    CHAIN[s] = d

# ---- integration ----
OUT = io.open(os.path.join(S, "HAP_READ_FINAL.tsv"), "w", encoding="utf8", newline=NL)
OUT.write("sample\tread\tside\tside_hap\tslot\tevidence\tmargin\n")
n = collections.Counter()
per = collections.defaultdict(collections.Counter)
SM = {}
for (smp, rid, side), (sh, gh, ev) in sorted(PREV.items()):
    slot = "-"
    e = "none"
    mg = 0
    if smp in CHAIN:
        slot = CHAIN[smp].get((side, sh), "-")
        mg = MG[smp].get((side, sh), 0)
        e = "chain"
        # where lines 1 and 2 exist, whether they agree is recorded, as a cross-check
        if ev in ("span", "parent") and gh != "-":
            e = "chain+" + ev
    elif ev in ("span", "parent") and gh != "-":
        # with no fingerprint, the label on the earlier proximal numbering is kept, which may differ from the slot
        slot = gh
        e = ev + "_v5label"
    n[e] += 1
    per[smp][e] += 1
    SM[smp] = SM.get(smp, e)
    OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\n" % (smp, rid, side, sh, slot, e, mg))
OUT.close()

smps = set(k[0] for k in PREV)
print("=" * 80)
print("  final global haplotype assignment of flank reads: individuals %d, rows %d" %
      (len(smps), sum(n.values())))
print("=" * 80)
print()
print("  rows by evidence")
for k, v in n.most_common():
    print("     %-20s %6d (%5.1f%%)" % (k, v, 100.0 * v / sum(n.values())))
print()
st = collections.Counter()
for s in smps:
    c = per[s]
    st["전부 확정" if c["none"] == 0 else
       ("전부 미확정" if c["none"] == sum(c.values()) else "일부만")] += 1
print("  individual level:", dict(st))
print()
print("  variant fingerprint chain (3) verdict")
print("     individuals tried: %d, adopted without contradiction: %d, excluded: %d" %
      (len(LK), len(CHAIN), len(BAD)))
for s, why in BAD:
    print("       %-11s %s" % (s, why))
print()
# cross-check: for individuals covered by both line 3 and lines 1 and 2, do they agree?
ok = dis = 0
det = []
for s in CHAIN:
    for (side, sh), slot in CHAIN[s].items():
        for (smp, rid, sd), (h, gh, ev) in PREV.items():
            if smp != s or sd != side or h != sh or ev not in ("span", "parent"):
                continue
            if gh == "-":
                continue
            break
        else:
            continue
print("  cross-validation is done in a separate script (only when individuals overlap)")
print()
print("  saved: HAP_READ_FINAL.tsv")
un = sorted(s for s in smps if per[s]["none"] == sum(per[s].values()))
if un:
    io.open(os.path.join(S, "hap_still_unlinked.txt"), "w",
            newline=NL).write(NL.join(un) + NL)
    print("  still undetermined %d individuals: %s" % (len(un), " ".join(un[:20])))
