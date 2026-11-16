# -*- coding: utf-8 -*-
# Global assignment of the proximal to distal pairing.
#   Step 1: fix the haplotypes.

# -*- coding: utf-8 -*-
# Global haplotype assignment for the flank reads.
#
#  The rule follows the canonical one in consensus_v2.py:
#    where an array-spanning read exists, that molecule fixes the proximal-to-distal pairing, and the labels must never be assumed to correspond.
#
#  The proximal and distal haplotype labels are independent of each other.
#    Measured: of 575 reads carrying both, the labels agree in 45.4%, which is chance.
#    Joining on the label would be wrong half the time; the link must come from a molecule, that is a spanning read, or from a parent.
#
#  Priority:
#    1) inheritance from a parent, which is globally consistent (measured: reads on both sides agree in 13 of 13);
#    2) the pairing fixed by a spanning read,
#       with elimination: there are only two haplotypes, so knowing one link determines the other
#         (without this step 139 individuals were left only partly resolved);
#    3) with neither, it is left unresolved and linked by the haplotype sequence at the next stage.
#
#  The global label is taken from the proximal numbering: proximal haplotypes 1 and 2 become global 1 and 2.
import io
import os
import collections

S = os.environ["RDNA_S"]
NL = "\n"


def load(p):
    R = collections.defaultdict(dict)
    for i, L in enumerate(io.open(os.path.join(S, p), encoding="utf8")):
        if i == 0:
            continue
        f = L.rstrip(NL).split("\t")
        if len(f) < 4:
            continue
        q = f[2].split("|")
        rid = q[1] if len(q) > 1 else f[2]
        R[(f[0], rid)][f[1]] = f[3]
    return R


V = load("hap_reads_v5.tsv")
P = load("parent_hap_full.tsv")

SPAN = set()
for i, L in enumerate(io.open(os.path.join(S, "units_all.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) >= 3 and f[2] == "SPAN":
        SPAN.add((f[0], f[1]))

# ---- the proximal-to-distal link, per individual ----
LINK = {}          # smp -> {prox_hap: dist_hap}
CONF = set()
raw = collections.defaultdict(collections.Counter)
for (s, r), v in V.items():
    if (s, r) not in SPAN or len(v) < 2:
        continue
    raw[s][(v["prox"], v["dist"])] += 1
for s, c in raw.items():
    d = {}
    bad = False
    for (ph, dh), n in c.most_common():
        if ph in d and d[ph] != dh:
            bad = True
        elif dh in d.values() and d.get(ph) != dh:
            bad = True
        else:
            d[ph] = dh
    if bad:
        CONF.add(s)
        continue
    # elimination: there are only two haplotypes
    if len(d) == 1:
        ph, dh = list(d.items())[0]
        op = "2" if ph == "1" else "1"
        od = "2" if dh == "1" else "1"
        d[op] = od
    LINK[s] = d

# ---- bring the parental labels onto the proximal numbering ----
# The parental labels are global but their numbering may differ, so the correspondence is established per individual and side by a majority vote.
PMAP = {}
for side in ("prox", "dist"):
    m = collections.defaultdict(collections.Counter)
    for (s, r), v in P.items():
        w = V.get((s, r))
        if not w or side not in v or side not in w:
            continue
        m[s][(w[side], v[side])] += 1
    for s, c in m.items():
        d = collections.defaultdict(collections.Counter)
        for (a, b), n in c.items():
            d[a][b] += n
        mp = {}
        pur = []
        for a in d:
            b, n = d[a].most_common(1)[0]
            mp[a] = b
            pur.append(n / float(sum(d[a].values())))
        PMAP[(s, side)] = (mp, min(pur) if pur else 0.0, sum(c.values()))

# linking through the parent: proximal label -> parental label -> distal label
PLINK = {}
for s in set(k[0] for k in PMAP):
    a = PMAP.get((s, "prox"))
    b = PMAP.get((s, "dist"))
    if not a or not b or a[1] < 0.8 or b[1] < 0.8 or a[2] < 4 or b[2] < 4:
        continue
    inv = {}
    ok = True
    for dl, pl in b[0].items():
        if pl in inv:
            ok = False
        inv[pl] = dl
    if not ok:
        continue
    d = {}
    for pl, par in a[0].items():
        if par in inv:
            d[pl] = inv[par]
    if len(d) == 1:
        ph, dh = list(d.items())[0]
        d["2" if ph == "1" else "1"] = "2" if dh == "1" else "1"
    if d:
        PLINK[s] = d

OUT = io.open(os.path.join(S, "HAP_READ_GLOBAL.tsv"), "w", encoding="utf8", newline=NL)
OUT.write("sample\tread\tside\tside_hap\tglobal_hap\tevidence\n")
n = collections.Counter()
per = collections.defaultdict(collections.Counter)
for (s, r), v in sorted(V.items()):
    lk = LINK.get(s)
    pk = PLINK.get(s)
    for side, h in sorted(v.items()):
        g = "-"
        ev = "none"
        use = lk if lk else pk
        src = "span" if lk else ("parent" if pk else None)
        if use:
            if side == "prox":
                g, ev = h, src
            else:
                inv = {dv: pv for pv, dv in use.items()}
                if h in inv:
                    g, ev = inv[h], src
        n[ev] += 1
        per[s][ev] += 1
        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (s, r, side, h, g, ev))
OUT.close()

smps = set(k[0] for k in V)
print("=" * 78)
print("  flank read 전역 하플로타입 배정 — 개체 %d · read x side %d" %
      (len(smps), sum(n.values())))
print("=" * 78)
print()
print("  개체별 연결 근거")
print("     관통(SPAN)으로 연결 : %d  (모순 %d 제외)" % (len(LINK), len(CONF)))
print("     부모로 연결        : %d  (관통 없는 개체 중 %d)" %
      (len(PLINK), len(set(PLINK) - set(LINK))))
print("     둘 중 하나       : %d / %d" % (len(set(LINK) | set(PLINK)), len(smps)))
print("     미확정           : %d" % len(smps - set(LINK) - set(PLINK)))
print()
print("  read x side 배정")
for k, v in n.most_common():
    print("     %-8s %6d (%5.1f%%)" % (k, v, 100.0 * v / sum(n.values())))
print()
st = collections.Counter()
for s in smps:
    c = per[s]
    st["전부 확정" if c["none"] == 0 else
       ("전부 미확정" if c["none"] == sum(c.values()) else "일부만")] += 1
print("  개체 단위 :", dict(st))
un = sorted(s for s in smps if per[s]["none"] == sum(per[s].values()))
io.open(os.path.join(S, "hap_unlinked.txt"), "w", newline=NL).write(NL.join(un) + NL)
print("  저장 : HAP_READ_GLOBAL.tsv · hap_unlinked.txt (%d 개체)" % len(un))
