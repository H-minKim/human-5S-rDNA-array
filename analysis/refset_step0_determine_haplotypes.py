# -*- coding: utf-8 -*-
# Fix exactly 2 per individual, writing hap_final2.tsv, which stage 3 reads.
#   Step 1: fix the haplotypes.

# -*- coding: utf-8 -*-
# Stage 1: fix the haplotypes, exactly 2 per individual.
#
#  Three things the inputs require:
#
#  1) The public-assembly table of 44 individuals has a different file name, asm44_all.tsv,
#     and is read as well: 270 + 44 = 314.
#     Note: these 44 have a public assembly only, with no hifiasm, spanning or recombined-spanning data.
#
#  2) There are two frames for the unit count.
#     The unit column of the source tables (pub_all, asm44_all, asm_units_all, hisum_v2, refseq_headers)
#     counts normal gaps of 2,100-2,400 bp, whereas cmp4 reports gaps + 1, the true number of units.
#     Measured: 595 of 610 public, 708 of 728 hifiasm and 178 of 210 spanning differ by exactly 1.
#     The source values are therefore all in the gap frame and are converted to the unit frame by adding 1.
#       The conversion happens in add_units() in this file alone.
#
#  3) Counting noise in the spanning reads must not be taken for a second haplotype.
#     The distribution of gaps between the strong spanning reads within an individual is clearly bimodal:
#       gap 0: 194, 1: 21, 2: 8, then 3: 3, 4: 2, 5: 1, 6: 2 and a flat tail.
#     A gap of 2 or less is counting noise on the same molecule; 3 or more is a genuinely different haplotype.
#     So the strong spanning reads are first clustered by single linkage at a gap of 2 or less and reduced to one median.
#
#  The priority: where an array-spanning read exists, it is the answer.
#    An array-spanning read is direct evidence that one molecule physically crossed from flank to flank.
#    A spanning read alone can still be wrong, so it ranks first only when another method supports it.
#      1) a strong spanning read, with both flanks confirmed, supported by another method: first
#      2) the public assembly
#      3) a strong spanning read alone
#      4) otherwise, by the most cross-support, then hifiasm, then spanning, then recombined-spanning
#
#    Why the public assembly is not given absolute priority: of the 12 cases where it lost,
#    the spanning reads supported the other answer in every one, and
#    HG00128's 81-unit public assembly array carries a 44,784 bp gap.
#
#  A strong spanning read is defined as:
#    kind == "SPAN", with hits in both the proximal and distal 50 kb flanks, chain continuity of at least 95%, and at least 19 normal gaps.
import io
import os
import re
import glob
import collections

S = os.environ["RDNA_S"]
D = os.environ["RDNA_ASM_COMPARE"]
U = os.environ["RDNA_UNIFIED_REFMAP"]
NL = "\n"
TOL = 1          # the difference in unit count within which two values count as the same haplotype
SNOISE = 2       # the width of the counting noise in spanning reads, the boundary measured in the bimodal distribution
PRI = {"공개조립본": 0, "우리hifiasm": 1, "관통": 2, "재조합관통": 3}

UNITS = collections.defaultdict(lambda: collections.defaultdict(set))
PUBCTG = collections.defaultdict(dict)      # sample -> contig -> unit count, per assembly contig


def add_units(smp, method, value, frame):
    """Conversion of the unit-count frame, done only in this function.
       frame='GAPS'  : number of regular intervals -> units = intervals + 1
       frame='UNITS' : already a unit count        -> as is"""
    UNITS[smp][method].add(value + 1 if frame == "GAPS" else value)


# a) the public assemblies, including asm44_all.tsv, in the gap frame
for p in (os.path.join(D, "asm_units_all.tsv"),
          os.path.join(S, "pub_all.tsv"),
          os.path.join(S, "asm44_all.tsv")):
    if not os.path.exists(p):
        continue
    for L in io.open(p, encoding="utf8"):
        f = L.rstrip(NL).split("\t")
        if len(f) >= 9 and f[8] == "OK":
            try:
                add_units(f[0], "공개조립본", int(f[7]), "GAPS")
                PUBCTG[f[0]][f[3]] = int(f[7]) + 1      # also kept per contig
            except ValueError:
                pass

# b) our hifiasm assemblies, in the gap frame
for L in io.open(os.path.join(S, "hisum_v2.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) >= 6 and f[1] == "OK":
        s = re.sub(r"_\d+$", "", f[0])                    # strip the trailing _N first
        for pat in (r"_(sw|f)[2-6]$", r"_l0$", r"_sub250$", r"_noSPAN$"):
            s = re.sub(pat, "", s)
        try:
            add_units(s, "우리hifiasm", int(f[5]), "GAPS")
        except ValueError:
            pass

# c) spanning and d) recombined-spanning reads, from the refseq headers, in the gap frame
for L in io.open(os.path.join(S, "refseq_headers.txt"), encoding="utf8"):
    f = L.strip().split("|")
    if len(f) < 4:
        continue
    m = re.match(r"^N(\d+)$", f[3])
    if not m:
        continue
    k = "관통" if f[2] == "tier1" else ("재조합관통" if f[2] == "tier2" else None)
    if k:
        add_units(f[0], k, int(m.group(1)), "GAPS")

# the value cmp4 measured directly from the sequence, in the unit frame
KN = {"pub": "공개조립본", "hifi": "우리hifiasm", "span": "관통"}
for L in io.open(os.path.join(S, "cmp4_all.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) < 9:
        continue
    for col, kk in ((3, f[1]), (4, f[2])):
        if kk not in KN or f[col] in ("-", ""):
            continue
        try:
            add_units(f[0], KN[kk], int(f[col]), "UNITS")
        except ValueError:
            pass

# ---- strong spanning reads, both flanks confirmed, in the gap frame, with the noise merged ----
RAW = collections.defaultdict(list)
for L in io.open(os.path.join(S, "verify3_all.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) < 10 or f[1] != "SPAN":
        continue
    try:
        px, dx, frac, gaps = int(f[4]), int(f[5]), float(f[8]), int(f[9])
    except ValueError:
        continue
    if px > 0 and dx > 0 and frac >= 95.0 and gaps >= 19:
        RAW[f[0]].append(gaps + 1)

STRONG = {}
NREAD = {}
for s, vs in RAW.items():
    vs.sort()
    grp = [[vs[0]]]
    for v in vs[1:]:
        if v - grp[-1][-1] <= SNOISE:      # a gap of 2 or less is counting noise on the same molecule
            grp[-1].append(v)
        else:
            grp.append([v])
    reps = [g[len(g) // 2] for g in grp]
    STRONG[s] = set(reps)
    NREAD[s] = dict((g[len(g) // 2], len(g)) for g in grp)
    for r in reps:
        add_units(s, "관통", r, "UNITS")

# ---- cohort ----
hprc = set()
for p in glob.glob(os.path.join(U, "refmap_hap_ip-*.tsv")):
    for i, L in enumerate(io.open(p, encoding="utf8")):
        if i == 0:
            continue
        f = L.rstrip(NL).split("\t")
        if len(f) >= 6:
            hprc.add(f[0])
EX44 = set(x.strip() for x in io.open(os.path.join(S, "ex44.txt"), encoding="utf8") if x.strip())
kgp = set(x.strip() for x in io.open(os.path.join(S, "L_kgp.txt"), encoding="utf8") if x.strip())
kgp -= (hprc | EX44)          # an individual present in both cohorts is not excluded
TEST = ("PARAM_", "D5_", "D50_", "D200_")
smps = sorted(s for s in UNITS if s not in kgp and not s.startswith(TEST))

TN = {0: "강한관통+교차지지", 1: "공개조립본", 2: "강한관통 단독", 3: "그 외(교차지지順)"}


def tier_of(c):
    if c["strong"] and c["n"] >= 2:
        return 0
    if "공개조립본" in c["m"]:
        return 1
    if c["strong"]:
        return 2
    return 3


def cluster(smp):
    pts = [(v, k) for k, vs in UNITS[smp].items() for v in vs]
    pts.sort(key=lambda z: -z[0])
    cl = []
    for v, k in pts:
        for c in cl:
            if abs(v - c["u"]) <= TOL:
                c["m"].add(k)
                c["v"].append(v)
                c["u"] = sorted(c["v"])[len(c["v"]) // 2]
                break
        else:
            cl.append({"u": v, "v": [v], "m": set([k])})
    for c in cl:
        c["n"] = len(c["m"])
        c["pri"] = min(PRI[x] for x in c["m"])
        c["strong"] = False
    # One representative value from the strong spanning reads attaches to the single nearest cluster.
    #   Marking by the tolerance of 2 alone would make two clusters 2 units apart both count as strongly spanned,
    #   letting one molecule occupy both haplotypes (as in HG03193 at 69 and 67).
    for x in STRONG.get(smp, set()):
        near = min(cl, key=lambda c: abs(c["u"] - x))
        if abs(near["u"] - x) <= SNOISE:
            near["strong"] = True
    cl.sort(key=lambda c: (tier_of(c), -c["n"], c["pri"], -c["u"]))
    return cl


def pick(smp):
    """Choose the two slots.

    Twin exception: in the truth set 4.4% of true haplotype pairs are within 1 unit of each other.
    Clustering (TOL=1) merges such a pair into one cluster, and the second slot would be filled by an unrelated value
    (3 measured cases: HG02040 truth 70,71 -> 71,83 / NA20805 75,76 -> 45,76 /
     HG01258 74,75 -> 72,74).
    Two different contigs in the public assembly are direct evidence of two molecules,
    so when the two contigs are close enough to fall in one cluster, their two values fill the two slots."""
    cl = cluster(smp)
    v = sorted(PUBCTG.get(smp, {}).values())
    if len(v) == 2 and v[1] - v[0] <= TOL and v[1] != v[0]:
        src = dict((c["u"], c) for c in cl)
        base = min(cl, key=lambda c: abs(c["u"] - v[1])) if cl else None
        if base is not None and abs(base["u"] - v[1]) <= TOL:
            out = []
            for x in (v[1], v[0]):
                d = dict(base)
                d["u"] = x
                d["twin"] = True
                out.append(d)
            return out, cl
    return cl[:2], cl


print("=" * 86)
print("  Stage 1 v2. Haplotype determination: individuals %d  (v1: 270)" % len(smps))
print("=" * 86)
print("  cohort  HPRC ONT %d, newly-extracted44 %d, other %d" %
      (len(set(smps) & hprc), len(set(smps) & EX44), len(set(smps) - hprc - EX44)))
print("  individuals with strong spanning (both flanks) %d, of which with 2 or more distinct spanning groups %d" %
      (len([s for s in smps if STRONG.get(s)]),
       len([s for s in smps if len(STRONG.get(s, ())) >= 2])))
print()

nc = collections.Counter()
sup = collections.Counter()
tiers = collections.Counter()
final = {}
ntwin = 0
for s in smps:
    got, cl = pick(s)
    nc[min(len(cl), 5)] += 1
    final[s] = got
    if any(c.get("twin") for c in got):
        ntwin += 1
    for c in got:
        sup[c["n"]] += 1
        tiers[tier_of(c)] += 1

print("  distribution of candidate cluster counts (before selection)")
for k in sorted(nc):
    print("     %-3s : %3d individuals (%5.1f%%)" %
          ("5+" if k == 5 else k, nc[k], 100.0 * nc[k] / len(smps)))
print()
n2 = sum(1 for s in smps if len(final[s]) >= 2)
tot = sum(tiers.values())
print("  individuals with 2 slots filled: %d / %d (%.1f%%), 1 slot %d, 0 slots %d" %
      (n2, len(smps), 100.0 * n2 / len(smps),
       sum(1 for s in smps if len(final[s]) == 1),
       sum(1 for s in smps if not final[s])))
print("  individuals rescued by the twin exception: %d" % ntwin)
print("  confirmed haplotypes: %d / %d slots = %.1f%%" %
      (tot, 2 * len(smps), 100.0 * tot / (2 * len(smps))))
print()
print("  basis of adoption")
for k in sorted(tiers):
    print("     %-20s %4d (%5.1f%%)" % (TN[k], tiers[k], 100.0 * tiers[k] / tot))
print()
print("  cross-support (how many methods pointed to the same unit count)")
for k in sorted(sup, reverse=True):
    print("     %d methods: %4d (%5.1f%%)" % (k, sup[k], 100.0 * sup[k] / tot))
print()
print("  source of haplotypes supported by a single method")
one = collections.Counter()
for s in smps:
    for c in final[s]:
        if c["n"] == 1:
            one[list(c["m"])[0]] += 1
for k, v in one.most_common():
    print("     %-12s %d" % (k, v))
print()

# ---- cases where the public assembly lost ----
print("  === cases where the public assembly was overruled (spanning won) ===")
print("  %-11s %-14s %-14s %s" % ("개체", "채택", "버린 공개조립본", "strong spanning (read count)"))
print("  " + "-" * 76)
lost = 0
for s in smps:
    got, cl = pick(s)
    if any(c.get("twin") for c in got):
        continue
    if len(cl) <= 2:
        continue
    drop = [c for c in cl[2:] if "공개조립본" in c["m"]]
    if not drop:
        continue
    lost += 1
    print("  %-11s %-14s %-14s %s" %
          (s, ",".join("%d" % c["u"] for c in got),
           ",".join("%d" % c["u"] for c in drop),
           ",".join("%d(%d)" % (x, NREAD[s][x]) for x in sorted(STRONG.get(s, ())))
           or "없음"))
print("  individuals where the public assembly was overruled: %d" % lost)
print()

# ---- cases where the spanning read and hifiasm disagree outright, for inspection by hand ----
print("  === cases where spanning alone beat hifiasm+recombined spanning (review needed) ===")
flag = 0
for s in smps:
    got, cl = pick(s)
    if not any(tier_of(c) == 2 for c in got):
        continue
    drop = [c for c in cl[2:] if c["n"] >= 2]
    if not drop:
        continue
    flag += 1
    print("  %-11s adopted %-12s discarded %s" %
          (s, ",".join("%d(%s)" % (c["u"], TN[tier_of(c)][:6]) for c in got),
           " · ".join("%d(%s)" % (c["u"], ",".join(sorted(c["m"]))) for c in drop)))
print("  to review: %d" % flag)

OUT = io.open(os.path.join(S, "hap_final2.tsv"), "w", encoding="utf8", newline=NL)
OUT.write("sample\tslot\tunits\tn_method\tstrong_span\tn_span_read\ttier\tmethods\n")
for s in smps:
    for i, c in enumerate(final[s]):
        nr = 0
        for x, k in NREAD.get(s, {}).items():
            if abs(c["u"] - x) <= SNOISE:
                nr = max(nr, k)
        OUT.write("%s\t%d\t%d\t%d\t%s\t%d\t%s\t%s\n" %
                  (s, i + 1, c["u"], c["n"], "Y" if c["strong"] else "N", nr,
                   TN[tier_of(c)], ",".join(sorted(c["m"]))))
OUT.close()
print()
print("  saved: hap_final2.tsv")
