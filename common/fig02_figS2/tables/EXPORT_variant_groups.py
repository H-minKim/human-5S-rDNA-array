# -*- coding: utf-8 -*-
"""
EXPORT: exports the v11 variant group tables **so that they can be read as they are**.

What is exported
  1) V11_VARIANT_GROUPS.tsv: one row per variant. Group membership (9 axes) and frequency and position statistics.
     The group definition uses `groups()` of `S4_zone_groups.py` **unchanged** (not rebuilt ad hoc).
     9 axes: class (_rep suffix), rep, elem, win (100bp), sub, sub_ctx (CpG/non), ti_tv, elem_sub, indel (4 size classes)
  2) V11_REPEAT_BLOCKS.tsv: coordinates of the 71 repeat blocks (on canonical unit 1-2240)
  3) V11_GROUP_MEMBERS.tsv: (axis, group) -> list of variants in that group, and site count
The set is primary (cancer, stage B/E, break slots and hg002_2 excluded, verify OK). For another set, pass a preset as argument.
"""
import sys, io, os, collections
import numpy as np
import unitdef_v11 as FC
import v11_io as G

PRESET = sys.argv[1] if len(sys.argv) > 1 else "primary"
OUTD = os.path.join(G.OUT, "export")
os.makedirs(OUTD, exist_ok=True)
U, A, M = G.load(PRESET)

def isti(r, a): return (r + a) in ("AG", "GA", "CT", "TC")

def groups(x):
    """Identical, character for character, to groups() in S4_zone_groups.py"""
    cl, pos, ref, alt, ln, blk, rep = FC.parse(x)
    inrep = bool(rep or blk)
    o = [("class", cl + ("_rep" if inrep else "")), ("rep", "in_repeat" if inrep else "unique")]
    if pos is not None:
        o.append(("elem", FC.element(pos)))
        o.append(("win", "w%02d" % ((pos - 1) // 100)))
    if cl == "snv" and ref in "ACGT" and alt in "ACGT":
        sub = "%s>%s" % (ref, alt)
        o += [("sub", sub), ("sub_ctx", "%s|%s" % (sub, "CpG" if FC.iscpg(pos) else "non")),
              ("ti_tv", "transition" if isti(ref, alt) else "transversion"),
              ("elem_sub", "%s|%s" % (FC.element(pos), sub))]
    if cl in ("del", "ins"):
        sz = ("1bp" if ln == 1 else ("2-3bp" if ln <= 3 else ("4-9bp" if ln <= 9 else "10bp+")))
        o.append(("indel", "%s_%s" % (cl, sz)))
    return o

# ---------- per-variant statistics ----------
nu = collections.Counter()          # number of carrier units
nh = collections.defaultdict(set)   # carrier haplotypes
apv = collections.defaultdict(list); adv = collections.defaultdict(list)
wfreq = collections.defaultdict(list)   # frequency within the array
raws = collections.defaultdict(collections.Counter)   # original G8 token forms
TOT = 0
for k, rs in U.items():
    n = len(rs)
    c = collections.Counter()
    for r in rs: c.update(r[4])
    for r in rs:
        TOT += 1
        for x in r[4]:
            nu[x] += 1; nh[x].add(k)
            apv[x].append(r[0]); adv[x].append(r[1]); wfreq[x].append(c[x] / float(n))
        for t in r[6]["raw"]:
            v = G.to_v4(t)
            if v is not None: raws[v][t] += 1
NHAP = len(U)
print("유닛 %d · 하플로타입 %d · 고유 변이 %d" % (TOT, NHAP, len(nu)))

p1 = os.path.join(OUTD, "V11_VARIANT_GROUPS.tsv")
g = io.open(p1, "w", encoding="utf8", newline="\n")
g.write("\t".join([
    "variant", "g8_token_example", "class", "pos", "ref", "alt", "len",
    "in_repeat", "block_id", "element", "window100", "sub", "sub_ctx", "ti_tv", "elem_sub",
    "indel_size", "is_cpg", "is_coding", "mature5S_pos",
    "n_units", "n_haplotypes", "unit_freq", "hap_freq",
    "mean_within_array_freq", "mean_ap", "mean_ad", "frac_prox8", "frac_dist8", "frac_mid"]) + "\n")
GM = collections.defaultdict(list)
for x in sorted(nu, key=lambda z: -nu[z]):
    cl, pos, ref, alt, ln, blk, rep = FC.parse(x)
    d = dict(groups(x))
    for ax, gr in groups(x): GM[(ax, gr)].append(x)
    a = np.array(apv[x]); b = np.array(adv[x]); w = np.array(wfreq[x])
    ex = raws[x].most_common(1)[0][0] if raws[x] else ""
    g.write("\t".join([
        x, ex, cl, ("" if pos is None else str(pos)), ref, alt, str(ln),
        ("Y" if bool(rep or blk) else "N"), ("" if blk is None else str(blk)),
        d.get("elem", ""), d.get("win", ""), d.get("sub", ""), d.get("sub_ctx", ""),
        d.get("ti_tv", ""), d.get("elem_sub", ""), d.get("indel", ""),
        ("Y" if (pos is not None and FC.iscpg(pos)) else "N"),
        ("Y" if (pos is not None and pos >= 2122) else "N"),
        ("" if (pos is None or pos < 2122) else str(pos - 2121)),
        str(nu[x]), str(len(nh[x])), "%.6f" % (nu[x] / float(TOT)), "%.4f" % (len(nh[x]) / float(NHAP)),
        "%.4f" % w.mean(), "%.2f" % a.mean(), "%.2f" % b.mean(),
        "%.4f" % (a <= 8).mean(), "%.4f" % (b <= 8).mean(),
        "%.4f" % ((a > 30) & (b > 30)).mean()]) + "\n")
g.close()
print("저장:", p1)

p2 = os.path.join(OUTD, "V11_REPEAT_BLOCKS.tsv")
g = io.open(p2, "w", encoding="utf8", newline="\n")
g.write("block_id\tstart\tend\tlength\tsequence\telement\n")
for i, (a, b) in enumerate(G.BLOCKS, 1):
    g.write("%d\t%d\t%d\t%d\t%s\t%s\n" % (i, a, b, b - a + 1, G.CAN[a - 1:b], FC.element(a)))
g.close()
print("저장: %s  (블록 %d 개 · %d bp)" % (p2, len(G.BLOCKS), sum(b - a + 1 for a, b in G.BLOCKS)))

p3 = os.path.join(OUTD, "V11_GROUP_MEMBERS.tsv")
g = io.open(p3, "w", encoding="utf8", newline="\n")
g.write("axis\tgroup\tn_variants\tn_variants_ge3units\tn_unit_occurrences\tmembers\n")
for (ax, gr) in sorted(GM):
    mem = GM[(ax, gr)]
    ge3 = [x for x in mem if nu[x] >= 3]
    g.write("%s\t%s\t%d\t%d\t%d\t%s\n" % (ax, gr, len(mem), len(ge3), sum(nu[x] for x in mem),
                                          ",".join(sorted(mem, key=lambda z: -nu[z])[:200])))
g.close()
print("저장: %s  (축 %d · 그룹 %d)" % (p3, len(set(a for a, _ in GM)), len(GM)))
print("")
print("축별 그룹 수 :", dict(collections.Counter(a for a, _ in GM)))
