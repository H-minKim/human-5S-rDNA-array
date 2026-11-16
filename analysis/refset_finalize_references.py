# -*- coding: utf-8 -*-
# tier 1 first, then tier 2, with duplicates removed and at most 2 per individual
#   Step: fix the references.
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# Fix the reference per haplotype, which sets the ruler the methylation map is built on.
#
#  The tier rules:
#    tier 1 is the array interval of an array-spanning read, used wherever one exists;
#    tier 2 is a hifiasm reassembly of the array, used where there is no spanning read;
#    tier 3 is neither, and no reference can be built.
#
#  Choosing between several hifiasm contigs:
#    measured in HG002, hifiasm produced three, with hap1 at 68 units (correct), hap2 at 77 (spurious) and pri at 103 (correct).
#    Choosing by the hap1 / hap2 label is therefore wrong.
#    The spurious 77 was a partial duplicate overlapping the real 103 at 99.84%.
#    The rule: take them in decreasing length, discarding any that overlaps an already chosen one substantially, at 90% coverage or more.
#    At most 2 are kept per individual, one per haplotype.
import io
import os
import sys
import subprocess
import collections

W = os.path.expanduser("~/refbuild")
os.makedirs(W, exist_ok=True)
T1 = os.path.expanduser("~/tier1ref")
HA = os.path.expanduser("~/hiasm/arr")
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
NL = "\n"


def read_fa1(p):
    L = io.open(p, encoding="utf8").read().split(NL)
    return L[0][1:], "".join(L[1:])


def dedup(cands):
    """in decreasing order of length; one that overlaps an already chosen one substantially is dropped"""
    keep = []
    for nm, sq, src in sorted(cands, key=lambda z: -len(z[1])):
        dup = False
        for knm, ksq, ksrc in keep:
            a = os.path.join(W, "a.fa")
            b = os.path.join(W, "b.fa")
            io.open(a, "w", newline=NL).write(">a" + NL + sq + NL)
            io.open(b, "w", newline=NL).write(">b" + NL + ksq + NL)
            paf = subprocess.run(["minimap2", "-x", "asm10", "-t", TH,
                                  "--secondary=no", b, a],
                                 capture_output=True, text=True)
            cov = 0
            iv = []
            for L in paf.stdout.splitlines():
                f = L.split("\t")
                if len(f) >= 11:
                    iv.append((int(f[2]), int(f[3])))
            iv.sort()
            mg = []
            for s0, e0 in iv:
                if mg and s0 <= mg[-1][1]:
                    mg[-1][1] = max(mg[-1][1], e0)
                else:
                    mg.append([s0, e0])
            cov = sum(e - s for s, e in mg)
            if cov >= 0.90 * len(sq):
                dup = True
                break
        if not dup:
            keep.append((nm, sq, src))
        if len(keep) >= 2:
            break
    return keep


OUT = io.open(os.path.join(W, "refs_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
for smp in sys.argv[1:]:
    # tier 1 candidates
    t1 = []
    for f in sorted(os.listdir(T1)) if os.path.isdir(T1) else []:
        if f.startswith("t1_%s_" % smp) and f.endswith(".fa"):
            nm, sq = read_fa1(os.path.join(T1, f))
            t1.append((nm, sq, "tier1"))
    # tier 2 candidates
    t2 = []
    if os.path.isdir(HA):
        for f in sorted(os.listdir(HA)):
            if f.startswith("hiarr_%s_" % smp) and f.endswith(".fa"):
                nm, sq = read_fa1(os.path.join(HA, f))
                t2.append((nm, sq, "tier2"))
    if not t1 and not t2:
        OUT.write("%s\t-\tNONE\t-\t-%s" % (smp, NL))
        OUT.flush()
        continue
    # Tier 1 first; where there is only one tier 1 candidate, the remaining slot is filled from tier 2.
    picked = dedup(t1) if t1 else []
    if len(picked) < 2 and t2:
        more = dedup(picked + t2)
        for x in more:
            if len(picked) >= 2:
                break
            if not any(x[0] == y[0] for y in picked):
                picked.append(x)
    od = os.path.join(W, "ref")
    os.makedirs(od, exist_ok=True)
    for k, (nm, sq, src) in enumerate(picked):
        p = os.path.join(od, "%s_h%d.fa" % (smp, k + 1))
        io.open(p, "w", newline=NL).write(
            ">%s_h%d|%s|%s|%d%s%s%s" % (smp, k + 1, src, nm, len(sq), NL, sq, NL))
        OUT.write("%s\th%d\t%s\t%s\t%d%s" % (smp, k + 1, src, nm, len(sq), NL))
    OUT.flush()
    sys.stderr.write("  %s : tier1후보 %d · tier2후보 %d -> 채택 %s\n"
                     % (smp, len(t1), len(t2),
                        " ".join("%s(%d bp)" % (x[2], len(x[1])) for x in picked)))
    sys.stderr.flush()
OUT.close()
print("done")
