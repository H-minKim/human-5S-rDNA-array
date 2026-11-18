# -*- coding: utf-8 -*-
# A read-level bootstrap interval plus the between-group difference, also bootstrapped, for the cluster of 6 chr4 arrays.
import os
import io
import random
import collections

random.seed(0)
NB = 5000
TAB = chr(9)
d = collections.defaultdict(list)
W = os.environ["NODE_WORK_DIR"]   # the compute node's working directory
for L in io.open(os.path.join(W, "GSC_PERREAD.tsv"), encoding="utf8"):
    f = L.rstrip("\n").split(TAB)
    if f[0] == "group":
        continue
    d[(f[0], f[1], f[2])].append((int(f[4]), int(f[5])))


def frac(v):
    c = sum(x[0] for x in v)
    return sum(x[1] for x in v) / c if c else float("nan")


def boot(v):
    out = []
    n = len(v)
    for _ in range(NB):
        s = [v[random.randrange(n)] for _ in range(n)]
        out.append(frac(s))
    return out          # No alignment is performed: a between-group difference must be taken between independent samples.


def ci(bs):
    t = sorted(bs)
    return t[int(0.025 * NB)], t[int(0.975 * NB)]


out = io.open(os.path.join(W, "GSC_STATS.tsv"), "w", encoding="utf8", newline="\n")
out.write(TAB.join(["group", "branch", "part", "n_read", "n_cpg", "frac", "lo95", "hi95", "minus_flank", "mf_lo95", "mf_hi95"]) + "\n")
BS = {}
groups = ["F_GSC_35to70", "M_GSC_35to70", "early_GSC_1to28", "F_Control_35to70", "M_Control_35to70",
          "early_Control_1to28", "F_GSC_35d", "F_GSC_40d", "F_GSC_45d", "F_GSC_50d", "F_GSC_70d"]
print("group | chr4 coding [95%] r | spacer [95%] r | outside [95%] r | coding−outside [95%]")
for g in groups:
    row = {}
    for p in ("코딩부", "스페이서", "바깥"):
        v = d.get((g, "chr4", p), [])
        if len(v) >= 2:
            bs = boot(v)
            BS[(g, p)] = bs
            row[p] = (frac(v), ci(bs), len(v), sum(x[0] for x in v))
        else:
            row[p] = None
    if row["코딩부"] and row["바깥"]:
        diff = sorted(a - b for a, b in zip(BS[(g, "코딩부")], BS[(g, "바깥")]))
        mf = (row["코딩부"][0] - row["바깥"][0], diff[int(0.025 * NB)], diff[int(0.975 * NB)])
    else:
        mf = None
    txt = []
    for p in ("코딩부", "스페이서", "바깥"):
        x = row[p]
        txt.append("%.3f [%.3f,%.3f] r%d" % (x[0], x[1][0], x[1][1], x[2]) if x else "-")
        out.write(TAB.join([g, "chr4", p] + ([str(x[2]), str(x[3]), "%.4f" % x[0], "%.4f" % x[1][0], "%.4f" % x[1][1]] if x else ["0", "0", "-", "-", "-"])
                           + (["%+.4f" % mf[0], "%+.4f" % mf[1], "%+.4f" % mf[2]] if (mf and p == "코딩부") else ["", "", ""])) + "\n")
    print("%s | %s | %s" % (g, " | ".join(txt), ("%+.3f [%+.3f,%+.3f]" % mf) if mf else "-"))

print()
print("between-group difference (chr4 coding, bootstrap 95%) :")
for a, b in [("F_GSC_35to70", "M_GSC_35to70"), ("F_GSC_35to70", "F_Control_35to70"),
             ("F_GSC_35to70", "early_GSC_1to28"), ("M_GSC_35to70", "M_Control_35to70")]:
    for p in ("코딩부",):
        if (a, p) in BS and (b, p) in BS:
            diff = sorted(x - y for x, y in zip(BS[(a, p)], BS[(b, p)]))
            m = frac(d[(a, "chr4", p)]) - frac(d[(b, "chr4", p)])
            ptwo = 2 * min(sum(1 for x in diff if x <= 0), sum(1 for x in diff if x >= 0)) / NB
            print("  %s − %s : %+.3f [%+.3f,%+.3f]  p≈%.3f" % (a, b, m, diff[int(0.025 * NB)], diff[int(0.975 * NB)], min(1, ptwo)))
# chr18
print()
for g in groups:
    for p in ("코딩부", "스페이서", "바깥"):
        v = d.get((g, "chr18", p), [])
        if v:
            print("chr18 %s %s r%d c%d frac %.3f" % (g, p, len(v), sum(x[0] for x in v), frac(v)))
out.close()
