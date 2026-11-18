# -*- coding: utf-8 -*-
#  The dispersed curve of Fig 4G recombined over a set of sites (the same 29 sites as S5)
#    Rule = prof() / group_get() of `$RDNA_WORK_ROOT/code/disp_03_analyze.py`:
#      for each bin (50 bp, 45 bins) the sum of meth over the sum of n (adding sites and halves 0+1); an empty bin where sum of n < 30
#    Input = `dispersed_vs_array_cut0508/DISP_SITE_BINS_ALL.tsv` (cut 0.508)
#    (1) whether the group rows of DISP_PROFILES.tsv are reproduced with the current 28 sites (cover_pattern = five segments)
#    (2) the curve over the same 29 sites as S5 -> DISP_PROFILES_S42SET29_2026-10-31.tsv (same shape; kind=group; name=S42SET29)
#  Usage: python -B f4g_profile_2026-10-31.py <output folder>
import collections
import hashlib
import io
import os
import sys

def _need(_v):
    """The deposited copy has no default paths. Stop when a root variable is missing."""
    _x = os.environ.get(_v)
    if not _x:
        sys.exit("environment variable %s is not set; see README.md, section 'Environment variables'" % _v)
    return _x


OUT = os.path.abspath(sys.argv[1])
os.makedirs(OUT, exist_ok=True)
D = os.path.join(_need("RDNA_WORK_ROOT"), "dispersed_vs_array_cut0508")
K = os.path.join(_need("DISPERSED_ROOT"), "05_paper_bundle_2026-10-09", "input_tables",
                 "DISP_SITES_TABLE_rule2.tsv")
NB, MINN = 45, 30


def md5(p):
    return hashlib.md5(io.open(p, "rb").read()).hexdigest()[:12]


def rd(p):
    L = io.open(p, encoding="utf8").read().rstrip("\n").split("\n")
    H = L[0].split("\t")
    return [dict(zip(H, x.split("\t"))) for x in L[1:]]


BINS = collections.defaultdict(lambda: [0, 0])
for r in rd(os.path.join(D, "DISP_SITE_BINS_ALL.tsv")):
    w = BINS[(r["site"], r["state"], int(r["bin"]))]
    w[0] += int(r["n"])
    w[1] += int(r["meth"])
S58 = {r["site"]: r for r in rd(os.path.join(D, "DISP_SITES_TABLE.tsv"))}
PR = {(r["kind"], r["name"], r["state"]): r for r in rd(os.path.join(D, "DISP_PROFILES.tsv"))}


def prof(sites, st):
    out = []
    for b in range(NB):
        n = sum(BINS.get((s, st, b), [0, 0])[0] for s in sites)
        k = sum(BINS.get((s, st, b), [0, 0])[1] for s in sites)
        out.append((k / float(n) if n >= MINN else None, n, k))
    return out


#  (1) reproduction
S28 = sorted(s for s, r in S58.items() if r["cover_pattern"] == "NTS1+GA+NTS2+prom+gene")
assert len(S28) == 28
worst = 0.0
for st in ("low", "high"):
    p = prof(S28, st)
    row = PR[("group", "NTS1+GA+NTS2+prom+gene", st)]
    for b in range(NB):
        a, t = p[b][0], row["b%d" % b]
        if a is None or t == "":
            assert a is None and t == "", (st, b, a, t)
            continue
        worst = max(worst, abs(round(a, 4) - float(t)))
print("(1) reproduction of the 28 positions: max |difference| per bin (after rounding to 4 decimals) %.4f" % worst)
assert worst == 0.0

#  (2) the 29 sites of S5
B = [r for r in rd(K) if r["site"] != "S0206" and len(r["cover_pattern"].split("+")) >= 3]
S29 = sorted(r["site"] for r in B)
assert len(S29) == 29 and set(S29) <= set(S58)
for s in S29:
    assert S58[s]["cover_pattern"] == [r for r in B if r["site"] == s][0]["cover_pattern"], s
print("(2) 29 positions, locus total of the 58 panels %d, five segments %d, overlap with 28 %d, only in 28 %s, only in 29 %s" % (
    sum(int(S58[s]["n_loci"]) for s in S29), sum(1 for s in S29 if S58[s]["cover_pattern"].count("+") == 4),
    len(set(S28) & set(S29)), sorted(set(S28) - set(S29)), sorted(set(S29) - set(S28))))
rows = ["kind\tname\tstate\t" + "\t".join("b%d" % b for b in range(NB))]
det = ["set\tstate\tbin\tbin_mid\tn\tmeth\tbeta"]
for nm, ss in (("S28CUR", S28), ("S42SET29", S29)):
    for st in ("low", "high"):
        p = prof(ss, st)
        if nm == "S42SET29":
            rows.append("group\t%s\t%s\t%s" % (nm, st, "\t".join("%.4f" % v if v is not None else "" for v, _n, _k in p)))
        for b, (v, n, k) in enumerate(p):
            det.append("%s\t%s\t%d\t%d\t%d\t%d\t%s" % (nm, st, b, b * 50 + 25, n, k, "%.6f" % v if v is not None else ""))
io.open(os.path.join(OUT, "DISP_PROFILES_S42SET29_2026-10-31.tsv"), "w", encoding="utf8", newline="\n").write("\n".join(rows) + "\n")
io.open(os.path.join(OUT, "DISP_PROFILE_BINS_28vs29_2026-10-31.tsv"), "w", encoding="utf8", newline="\n").write("\n".join(det) + "\n")
io.open(os.path.join(OUT, "S42SET29_sites.tsv"), "w", encoding="utf8", newline="\n").write(
    "site\tchr\tstart\tcover_pattern\tn_loci_58판\n" + "\n".join(
        "%s\t%s\t%s\t%s\t%s" % (s, S58[s]["chr"], S58[s]["start"], S58[s]["cover_pattern"], S58[s]["n_loci"]) for s in S29) + "\n")
#  Summary: the mean per segment (mean over bins; 28 against 29)
EL = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA", 1179, 1551), ("NTS2", 1552, 2031), ("prom", 2032, 2121), ("gene", 2122, 2240)]
for st in ("low", "high"):
    a, b_ = prof(S28, st), prof(S29, st)
    out = []
    for nm, x0, x1 in EL:
        ib = [b for b in range(NB) if x0 <= b * 50 + 25 <= x1]
        va = [a[b][0] for b in ib if a[b][0] is not None]
        vb = [b_[b][0] for b in ib if b_[b][0] is not None]
        out.append("%s %s→%s" % (nm, "%.3f" % (sum(va) / len(va)) if va else "-", "%.3f" % (sum(vb) / len(vb)) if vb else "-"))
    diff = [abs(a[b][0] - b_[b][0]) for b in range(NB) if a[b][0] is not None and b_[b][0] is not None]
    print("  %-4s 28->29: %s, max bin |difference| %.3f, empty bins %d->%d" % (
        st, " · ".join(out), max(diff), sum(v is None for v, _n, _k in a), sum(v is None for v, _n, _k in b_)))
print("Σn (low+high) 28 %d · 29 %d" % (sum(p[1] for st in ("low", "high") for p in prof(S28, st)),
                                       sum(p[1] for st in ("low", "high") for p in prof(S29, st))))
print("input md5 BINS %s, SITES %s, PROFILES %s, rule2 %s" % (
    md5(os.path.join(D, "DISP_SITE_BINS_ALL.tsv")), md5(os.path.join(D, "DISP_SITES_TABLE.tsv")),
    md5(os.path.join(D, "DISP_PROFILES.tsv")), md5(K)))
