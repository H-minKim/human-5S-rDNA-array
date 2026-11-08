# -*- coding: utf-8 -*-
r"""
W15 - not the 'terminal-type marker' but whether every kind of unit distinguishable within a haplotype shows a terminal-specific arrangement (v11-fix2, primary)

User question: D1096:1 is not the point. Do units that are distinguishable from the other units within one array (inherited or anything else)
              appear 'spaced out' at the distal end and 'clustered' at the proximal end?
Definitions (judged within one array only, cohort frequency not used)
  - minority unit kind: when the 15 end units are grouped at Jaccard 0.85 (= W8 rule), a group with 2 or more members
    and at most half of that end. The list of positions (depths) of that group is analyzed.
  - minority variant: a variant token carried by 2 or more and at most half of the 15 end units (the same idea on a different scale)
Test: per haplotype, keep the number of positions and redraw the depths at random (1,000 sets) as the null -> compare fractions of spacing 1, 2 and 3
Split: reported separately for terminal-type groups (including D1096:1), non-terminal-type groups and proximal groups.
Output: $OUTDIR\W15_minority_group_gaps.tsv, W15_minority_variant_gaps.tsv
"""
import io, os, sys, collections
import numpy as np
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir)))
import v11_io as G
import termdef as TD   # terminal-type definition (10 or more of 12 diagnostic variants)
OUT = (os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"])
rng = np.random.default_rng(20260924)
NPERM = 1000; TH = 0.85; DMAX = 15; MARK = "D1096:1"
U, A_, M_ = G.load("primary", verbose=False)
HK = sorted(U)
def jac(a, b): return len(a & b) / max(len(a | b), 1)
END = {}
for k in HK:
    p = {}; d = {}
    for r in U[k]:
        if r[0] <= DMAX: p[r[0]] = r[4]
        if r[1] <= DMAX: d[r[1]] = r[4]
    END[k] = {"근위": p, "원위": d}
def cluster(dd):
    ds = sorted(dd); lab = list(range(len(ds)))
    for i in range(len(ds)):
        for j in range(i + 1, len(ds)):
            if jac(dd[ds[i]], dd[ds[j]]) >= TH:
                a, b = lab[i], lab[j]
                if a != b: lab = [a if x == b else x for x in lab]
    grp = collections.defaultdict(list)
    for i, c in enumerate(lab): grp[c].append(ds[i])
    return grp
# ---------- by group ----------
cases = collections.defaultdict(lambda: {"obs": collections.Counter(), "sets": []})
for k in HK:
    for side in ("근위", "원위"):
        dd = END[k][side]
        if len(dd) < 10: continue
        avail = sorted(dd)
        for c, members in cluster(dd).items():
            if not (2 <= len(members) <= len(dd) / 2): continue
            isterm = side == "원위" and all(TD.is_term(dd[m]) for m in members)
            tag = ("원위 말단형" if isterm else "원위 그 밖의 소수형") if side == "원위" else "근위 소수형"
            cs = cases[tag]
            for i in range(len(members) - 1): cs["obs"][members[i + 1] - members[i]] += 1
            cs["sets"].append((len(members), avail))
g = io.open(os.path.join(OUT, "W15_minority_group_gaps.tsv"), "w", encoding="utf8", newline="\n")
g.write("kind\tn_groups\tn_pairs\tgap\tfrac_obs\tfrac_null\tenrich\tp\n")
print("=== 소수 유닛 '그룹' 의 자리 간격 (하플로타입 안 · 끝 15 유닛 · 자카드 %.2f) ===" % TH)
print("  무리                 그룹   짝    간격1 관측/귀무 (배수)      간격2 관측/귀무      간격3 관측/귀무")
for tag in ("원위 말단형", "원위 그 밖의 소수형", "근위 소수형"):
    cs = cases[tag]; tot = sum(cs["obs"].values())
    if tot < 20: continue
    nul = np.zeros((NPERM, DMAX + 1))
    for j in range(NPERM):
        for n, avail in cs["sets"]:
            pick = np.sort(rng.choice(avail, n, replace=False))
            for i in range(n - 1):
                gp = pick[i + 1] - pick[i]
                if gp <= DMAX: nul[j, gp] += 1
    nf = nul / np.maximum(nul.sum(1, keepdims=True), 1)
    out = []
    for gp in (1, 2, 3):
        fo = cs["obs"].get(gp, 0) / tot; fn = float(nf[:, gp].mean())
        p = (np.sum(np.abs(nf[:, gp] - fn) >= abs(fo - fn)) + 1) / (NPERM + 1.0)
        g.write("%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.3f\t%.4g\n" % (tag, len(cs["sets"]), tot, gp, fo, fn, fo / max(fn, 1e-9), p))
        out.append("%.3f/%.3f (%.2f배, p %.3g)" % (fo, fn, fo / max(fn, 1e-9), p))
    print("  %-18s %5d %5d   %s  %s  %s" % (tag, len(cs["sets"]), tot, out[0], out[1], out[2]))
g.close()
# ---------- by variant token ----------
cases2 = collections.defaultdict(lambda: {"obs": collections.Counter(), "sets": []})
for k in HK:
    for side in ("근위", "원위"):
        dd = END[k][side]
        if len(dd) < 10: continue
        avail = sorted(dd)
        cnt = collections.Counter()
        for d in avail: cnt.update(dd[d])
        for t, n in cnt.items():
            if not (2 <= n <= len(avail) / 2): continue
            pos = [d for d in avail if t in dd[d]]
            tag = ("원위 말단형 표지" if t == MARK else "원위 그 밖") if side == "원위" else "근위"
            cs = cases2[tag]
            for i in range(len(pos) - 1): cs["obs"][pos[i + 1] - pos[i]] += 1
            cs["sets"].append((len(pos), avail))
g = io.open(os.path.join(OUT, "W15_minority_variant_gaps.tsv"), "w", encoding="utf8", newline="\n")
g.write("kind\tn_variant_instances\tn_pairs\tgap\tfrac_obs\tfrac_null\tenrich\tp\n")
print("\n=== 같은 것을 '변이 토큰' 잣대로 (끝 15 유닛 중 2 개 이상·절반 이하가 가진 변이) ===")
print("  무리               사례    짝     간격1 관측/귀무 (배수)      간격2               간격3")
for tag in ("원위 말단형 표지", "원위 그 밖", "근위"):
    cs = cases2[tag]; tot = sum(cs["obs"].values())
    if tot < 20: continue
    sets = cs["sets"]
    if len(sets) > 4000: sets = [sets[i] for i in rng.choice(len(sets), 4000, replace=False)]
    nul = np.zeros((200, DMAX + 1))
    for j in range(200):
        for n, avail in sets:
            pick = np.sort(rng.choice(avail, n, replace=False))
            for i in range(n - 1):
                gp = pick[i + 1] - pick[i]
                if gp <= DMAX: nul[j, gp] += 1
    nf = nul / np.maximum(nul.sum(1, keepdims=True), 1)
    out = []
    for gp in (1, 2, 3):
        fo = cs["obs"].get(gp, 0) / tot; fn = float(nf[:, gp].mean())
        p = (np.sum(np.abs(nf[:, gp] - fn) >= abs(fo - fn)) + 1) / 201.0
        g.write("%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.3f\t%.4g\n" % (tag, len(cs["sets"]), tot, gp, fo, fn, fo / max(fn, 1e-9), p))
        out.append("%.3f/%.3f (%.2f배)" % (fo, fn, fo / max(fn, 1e-9)))
    print("  %-16s %6d %6d   %s  %s  %s" % (tag, len(cs["sets"]), tot, out[0], out[1], out[2]))
g.close()
print("\n저장 :", OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('W15_minority_group_gaps.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
