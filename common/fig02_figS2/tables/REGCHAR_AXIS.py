# -*- coding: utf-8 -*-
r"""Lay out the regional characters along the array axis (proximal 1 ... middle ... distal 1)

User request
  "for each position proximal 1,2,3... distal ...3,2,1 within each array, show how much proximal character
   and how much distal character it has on average", by length stratum.
  The cut is placed where the transition to the middle type ends in every length stratum.

What it does (the character vectors are rebuilt with the same method as REGCHAR.py)
  1) per unit variant feature vector (S5 classifier) -> standardized against middle units
  2) mean z difference of the proximal end (depth<=3) and distal end (depth<=3) relative to the middle = character vector (unit vector)
     at most 3 middle units are drawn per haplotype so that long arrays do not dominate the reference
  3) a proximal score and a distal score per unit, averaged per position (P1..P40, MID, D40..D1)
  4) depth where the transition ends: the first depth at which the position mean enters the 95% confidence band of the middle mean
     and stays there afterwards (per length stratum and per character)
  5) relative multiple by length: per N bin (mean count of that kind / mean count overall). Proportional means it grows with N, fixed means flat
Output: REPRO/REGCHAR_axis.tsv, _transition.tsv, _relcount.tsv
"""
import io, os, sys, collections
import numpy as np

CODE = __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir))
sys.path.insert(0, CODE)
import v11_io as G

OUT = os.path.join((os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"]), "REPRO")
rng = np.random.default_rng(20260927)
NB = 400
DMAX = 40

_src = io.open(os.path.join(CODE, "S5_region_profile.py"), encoding="utf8").read()
_head = _src.split("REG=[")[0]
_head = chr(10).join(l for l in _head.split(chr(10)) if "G.load(" not in l and "os.makedirs" not in l)
ns = {"np": np, "io": io, "os": os, "sys": sys, "collections": collections}
exec(compile(_head, "S5_head", "exec"), ns)
counts, MET = ns["counts"], ns["MET"]

U, A_, M_ = G.load("primary", verbose=False)
HK = sorted(U)
HID = {h: i for i, h in enumerate(HK)}
N = {h: len(U[h]) for h in HK}

X, AP, AD, HI = [], [], [], []
for h in HK:
    for r in U[h]:
        X.append(counts(r[4])); AP.append(r[0]); AD.append(r[1]); HI.append(HID[h])
X = np.asarray(X, float); AP = np.array(AP); AD = np.array(AD); HI = np.array(HI)
MIDM = (AP > 30) & (AD > 30)
mu, sd = X[MIDM].mean(0), X[MIDM].std(0)
Z = (X - mu) / np.where(sd > 0, sd, 1.0)

mid = []
for h in HK:
    c = np.where(MIDM & (HI == HID[h]))[0]
    if len(c):
        mid += list(rng.choice(c, min(3, len(c)), replace=False))
mid = np.array(mid, int)
base = Z[mid].mean(0)
CV = {}
for nm, m in (("근위다움", AP <= 3), ("원위다움", AD <= 3)):
    w = Z[m].mean(0) - base
    CV[nm] = w / (np.linalg.norm(w) or 1.0)
S = {nm: Z @ CV[nm] for nm in CV}
print("유닛 %d · 하플로타입 %d · 유닛당 변이 %.3f"
      % (len(X), len(HK), float(np.mean([len(r[4]) for h in HK for r in U[h]]))))

QS = np.percentile([N[h] for h in HK], [25, 50, 75])
def stratum(n):
    return ("N<=%d" % int(QS[0]) if n <= QS[0] else
            "N %d-%d" % (int(QS[0]) + 1, int(QS[1])) if n <= QS[1] else
            "N %d-%d" % (int(QS[1]) + 1, int(QS[2])) if n <= QS[2] else "N>%d" % int(QS[2]))
ORD = ["N<=%d" % int(QS[0]), "N %d-%d" % (int(QS[0]) + 1, int(QS[1])),
       "N %d-%d" % (int(QS[1]) + 1, int(QS[2])), "N>%d" % int(QS[2])]
UST = np.array([stratum(N[HK[i]]) for i in HI], dtype=object)

LAB = ["P%d" % d for d in range(1, DMAX + 1)] + ["MID"] + ["D%d" % d for d in range(DMAX, 0, -1)]
def sel_of(lab, st):
    m = (UST == st)
    if lab == "MID":
        return m & MIDM
    d = int(lab[1:])
    return m & ((AP == d) if lab[0] == "P" else (AD == d))

# threshold = 95th percentile of middle unit scores (5% pass it by background)
THR = {nm: float(np.percentile(S[nm][mid], 95)) for nm in CV}
g = io.open(os.path.join(OUT, "REGCHAR_axis.tsv"), "w", encoding="utf8", newline="\n")
g.write("특징\t길이층\t자리\t쪽\t깊이\t유닛수\t평균점수\tci_lo\tci_hi\t문턱넘는몫\t초과몫\n")
AXV = {}
for nm in CV:
    for st in ORD:
        for lab in LAB:
            m = sel_of(lab, st)
            k = int(m.sum())
            if k < 10:
                continue
            v = S[nm][m]
            bs = [v[rng.integers(0, k, k)].mean() for _ in range(NB)]
            lo, hi = np.percentile(bs, [2.5, 97.5])
            side = "MID" if lab == "MID" else ("prox" if lab[0] == "P" else "dist")
            dep = 0 if lab == "MID" else int(lab[1:])
            AXV[(nm, st, lab)] = (float(v.mean()), float(lo), float(hi), k)
            _p = float((v >= THR[nm]).mean())
            g.write("%s\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n"
                    % (nm, st, lab, side, dep, k, v.mean(), lo, hi, _p, _p - 0.05))
g.close()

# ---- depth where the transition ends (10% decay criterion) ----
FRAC, HOLD = 0.10, 3
g = io.open(os.path.join(OUT, "REGCHAR_transition.tsv"), "w", encoding="utf8", newline="\n")
g.write("특징\t길이층\t쪽\t깊이1값\t가운데값\t문턱(10%안)\t전환끝깊이\n")
print("\n=== 가운데형 전환이 끝나는 깊이 (깊이1~가운데 거리의 10% 안 · 3 칸 유지) ===")
TRA = collections.defaultdict(dict)
for nm in CV:
    for st in ORD:
        mv = AXV.get((nm, st, "MID"))
        if mv is None:
            continue
        m0 = mv[0]
        for side, pre in (("prox", "P"), ("dist", "D")):
            ds = [d for d in range(1, DMAX + 1) if (nm, st, pre + str(d)) in AXV]
            if not ds:
                continue
            v1 = AXV[(nm, st, pre + str(ds[0]))][0]
            thr = abs(v1 - m0) * FRAC
            end = None
            for i, d in enumerate(ds):
                nxt = ds[i:i + HOLD + 1]
                if len(nxt) < min(HOLD + 1, len(ds) - i) :
                    pass
                if all(abs(AXV[(nm, st, pre + str(dd))][0] - m0) <= thr for dd in nxt):
                    end = d; break
            TRA[nm][(st, side)] = end
            g.write("%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%s\n"
                    % (nm, st, side, v1, m0, thr, "" if end is None else end))
    print("  %s" % nm)
    for side in ("prox", "dist"):
        print("    %-5s %s" % (side, " · ".join("%s %s" % (st, TRA[nm].get((st, side)))
                                                for st in ORD)))
g.close()
alld = [v for nm in CV for v in TRA[nm].values() if v]
CUT = int(max(alld)) if alld else 30
print("\n모든 길이 층·두 특징에서 전환이 끝나는 가장 깊은 자리 = %d" % CUT)
io.open(os.path.join(OUT, "REGCHAR_cut.tsv"), "w", encoding="utf8", newline="\n").write(
    "stat\tvalue\n전환끝_최대깊이\t%d\n기준\t깊이1~가운데 거리의 10%% 안 · %d 칸 유지\n" % (CUT, HOLD))

# ---- relative multiple by length ----
ITEM = {}
for nm, dep in (("근위다움", AP), ("원위다움", AD)):
    ITEM["%s 영역(깊이 10 안)" % nm] = np.array(
        [int(((HI == HID[h]) & (S[nm] >= THR[nm]) & (dep <= 10)).sum()) for h in HK], float)
# the middle changes a lot with the cut depth, so three definitions are reported together
for _c in (30, 10, 5):
    _m = (AP > _c) & (AD > _c)
    ITEM["가운데 유닛(양쪽 깊이>%d)" % _c] = np.array(
        [int(((HI == HID[h]) & _m).sum()) for h in HK], float)
TOK = {"비례 개별변이 2205TC(U84C)": "2205TC", "고정 개별변이 2235TG": "2235TG",
       "고정 개별변이 D1096:1(말단형 진단)": "D1096:1"}
for nm, tk in TOK.items():
    ITEM[nm] = np.array([sum(1 for r in U[h] if tk in r[4]) for h in HK], float)
x = np.array([N[h] for h in HK], float)
EDG = np.percentile(x, [0, 20, 40, 60, 80, 100])
g = io.open(os.path.join(OUT, "REGCHAR_relcount.tsv"), "w", encoding="utf8", newline="\n")
g.write("항목\tN칸\tN평균\tN상대\t개수평균\t개수상대\t하플로수\n")
print("\n=== 길이 칸별 상대 배수 (1 = 전체 평균) ===")
for nm, y in ITEM.items():
    row = []
    for i in range(5):
        m = (x >= EDG[i]) & (x <= EDG[i + 1]) if i == 4 else (x >= EDG[i]) & (x < EDG[i + 1])
        if m.sum() < 10:
            continue
        g.write("%s\t%d\t%.2f\t%.4f\t%.4f\t%.4f\t%d\n"
                % (nm, i + 1, x[m].mean(), x[m].mean() / x.mean(), y[m].mean(),
                   (y[m].mean() / y.mean()) if y.mean() else np.nan, int(m.sum())))
        row.append("%.2f" % ((y[m].mean() / y.mean()) if y.mean() else np.nan))
    print("  %-34s %s" % (nm, " · ".join(row)))
g.close()
print("\n저장 REPRO/REGCHAR_axis.tsv · _transition.tsv · _cut.tsv · _relcount.tsv")

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('REGCHAR_axis.tsv', 'REGCHAR_relcount.tsv'): _taben.to_en_file(_o.path.join(OUT, _f))
