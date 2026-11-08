# -*- coding: utf-8 -*-
r"""Main Fig 2 cell e: only figure 1 of MAJREF.py (fig_T0_majref), at cell size (191.4 x 82 pt)
   Data and computation: lines 20-153 of MAJREF.py moved over without changing a single character (only the value table names became CELL_F2e_*, so the original value tables are not overwritten)
   Like MAJREF, run with hook_tieronly / hook_tieronly_en (with the ordinary hook the majority conversion is applied twice)
   Only the appearance changed: figure size, font size (common to main-text cells), line width, annotation text positions, no figure title or large title (the composite removes them anyway)
   Content check: CMPART.py (old fig_T0_majref axes 0 vs this figure)
Output: $FIGOUT\cellF2_e_majref.png/svg, CELL_F2e_source_values.tsv, CELL_F2e_flips.tsv
"""
import io, os, sys, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import matplotlib.font_manager as _fm
for _p in ((__import__("os").environ.get("RDNA_FONT_DIR", __import__("os").path.join(__import__("os").environ.get("SystemRoot", ""), "Fonts")) + "/malgun.ttf"), (__import__("os").environ.get("RDNA_FONT_DIR", __import__("os").path.join(__import__("os").environ.get("SystemRoot", ""), "Fonts")) + "/gulim.ttc")):
    if os.path.exists(_p):
        _fm.fontManager.addfont(_p)
        plt.rcParams["font.family"] = _fm.FontProperties(fname=_p).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["svg.fonttype"] = "none"
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir, "common", "fig02_figS2")))
import v11_io as G

# ---- element names: code and tables unchanged, only the figure text is changed ----
_DISP = {"AluY": "AluY element", "NTS1": "spacer-1", "NTS2": "spacer-2",
         "GArich": "GA-rich", "GA-rich": "GA-rich",
         "promoter": "5′-flanking", "prom": "5′-flanking", "프로모터": "5′-flanking",
         "gene": "5S gene", "coding": "5S gene"}
def disp(s):
    """Element name of the table or ELEM -> name to write in the figure (stored values are not changed)"""
    return _DISP.get(s, s)


OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
BLOCKS = G.BLOCKS
ELEM = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA-rich", 1179, 1551),
        ("NTS2", 1552, 2031), ("promoter", 2032, 2121), ("gene", 2122, 2240)]
UL = 2240; DMAX = 30

def pos_of(t):
    if t.startswith("DB") or t.startswith("IB"):
        b = int(t[2:].split(":")[0]); return BLOCKS[b - 1][0] if 1 <= b <= len(BLOCKS) else None
    if t[0] == "D": return int(t[1:].split(":")[0])
    if t[0] == "I": return int(t[1:].split(":")[0])
    if t[0] == "R": return int(t[1:-2])
    try: return int(t[:-2])
    except Exception: return None

U, _, _ = G.load("primary", verbose=False)
HK = sorted(U)
# ---- per-position alleles for each unit ----
units = []        # (hap, array_pos, ap, ad, {position: token})
site_alle = collections.defaultdict(collections.Counter)
for k in HK:
    for r in U[k]:
        d = {}
        for t in r[4]:
            p = pos_of(t)
            if p is None or not (1 <= p <= UL): continue
            if p in d: d[p] = min(d[p], t)        # 2 at the same position (4 cases): the first in lexicographic order
            else: d[p] = t
        units.append((k, r[6]["array_pos"], r[0], r[1], d))
        for p, t in d.items(): site_alle[p][t] += 1
NU = len(units)
for p in site_alle: site_alle[p]["REF"] = NU - sum(site_alle[p].values())
# the reference choice was unified from mode -> majority (c/NU > 0.5).
#   hook_tier (= the G6 and lambda measurement frame) uses the same rule. The mode rule takes a non-majority token as the reference,
#   as at position 1476 (29 allele types, mode 0.4233), so the reference shifts with the sample.
MAJ = {}
for p in site_alle:
    _t, _c = site_alle[p].most_common(1)[0]
    MAJ[p] = _t if (_c / NU) > 0.5 else "REF"      # without a majority the CHM13 reference is kept

FLIP = sorted([p for p in MAJ if MAJ[p] != "REF"])
print("유닛 %d · 자리 %d · 다수형이 CHM13 과 다른 자리 %d 개" % (NU, len(site_alle), len(FLIP)))

def el(p):
    for n_, a, b in ELEM:
        if a <= p <= b: return n_
    return "?"
fw = io.open(os.path.join(OUT, "CELL_F2e_flips.tsv"), "w", encoding="utf8", newline="\n")
fw.write("bp\t구획\tCHM13기준일때_변이유닛비율\t다수형토큰\t다수형비율\t다수형기준일때_변이유닛비율\t그자리대립형질수\n")
print("\n기준이 바뀌는 자리")
print("  bp    구획       CHM13 기준 변이율   다수형 토큰        다수형 비율   다수형 기준 변이율")
for p in FLIP:
    old = (NU - site_alle[p]["REF"]) / NU
    m = MAJ[p]; mf = site_alle[p][m] / NU; new = 1 - mf
    print("  %4d  %-9s %16.4f   %-16s %10.4f %16.4f" % (p, el(p), old, m, mf, new))
    fw.write("%d\t%s\t%.6f\t%s\t%.6f\t%.6f\t%d\n" % (p, el(p), old, m, mf, new, len(site_alle[p])))
fw.close()

# ---- number of new variants per unit ----
FS = set(FLIP)
rows = []
for k, apos, ap, ad, d in units:
    old = len(d)
    same = sum(1 for p in FS if d.get(p, "REF") == MAJ[p])
    new = sum(1 for p in d if p not in FS) + (len(FS) - same)
    rows.append((k, apos, ap, ad, old, new, d))
O = np.array([r[4] for r in rows], float); N_ = np.array([r[5] for r in rows], float)
print("\n유닛당 변이 수 : CHM13 기준 %.3f  →  다수형 기준 %.3f  (%.1f%% 감소)"
      % (O.mean(), N_.mean(), 100 * (1 - N_.mean() / O.mean())))

# ---- per-position profile (zone 5) ----
def zone5(ap, ad):
    if ap <= ad:
        return "근위1-10" if ap <= 10 else ("근위11-30" if ap <= 30 else "가운데")
    return "원위1-10" if ad <= 10 else ("원위11-30" if ad <= 30 else "가운데")
Z5 = ["근위1-10", "근위11-30", "가운데", "원위11-30", "원위1-10"]
C5 = {"근위1-10": "#c0392b", "근위11-30": "#e59866", "가운데": "#566573",
      "원위11-30": "#7fb3d5", "원위1-10": "#2471a3"}
SITE = {z: collections.Counter() for z in Z5}; NZ = collections.Counter()
REG = collections.defaultdict(list)
for k, apos, ap, ad, old, new, d in rows:
    z = zone5(ap, ad); NZ[z] += 1
    for p in d:
        if p not in FS: SITE[z][p] += 1       # the majority form at these positions is REF, so a token means a variant
    for p in FS:
        if d.get(p, "REF") != MAJ[p]: SITE[z][p] += 1
    rg = ("P%d" % ap) if ap <= DMAX else (("D%d" % ad) if ad <= DMAX else ("MID" if ap > DMAX and ad > DMAX else None))
    if rg: REG[rg].append(new)
LAB = ["P%d" % d_ for d_ in range(1, DMAX + 1)] + ["MID"] + ["D%d" % d_ for d_ in range(DMAX, 0, -1)]
# the x axis notation is matched to the same form as fig_VTORG2 (main g and h)
_TK = [("P1", "근위 1"), ("P10", "10"), ("P20", "20"), ("D20", "20"), ("D10", "10"), ("D1", "원위 1")]
def _xticks(ax, lab):
    li = {l: i for i, l in enumerate(lab)}
    tk = [(li[a], b) for a, b in _TK if a in li]
    ax.set_xticks([i for i, _ in tk])
    ax.set_xticklabels([b for _, b in tk], fontsize=8.2, rotation=0)
LAB = [l for l in LAB if len(REG.get(l, [])) >= 50]

SRC = io.open(os.path.join(OUT, "CELL_F2e_source_values.tsv"), "w", encoding="utf8", newline="\n")
SRC.write("종류\t키\t값\t비고\n")
SRC.write("요약\t유닛수\t%d\t\n요약\tCHM13기준_유닛당\t%.5f\t\n요약\t다수형기준_유닛당\t%.5f\t\n"
          % (NU, O.mean(), N_.mean()))
for l in LAB:
    v = np.array(REG[l], float)
    SRC.write("array자리\t%s\t%.5f\t유닛 %d\n" % (l, v.mean(), len(v)))
for z in Z5:
    for p in sorted(SITE[z]):
        SRC.write("유닛자리\t%s|%d\t%.6f\t유닛 %d\n" % (z, p, SITE[z][p] / NZ[z], NZ[z]))
SRC.close()

# ================= cell figure (191.4 x 82 pt) =================
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cellstyle as CS
_EN = bool(os.environ.get("EN_DICT"))                         # in the English version the horizontal axis title takes two lines (EN_WRAP_MAP)
fig, ax = CS.cell_at("e")   # margins and size from the single cellstyle.ROWS table (plot frame alignment)
XM = LAB.index("MID")
ax.axvspan(-0.5, XM - 0.5, color="#eaf2f8", zorder=0)
ax.axvspan(XM - 0.5, XM + 0.5, color="#d5d8dc", zorder=0)
ax.axvspan(XM + 0.5, len(LAB) - 0.5, color="#fdedec", zorder=0)
y = np.array([np.mean(REG[l]) for l in LAB], float)
se = np.array([np.std(REG[l], ddof=1) / np.sqrt(len(REG[l])) for l in LAB], float)
ax.fill_between(range(len(LAB)), y - 1.96 * se, y + 1.96 * se, color="#aab7b8", alpha=0.45, lw=0, zorder=2)
ax.plot(range(XM), y[:XM], "-o", ms=1.3, lw=0.7, color="#2471a3", zorder=3, label="근위 쪽")
ax.plot(range(XM, len(LAB)), y[XM:], "-o", ms=1.3, lw=0.7, color="#c0392b", zorder=3, label="원위 쪽")
ax.plot([XM], [y[XM]], "o", ms=3.0, color="#566573", zorder=4, label="가운데")
# only the annotation text positions are fitted to the cell (the points and text they refer to are unchanged), (offset in points, horizontal alignment)
_OFF = {"P1": ((2, -6), "left"), "P2": ((2, 12), "left"), "P3": ((5, 1), "left"),
        "D7": ((-2, 5), "right"), "D5": ((0, 10), "center"), "D2": ((-3, 0), "right"), "D1": ((-3, -15), "right")}   # D1 goes down inside the cell (right frame alignment)
# peak value marks (P1, P2, P3, D7, D5, D2, D1) are removed from the figure (values are in the per-position rows of the value table)
for l in ():
    if l not in LAB: continue
    i = LAB.index(l)
    (_dx, _dy), _ha = _OFF[l]
    ax.annotate("%s %.2f" % (l, y[i]), xy=(i, y[i]), xytext=(_dx, _dy), textcoords="offset points", ha=_ha, va="center",
                fontsize=CS.FS_ANN, weight="bold", annotation_clip=False,
                arrowprops=dict(arrowstyle="-", color="#17202a", lw=0.35, shrinkA=0, shrinkB=0))
ax.set_xlim(-0.5, len(LAB) - 0.5)
_xticks(ax, LAB)
# tick label colour and weight removed (old proximal #2471a3, distal #c0392b, middle #566573 bold), L (cellF2_m_cpg) is already black
for t_ in ax.get_xticklabels(): t_.set_fontsize(CS.FS_TICK)
_yl = ax.get_ylim()                                            # vertical range set by the data, as in the original
ax.set_yticks([6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5])      # original ticks unchanged
ax.set_ylim(_yl)                                               # fixing the ticks widens the range, so restore it
ax.set_ylabel(CS.L("유닛당 변이 수", "Variants per unit"))   # axis title on one line, short (old: number of variants per unit (majority reference))
# the horizontal axis title is removed and goes to the Legend text (the ticks read 'proximal 1 ... distal 1')
ax.legend(loc="upper center", ncol=3, frameon=False)
CS.save(fig, OUT, "cellF2_e_majref")
print("유닛당 %.2f 개 (CHM13 기준 %.2f) · 바뀐 자리 %d" % (N_.mean(), O.mean(), len(FLIP)))
