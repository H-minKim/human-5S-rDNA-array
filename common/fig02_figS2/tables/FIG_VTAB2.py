# -*- coding: utf-8 -*-
r"""Version without inherited variants (multi-individual shared type), same 6-panel layout

  A1 absolute, all 23 groups      A2 absolute, only the groups that change most in a zone
  B1 relative, all (sum = 1)      B2 relative, renormalized among the chosen groups (sum = 1)
  Ca ratio to middle, all 23 groups   Cb chosen groups only
  In this version Ca and Cb are also drawn with 23 groups (Ca and Cb of the main figure use the old 11 types).

Layer: ORIGIN_KEEP = notshared (default, multi-individual shared type removed) / single (single-unit type only)
This layer is not 'newly arisen variants'. See the warning in origin_filter.py.
"""
import io, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CODE = __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir))
sys.path.insert(0, CODE)
import v11_io as G
from vgroup import GROUPS, GI, COLOR, vgroup, iscpg_factory
import origin_filter as OF

OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
KEEP = os.environ.get("ORIGIN_KEEP", "notshared")
KTXT = {"notshared": "다개체 공유형(물려받은 것) 뺌", "single": "단일유닛형만",
        "all": "전부"}[KEEP]
NAME = "fig_VTORG_%s_region_type_composition" % KEEP
ZCUT = float(os.environ.get("ZCUT", "5"))
plt.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["svg.fonttype"] = "none"
iscpg = iscpg_factory(G.CAN)

U, _, _ = G.load("primary", verbose=False)
HK = sorted(U)
NH = len(HK)
NG = len(GROUPS)
hapf, kcnt, NS, _ = OF.build(U, HK)
OKC = {"notshared": {"배열내 다중형", "단일유닛형"}, "single": {"단일유닛형"},
       "all": {"다개체 공유형", "배열내 다중형", "단일유닛형"}}[KEEP]

REG = ["P%d" % d for d in range(1, 31)] + ["MID"] + ["D%d" % d for d in range(30, 0, -1)]
RI = {r: i for i, r in enumerate(REG)}
NR = len(REG)
S = np.zeros((NH, NR, NG)); C = np.zeros((NH, NR))
kept = drop = 0
for hi, k in enumerate(HK):
    for row in U[k]:
        ap, ad = row[0], row[1]
        v = np.zeros(NG)
        for x in row[4]:
            if OF.klass(x, k, hapf, kcnt) in OKC:
                v[GI[vgroup(x, iscpg)]] += 1; kept += 1
            else:
                drop += 1
        for r in ([("P%d" % ap)] if ap <= 30 else []) + ([("D%d" % ad)] if ad <= 30 else []) \
                + (["MID"] if (ap > 30 and ad > 30) else []):
            S[hi, RI[r]] += v; C[hi, RI[r]] += 1
Cn = np.where(C > 0, C, np.nan)
MG = np.nanmean(S / Cn[:, :, None], 0)
NHpos = np.sum(C > 0, 0)
TOT = MG.sum(1)
FR = MG / TOT[:, None]
print("층 %s · 남은 변이 출현 %d / %d (%.4f) · 유닛당 %.3f"
      % (KEEP, kept, kept + drop, kept / (kept + drop), TOT[RI["MID"]]))

ZP = os.path.join(HERE, "REPRO", "VTZONE3_%s_zones.tsv" % KEEP)
ZL = [l.rstrip("\n").split("\t") for l in io.open(ZP, encoding="utf8")]
ZH = {c: i for i, c in enumerate(ZL[0])}
ZCOL = ZL[0][1:5]
ZR = {r[0]: dict(zone=[float(x) for x in r[1:5]], z=float(r[ZH["최대_귀무z"]]),
                 far=r[ZH["가장_먼_구역"]]) for r in ZL[1:]}
PICK = sorted([g for g in GROUPS if g in ZR and abs(ZR[g]["z"]) >= ZCUT],
              key=lambda g: -abs(ZR[g]["z"]))
print("고른 무리 %d : %s" % (len(PICK), " · ".join(PICK)))

src = io.open(os.path.join(OUT, NAME + "_source_values.tsv"), "w", encoding="utf8", newline="\n")
src.write("자리\t쪽\t접합부로부터_유닛\t하플로타입수\t유닛당_변이수_합\t"
          + "\t".join("절대_" + g for g in GROUPS) + "\t"
          + "\t".join("몫_" + g for g in GROUPS) + "\t층\t출처\n")
for r in REG:
    i = RI[r]
    src.write("%s\t%s\t%d\t%d\t%.5f\t%s\t%s\t%s\tv11_io.load('primary') · 2단계 평균\n"
              % (r, "MID" if r == "MID" else ("prox" if r[0] == "P" else "dist"),
                 -1 if r == "MID" else int(r[1:]), NHpos[i], TOT[i],
                 "\t".join("%.5f" % x for x in MG[i]), "\t".join("%.5f" % x for x in FR[i]), KEEP))
src.close()

X = np.arange(NR, dtype=float)
MIDX = RI["MID"]
TICK = [("P1", "근위 1"), ("P10", "10"), ("P20", "20"), ("D20", "20"), ("D10", "10"), ("D1", "원위 1")]
TICKC = [("P1", "근위 1"), ("P5", "5"), ("P10", "10"), ("P20", "20"), ("P28", "28"),
         ("D28", "28"), ("D20", "20"), ("D10", "10"), ("D5", "5"), ("D1", "원위 1")]
ZEDGE = [RI["P3"] + 0.5, RI["P17"] + 0.5, RI["D11"] - 0.5]
ZSPAN = [(0, RI["P3"], "근위 1-3"), (RI["P4"], RI["P17"], "근위 4-17"),
         (RI["P18"], RI["D11"], "가운데"), (RI["D10"], NR - 1, "원위 1-10")]
ZFC = {"근위 1-3": "#fdf2e9", "근위 4-17": "#fef9e7", "가운데": "#ffffff", "원위 1-10": "#fdedec"}

fig = plt.figure(figsize=(17.6, 20.4))
gs = fig.add_gridspec(4, 2, hspace=0.33, wspace=0.14, left=0.062, right=0.800,
                      top=0.950, bottom=0.044)


def stack(ax, Y, keys, ylim, ttl, yl, lab=False):
    base = np.zeros(NR)
    for g in GROUPS:
        j = GI[g]
        if g not in keys:
            continue
        top = base + Y[:, j]
        ax.fill_between(X, base, top, color=COLOR[g], lw=0, zorder=3, label=g if lab else None)
        ax.plot(X, top, color="white", lw=0.5, zorder=4)
        base = top
    ax.plot(X, base, color="#2c3e50", lw=1.0, zorder=5)
    for xx in (MIDX - 0.5, MIDX + 0.5):
        ax.axvline(xx, color="#566573", lw=0.8, ls=(0, (4, 3)), zorder=6)
    ax.set_xlim(0, NR - 1); ax.set_ylim(0, ylim)
    ax.set_xticks([RI[a] for a, _ in TICK])
    ax.set_xticklabels([b for _, b in TICK], fontsize=8.2)
    ax.set_ylabel(yl, fontsize=9.5); ax.set_title(ttl, loc="left", fontsize=11.4)
    ax.tick_params(labelsize=8.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


ALL = set(GROUPS); PK = set(PICK)
TOPA = TOT.max() * 1.06
stack(fig.add_subplot(gs[0, 0]), MG, ALL, TOPA, "A1  절대 · 23 무리 전부", "유닛당 변이 수", lab=True)
stack(fig.add_subplot(gs[0, 1]), MG, PK, TOPA,
      "A2  절대 · 구역에서 가장 많이 변하는 %d 무리만 (|귀무 z| ≥ %.0f)" % (len(PICK), ZCUT),
      "유닛당 변이 수")
stack(fig.add_subplot(gs[1, 0]), FR, ALL, 1.0, "B1  상대 · 23 무리 전부 (합 = 1)", "몫")
PKJ = [GI[g] for g in PICK]
PKS = MG[:, PKJ].sum(1)
ax = fig.add_subplot(gs[1, 1])
stack(ax, MG / PKS[:, None], PK, 1.0,
      "B2  상대 · 고른 %d 무리끼리 다시 정규화 (합 = 1)" % len(PICK), "몫")
ax.text(0.5, -0.135, "분모 = 고른 %d 무리의 합 · 이들이 이 층 전체에서 차지하는 몫 %.3f–%.3f"
        % (len(PICK), (PKS / TOT).min(), (PKS / TOT).max()),
        transform=ax.transAxes, ha="center", va="top", fontsize=8.4, color="#566573")

mid = MG[RI["MID"]]
for ci, (which, ttl) in enumerate(
        [("all", "Ca  가운데 대비 비 · 23 무리 전부"),
         ("pick", "Cb  가운데 대비 비 · 고른 %d 무리만" % len(PICK))]):
    ax = fig.add_subplot(gs[2 + ci, :])
    for a, b, nm in ZSPAN:
        ax.axvspan(a - 0.5, b + 0.5, color=ZFC[nm], zorder=0)
        ax.text((a + b) / 2, -2.55, nm, ha="center", va="bottom", fontsize=8.4,
                color="#566573", zorder=8)
    for xx in ZEDGE:
        ax.axvline(xx, color="#b03a2e", lw=1.2, ls=(0, (5, 3)), zorder=7)
    for g in GROUPS:
        j = GI[g]
        if mid[j] <= 0:
            continue
        y = np.log2(np.where(MG[:, j] > 0, MG[:, j], np.nan) / mid[j])
        if which == "all":
            ax.plot(X, y, color=COLOR[g], lw=1.4, zorder=5, label=g)
        elif g in PICK:
            ax.plot(X, y, color=COLOR[g], lw=1.9, zorder=5,
                    label="%s  (z %+.0f · %s)" % (g, ZR[g]["z"], ZR[g]["far"]))
        else:
            ax.plot(X, y, color="#e5e8e8", lw=0.8, zorder=2)
    if which == "pick":
        zi = {z_: k for k, z_ in enumerate(ZCOL)}
        for g in PICK:
            m0 = ZR[g]["zone"][zi["가운데"]]
            for a, b, nm in ZSPAN:
                v = ZR[g]["zone"][zi[nm]]
                if v > 0 and m0 > 0:
                    ax.plot([a - 0.5, b + 0.5], [np.log2(v / m0)] * 2, color=COLOR[g],
                            lw=3.4, alpha=0.40, solid_capstyle="butt", zorder=4)
    ax.axhline(0, color="#2c3e50", lw=1.0, zorder=6)
    ax.set_xlim(0, NR - 1); ax.set_ylim(-2.7, 8.2)
    ax.set_yticks([-2, -1, 0, 1, 2, 3, 4])
    ax.set_yticklabels(["¼", "½", "1", "2", "4", "8", "16"])
    ax.set_xticks([RI[a] for a, _ in TICKC])
    ax.set_xticklabels([b for _, b in TICKC], fontsize=8.6)
    ax.set_ylabel("가운데 대비 비 (log₂ 눈금)", fontsize=10)
    ax.set_title(ttl, loc="left", fontsize=11.4)
    ax.tick_params(labelsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8.0, loc="upper left", frameon=False, ncol=4 if which == "all" else 3)
    if ci == 1:
        ax.set_xlabel("유닛 자리 — 왼쪽이 근위 접합부, 오른쪽이 원위 접합부 "
                      "(가운데는 양쪽 깊이 >30 을 한 칸으로 압축) · 붉은 점선 = S24 구역 경계",
                      fontsize=10)

h, lb = fig.axes[0].get_legend_handles_labels()
fig.legend(h[::-1], lb[::-1], loc="upper left", bbox_to_anchor=(0.812, 0.950), fontsize=8.6,
           frameon=False, title="변이 무리 (한 변이는 한 무리에만)", title_fontsize=9.2)
# The explanation box was removed on 2026-10-07. Its content moved to the Legend and Results (paper figure rule)
# keep only the layer name: single and notshared are nearly the same figure, so an identifier is needed
fig.suptitle("변이 기원 층 : %s" % KTXT, fontsize=13.0, y=0.985)

fig.savefig(os.path.join(OUT, NAME + ".png"), dpi=150)
fig.savefig(os.path.join(OUT, NAME + ".svg"))
print("저장 %s" % os.path.join(OUT, NAME + ".png"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('fig_VTORG_single_region_type_composition_source_values.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
