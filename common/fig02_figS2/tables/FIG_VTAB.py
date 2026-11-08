# -*- coding: utf-8 -*-
r"""Variant composition from proximal 1 -> middle (compressed) -> distal 1, panels A and B each in two versions

  A1  absolute, all 23 groups unfiltered (sum of layers = variants per unit at that position)
  A2  absolute, only the groups that change most between zones (|null z| >= ZCUT)
  B1  relative, all 23 groups (variant count at that position set to 1)
  B2  relative, selected groups renormalized among themselves (their sum = 1, the denominator is not the total variant count)
  Ca  ratio to the middle for the unfiltered 11 types      <- left unchanged
  Cb  only the groups that change most across zones        <- left unchanged

Cohort v11 primary 564, two-stage mean, zone boundaries = S24 (proximal 3|4, 17|18, distal 10|11)
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

OUT = (os.environ.get("FIGOUT") or os.environ["RDNA_OUT_DIR"])
NAME = "fig_VTPD_region_type_composition"
plt.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["svg.fonttype"] = "none"
iscpg = iscpg_factory(G.CAN)
ZCUT = 20.0

U, _, _ = G.load("primary", verbose=False)
HK = sorted(U)
NH = len(HK)
NG = len(GROUPS)

# ---------- coarse 11 types (panel C only, do not change) ----------
CTY = [("1bp 결실", "#c0392b"), ("반복블록 결실", "#e67e22"), ("산화 C>A·G>T", "#f1c40f"),
       ("탈아미노 비CpG", "#2471a3"), ("탈아미노 CpG", "#5dade2"),
       ("그 밖 전이 T>C·A>G", "#117a65"), ("그 밖 전환", "#aab7b8"), ("긴 결실", "#7d3c98"),
       ("1bp 삽입", "#48c9b0"), ("긴 삽입", "#1abc9c"), ("반복블록 삽입", "#85929e")]
CTI = {t: i for i, (t, _) in enumerate(CTY)}


def ctype(x):
    if x.startswith("DB"):
        return "반복블록 결실"
    if x.startswith("IB"):
        return "반복블록 삽입"
    if x[0] == "D":
        return "1bp 결실" if int(x[1:].split(":")[1]) == 1 else "긴 결실"
    if x[0] == "I":
        return "1bp 삽입" if int(x[1:].split(":")[1]) == 1 else "긴 삽입"
    p = int(x[1:-2]) if x[0] == "R" else int(x[:-2])
    s = x[-2] + ">" + x[-1]
    if s in ("C>T", "G>A"):
        return "탈아미노 CpG" if iscpg(p) else "탈아미노 비CpG"
    if s in ("C>A", "G>T"):
        return "산화 C>A·G>T"
    if s in ("T>C", "A>G"):
        return "그 밖 전이 T>C·A>G"
    return "그 밖 전환"


# ---------- two-stage mean per position (both classifications together) ----------
REG = ["P%d" % d for d in range(1, 31)] + ["MID"] + ["D%d" % d for d in range(30, 0, -1)]
RI = {r: i for i, r in enumerate(REG)}
NR = len(REG)
Sg = np.zeros((NH, NR, NG)); Sc = np.zeros((NH, NR, len(CTY))); C = np.zeros((NH, NR))
for hi, k in enumerate(HK):
    for row in U[k]:
        ap, ad = row[0], row[1]
        vg = np.zeros(NG); vc = np.zeros(len(CTY))
        for x in row[4]:
            vg[GI[vgroup(x, iscpg)]] += 1
            vc[CTI[ctype(x)]] += 1
        for r in ([("P%d" % ap)] if ap <= 30 else []) + ([("D%d" % ad)] if ad <= 30 else []) \
                + (["MID"] if (ap > 30 and ad > 30) else []):
            Sg[hi, RI[r]] += vg; Sc[hi, RI[r]] += vc; C[hi, RI[r]] += 1
Cn = np.where(C > 0, C, np.nan)
MG = np.nanmean(Sg / Cn[:, :, None], 0)
MC = np.nanmean(Sc / Cn[:, :, None], 0)
NHpos = np.sum(C > 0, 0)
TOT = MG.sum(1)
FR = MG / TOT[:, None]
print("유닛당 변이 — 두 분류 합이 같은가 : 최대 차 %.3g" % np.nanmax(np.abs(MG.sum(1) - MC.sum(1))))

# ---------- read zone outputs (selection rule) ----------
ZP = os.path.join(HERE, "REPRO", "VTZONE2_zones.tsv")
ZL = [l.rstrip("\n").split("\t") for l in io.open(ZP, encoding="utf8")]
ZH = {c: i for i, c in enumerate(ZL[0])}
ZCOL = ZL[0][1:5]
ZR = {r[0]: dict(zone=[float(x) for x in r[1:5]], z=float(r[ZH["최대_귀무z"]]),
                 far=r[ZH["가장_먼_구역"]]) for r in ZL[1:]}
PICK = [g for g in GROUPS if g in ZR and abs(ZR[g]["z"]) >= ZCUT]
PICK.sort(key=lambda g: -abs(ZR[g]["z"]))
print("고른 무리 %d / %d : %s" % (len(PICK), NG, " · ".join(PICK)))

ZP2 = os.path.join(HERE, "REPRO", "VTZONE_zones.tsv")      # coarse 11 types (for panel Cb)
import taben  # English tables are read back with their original labels
ZL2 = [l.rstrip("\n").split("\t") for l in taben.open_ko(ZP2, encoding="utf8")]
ZH2 = {c: i for i, c in enumerate(ZL2[0])}
ZR2 = {r[0]: dict(zone=[float(x) for x in r[1:5]], z=float(r[ZH2["최대_귀무z"]]),
                  far=r[ZH2["가장_먼_구역"]]) for r in ZL2[1:]}
PICK2 = sorted([t for t, _ in CTY if t in ZR2 and abs(ZR2[t]["z"]) >= 10.0],
               key=lambda t: -abs(ZR2[t]["z"]))

# ---------- source table ----------
src = io.open(os.path.join(OUT, NAME + "_source_values.tsv"), "w", encoding="utf8", newline="\n")
src.write("자리\t쪽\t접합부로부터_유닛\t하플로타입수\t유닛당_변이수_합\t"
          + "\t".join("절대_" + g for g in GROUPS) + "\t"
          + "\t".join("몫_" + g for g in GROUPS) + "\t"
          + "\t".join("거친_" + t for t, _ in CTY) + "\t출처\n")
for r in REG:
    i = RI[r]
    src.write("%s\t%s\t%d\t%d\t%.5f\t%s\t%s\t%s\tv11_io.load('primary') · 2단계 평균\n"
              % (r, "MID" if r == "MID" else ("prox" if r[0] == "P" else "dist"),
                 -1 if r == "MID" else int(r[1:]), NHpos[i], TOT[i],
                 "\t".join("%.5f" % x for x in MG[i]), "\t".join("%.5f" % x for x in FR[i]),
                 "\t".join("%.5f" % x for x in MC[i])))
src.close()

# ---------- figure ----------
X = np.arange(NR, dtype=float)
MIDX = RI["MID"]
TICK = [("P1", "근위 1"), ("P10", "10"), ("P20", "20"), ("D20", "20"), ("D10", "10"),
        ("D1", "원위 1")]
TICKC = [("P1", "근위 1"), ("P5", "5"), ("P10", "10"), ("P20", "20"), ("P28", "28"),
         ("D28", "28"), ("D20", "20"), ("D10", "10"), ("D5", "5"), ("D1", "원위 1")]
ZEDGE = [RI["P3"] + 0.5, RI["P17"] + 0.5, RI["D11"] - 0.5]
ZSPAN = [(0, RI["P3"], "근위 1-3"), (RI["P4"], RI["P17"], "근위 4-17"),
         (RI["P18"], RI["D11"], "가운데"), (RI["D10"], NR - 1, "원위 1-10")]
ZFC = {"근위 1-3": "#fdf2e9", "근위 4-17": "#fef9e7", "가운데": "#ffffff", "원위 1-10": "#fdedec"}

fig = plt.figure(figsize=(17.6, 20.4))
gs = fig.add_gridspec(4, 2, hspace=0.33, wspace=0.14,
                      left=0.062, right=0.800, top=0.957, bottom=0.044,
                      height_ratios=[1.0, 1.0, 1.0, 1.0])


def stack(ax, Y, keys, ylim, ttl, yl, lab=False):
    base = np.zeros(NR)
    for g in GROUPS:
        j = GI[g]
        if g not in keys:
            continue
        top = base + Y[:, j]
        ax.fill_between(X, base, top, color=COLOR[g], lw=0, zorder=3,
                        label=g if lab else None)
        ax.plot(X, top, color="white", lw=0.5, zorder=4)
        base = top
    ax.plot(X, base, color="#2c3e50", lw=1.0, zorder=5)
    for xx in (MIDX - 0.5, MIDX + 0.5):
        ax.axvline(xx, color="#566573", lw=0.8, ls=(0, (4, 3)), zorder=6)
    ax.set_xlim(0, NR - 1); ax.set_ylim(0, ylim)
    ax.set_xticks([RI[a] for a, _ in TICK])
    ax.set_xticklabels([b for _, b in TICK], fontsize=8.2)
    ax.set_ylabel(yl, fontsize=9.5)
    ax.set_title(ttl, loc="left", fontsize=11.4)
    ax.tick_params(labelsize=8.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


ALL = set(GROUPS)
PK = set(PICK)
TOPA = TOT.max() * 1.06
stack(fig.add_subplot(gs[0, 0]), MG, ALL, TOPA,
      "A1  절대 · 거르지 않은 %d 무리 전부" % NG, "유닛당 변이 수", lab=True)
stack(fig.add_subplot(gs[0, 1]), MG, PK, TOPA,
      "A2  절대 · 구역에서 가장 많이 변하는 %d 무리만 (|귀무 z| ≥ %.0f)" % (len(PICK), ZCUT),
      "유닛당 변이 수")
stack(fig.add_subplot(gs[1, 0]), FR, ALL, 1.0,
      "B1  상대 · %d 무리 전부 (합 = 1)" % NG, "몫")
# B2: renormalize among the selected groups only (their sum = 1)
PKJ = [GI[g] for g in PICK]
PKS = MG[:, PKJ].sum(1)
FRP = MG / PKS[:, None]
ax = fig.add_subplot(gs[1, 1])
stack(ax, FRP, PK, 1.0,
      "B2  상대 · 고른 %d 무리끼리 다시 정규화 (그들끼리 합 = 1)" % len(PICK), "몫")
ax.text(0.5, -0.135,
        "분모 = 고른 %d 무리의 합 (전체 변이 수가 아니다) · "
        "이 무리들이 전체 변이에서 차지하는 몫은 %.3f–%.3f"
        % (len(PICK), (PKS / TOT).min(), (PKS / TOT).max()),
        transform=ax.transAxes, ha="center", va="top", fontsize=8.4, color="#566573")

# ---------- Ca and Cb: unchanged ----------
midc = MC[RI["MID"]]
for ci, (which, ttl) in enumerate(
        [("all", "Ca  거르지 않고 — 변이 유형 전부 (11 종, 거친 분류). 가운데 대비 비"),
         ("pick", "Cb  구역(S24 경계)에서 가장 많이 변하는 무리만 — |귀무 z| ≥ 10 (거친 분류)")]):
    ax = fig.add_subplot(gs[2 + ci, :])
    for a, b, nm in ZSPAN:
        ax.axvspan(a - 0.5, b + 0.5, color=ZFC[nm], zorder=0)
        ax.text((a + b) / 2, -1.62, nm, ha="center", va="bottom", fontsize=8.4,
                color="#566573", zorder=8)
    for xx in ZEDGE:
        ax.axvline(xx, color="#b03a2e", lw=1.2, ls=(0, (5, 3)), zorder=7)
    for j_, (t_, c_) in enumerate(CTY):
        if midc[j_] <= 0:
            continue
        y = np.log2(np.where(MC[:, j_] > 0, MC[:, j_], np.nan) / midc[j_])
        if which == "all":
            ax.plot(X, y, color=c_, lw=1.5, zorder=5, label=t_)
        elif t_ in PICK2:
            ax.plot(X, y, color=c_, lw=1.9, zorder=5,
                    label="%s  (z %+.0f · %s)" % (t_, ZR2[t_]["z"], ZR2[t_]["far"]))
        else:
            ax.plot(X, y, color="#e5e8e8", lw=0.8, zorder=2)
    if which == "pick":
        zi = {z_: k for k, z_ in enumerate(ZCOL)}
        for t_, c_ in CTY:
            if t_ not in PICK2:
                continue
            m0 = ZR2[t_]["zone"][zi["가운데"]]
            for a, b, nm in ZSPAN:
                v = ZR2[t_]["zone"][zi[nm]]
                if v > 0 and m0 > 0:
                    ax.plot([a - 0.5, b + 0.5], [np.log2(v / m0)] * 2, color=c_, lw=3.4,
                            alpha=0.40, solid_capstyle="butt", zorder=4)
    ax.axhline(0, color="#2c3e50", lw=1.0, zorder=6)
    ax.set_xlim(0, NR - 1); ax.set_ylim(-1.75, 6.3)
    ax.set_yticks([-1, 0, 1, 2, 3, 4]); ax.set_yticklabels(["½", "1", "2", "4", "8", "16"])
    ax.set_xticks([RI[a] for a, _ in TICKC])
    ax.set_xticklabels([b for _, b in TICKC], fontsize=8.6)
    ax.set_ylabel("가운데 대비 비 (log₂ 눈금)", fontsize=10)
    ax.set_title(ttl, loc="left", fontsize=11.4)
    ax.tick_params(labelsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=8.2, loc="upper left", frameon=False, ncol=3 if which == "all" else 2)
    if ci == 1:
        ax.set_xlabel("유닛 자리 — 왼쪽이 근위 접합부, 오른쪽이 원위 접합부 "
                      "(가운데는 양쪽 깊이 >30 을 한 칸으로 압축) · 붉은 점선 = S24 구역 경계",
                      fontsize=10)

h, lb = fig.axes[0].get_legend_handles_labels()
fig.legend(h[::-1], lb[::-1], loc="upper left", bbox_to_anchor=(0.812, 0.955), fontsize=8.6,
           frameon=False, title="A·B 의 변이 무리 — 한 변이는 한 무리에만\n(치환 12 종을 뭉치지 않는다)",
           title_fontsize=9.2)

TXT = ["A·B 는 %d 무리 · 판 C 는 옛 11 유형 그대로다." % NG,
       # the number is computed, not hard-coded: tokens are 4,051, while 360,904 is the count of variant occurrences
       "전수 조사 : 변이 출현 %s 개 · '기타' 0 개 (거른 것 없음)." % format(int(TOT.sum()), ","),
       "",
       "뭉치면 지워진다 — 옛 '그 밖 전환' 은 z +2.9 로 평평했지만,",
       "   풀어 보면 안에서 서로 반대로 움직이고 있었다 :"]
for g in PICK:
    TXT.append("   %-10s z %+5.1f  (%s)" % (g, ZR[g]["z"], ZR[g]["far"]))
TXT += ["", "귀무 = 하플로타입 안에서 유닛을 섞음 200 벌",
        "하플로타입 %d (자리마다 %d–%d)" % (NH, NHpos.min(), NHpos.max()),
        "2 단계 평균 : 하플로타입 안 유닛당 평균",
        "              → 하플로타입 사이 평균"]
fig.text(0.812, 0.455, "\n".join(TXT), fontsize=8.4, ha="left", va="top", linespacing=1.55)

fig.savefig(os.path.join(OUT, NAME + ".png"), dpi=150)
fig.savefig(os.path.join(OUT, NAME + ".svg"))
print("저장 %s" % os.path.join(OUT, NAME + ".png"))

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('fig_VTPD_region_type_composition_source_values.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
