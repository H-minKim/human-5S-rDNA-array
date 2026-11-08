# -*- coding: utf-8 -*-
# 대주제 A 논문용 그림 ② — 종마다 5S array 가 염색체 어디에 있나 + 반복 유닛 모식도(길이와 유전자만)
#   자료(전부 읽기만 · 타종 자료 루트) :
#     · 염색체(콘티그) 길이·array 좌표 : 00_catalog/5S_array_structure.tsv (contig_len · array_start/end · period_med · gene_len_med)
#     · 전체 array 목록(분석에 안 쓴 것 포함) : 00_catalog/5S_array_catalog.tsv
#     · 염색체 이름 : 제브라 00_catalog/ZEBRAFISH_ALL_ARRAYS.tsv(chrom_name) · 금화조 06_aux/zfinch_5S_array_meth_summary.tsv(chrom_bigwig)
#       유인원·애기장대는 이름표 자료가 없어 접근번호 그대로 적는다(추정하지 않음).
# · arrays used = tables/FIGDATA_A_array_summary.tsv (53, siamang excluded)
#   그리지 않는 것 : 중심절(센트로미어) — 이 자료에 위치가 없다.
#   산출 : figures/FigA1_array_chromosome_map.{png,svg,pdf} · tables/FIGDATA_A1_chrom_map.tsv
import io, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec
from matplotlib.patches import FancyBboxPatch, Rectangle

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


NEW = _need("RDNA_OUT_DIR")
B = _need("RDNA_DATA_ROOT")
SVGDIR = _need("RDNA_SVG_DIR")
sys.path.insert(0, NEW + "/code")
import arch  # noqa
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
rcParams["font.family"] = "Malgun Gothic"
rcParams["axes.unicode_minus"] = False
rcParams["figure.dpi"] = 120
rcParams["svg.fonttype"] = "none"
rcParams["pdf.fonttype"] = 42
LOG = io.open(NEW + "/logs/X29_CHROMMAP.txt", "w", encoding="utf8", newline="\n")
REASON = "대주제 A — array 염색체 위치 지도 새로 만듦(2026-10-01)"


def P(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


ST = pd.read_csv(B + "/00_catalog/5S_array_structure.tsv", sep="\t")
CA = pd.read_csv(B + "/00_catalog/5S_array_catalog.tsv", sep="\t")
ZF = pd.read_csv(B + "/00_catalog/ZEBRAFISH_ALL_ARRAYS.tsv", sep="\t")
ZNAME = dict(zip(ZF.array_id, ZF.chrom_name))
FN = pd.read_csv(B + "/06_aux/zfinch_5S_array_meth_summary.tsv", sep="\t")
FNAME = {(t.chrom_genbank, int(t.array_start), int(t.array_end)): t.chrom_bigwig for t in FN.itertuples(index=False)}
USED = set(pd.read_csv(NEW + "/tables/FIGDATA_A_array_summary.tsv", sep="\t").array_id)

SP = [("고릴라", "gorilla_mGorGor1"), ("보노보", "bonobo_mPanPan1"), ("침팬지", "chimp_mPanTro3"),
      ("보르네오오랑우탄", "orang_bornean_mPonPyg2"), ("수마트라오랑우탄", "orang_sumatran_mPonAbe1"),
      ("애기장대", "arabidopsis_ColCC"), ("금화조", "zfinch_v03_dip"), ("제브라피시", "zebrafish_GRCz12tu")]
GRP = {"고릴라": "유인원 5 종 — array 안에서 유닛이 이봉", "보노보": "", "침팬지": "", "보르네오오랑우탄": "", "수마트라오랑우탄": "",
       "애기장대": "애기장대 — 전부 고메틸", "금화조": "금화조·제브라피시 — array 수준에서 갈림", "제브라피시": ""}


def cname(aid, acc):
    if aid in ZNAME and isinstance(ZNAME[aid], str):
        return ZNAME[aid]
    tail = aid.split("__")[1]
    a, s_, e_ = tail.rsplit("_", 2)
    k = (a, int(s_), int(e_))
    if k in FNAME:
        return str(FNAME[k])
    return acc                      # 이름 자료 없음


rows = []
for ko, sp in SP:
    g = ST[ST.species == sp]
    for t in g.itertuples(index=False):
        rows.append(dict(종=ko, 종키=sp, array_id=t.array_id, 염색체=t.chrom, 염색체이름=cname(t.array_id, t.chrom),
                         염색체길이=int(t.contig_len), 시작=int(t.array_start), 끝=int(t.array_end),
                         사본수=int(t.copies), 유닛길이=int(t.period_med), 유전자길이=int(t.gene_len_med),
                         분석에씀=("Y" if t.array_id in USED else "N")))
M = pd.DataFrame(rows)
P("array %d (분석 %d · 미분석 %d) · 종 %d" % (len(M), (M.분석에씀 == "Y").sum(), (M.분석에씀 == "N").sum(), M.종.nunique()))
for ko in M.종.unique():
    g = M[M.종 == ko]
    P("  %-10s 염색체 %d · array %d(분석 %d) · 유닛 %s bp · 유전자 %s bp · 가장 긴 염색체 %.1f Mb"
      % (ko, g.염색체.nunique(), len(g), (g.분석에씀 == "Y").sum(),
         "/".join(map(str, sorted(g.유닛길이.unique()))), "/".join(map(str, sorted(g.유전자길이.unique()))),
         g.염색체길이.max() / 1e6))

TMP = os.environ.get("TEMP", ".") + "/x29"
os.makedirs(TMP, exist_ok=True)
M.to_csv(TMP + "/FIGDATA_A1_chrom_map.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")

# ---------------------------------------------------------------- 그림
NCH = {ko: M[M.종 == ko].염색체.nunique() for ko, _ in SP}
H = [max(0.85, NCH[ko] * 0.46) for ko, _ in SP]
fig = plt.figure(figsize=(15.5, 1.05 * sum(H) + 2.6))
gs = gridspec.GridSpec(len(SP), 2, figure=fig, width_ratios=[3.1, 1.0], height_ratios=H,
                       hspace=1.05, wspace=0.14, top=0.915, bottom=0.07, left=0.085, right=0.975)
MBMAX = M.염색체길이.max() / 1e6
UBMAX = M.유닛길이.max()
CG = {"유인원": "#B03A2E", "애기장대": "#1E8449", "금화조": "#B7950B", "제브라피시": "#2471A3"}
for i, (ko, sp) in enumerate(SP):
    g = M[M.종 == ko]
    col = CG.get(ko, CG["유인원"])
    ax = fig.add_subplot(gs[i, 0])
    chroms = sorted(g.염색체.unique(), key=lambda c: -g[g.염색체 == c].염색체길이.iloc[0])
    for j, ch in enumerate(chroms):
        gg = g[g.염색체 == ch]
        L = gg.염색체길이.iloc[0] / 1e6
        y = -j
        ax.add_patch(FancyBboxPatch((0, y - 0.16), L, 0.32, boxstyle="round,pad=0,rounding_size=0.25",
                                    fc="#ECEFF1", ec="#90A4AE", lw=0.8, zorder=2))
        for t in gg.itertuples(index=False):
            x = (t.시작 + t.끝) / 2e6
            used = t.분석에씀 == "Y"
            ax.plot([x, x], [y - 0.30, y + 0.30], "-", color=(col if used else "#B0BEC5"),
                    lw=2.0 if used else 1.1, zorder=4, solid_capstyle="butt")
        ax.text(-0.012 * MBMAX, y, "%s" % gg.염색체이름.iloc[0], fontsize=8.0, ha="right", va="center")
        u = gg[gg.분석에씀 == "Y"]
        ax.text(L + 0.010 * MBMAX, y, "%.0f Mb  ·  array %d%s"
                % (L, len(gg), ("" if not len(u) else " · 최대 %s 사본" % format(int(u.사본수.max()), ","))),
                fontsize=7.2, ha="left", va="center", color="#777")
    ax.set_xlim(-0.14 * MBMAX, MBMAX * 1.10)
    ax.set_ylim(-len(chroms) + 0.35, 0.8)
    ax.axis("off")
    nused, ntot = (g.분석에씀 == "Y").sum(), len(g)
    ttl = "%s  —  array %d 개 (분석 %d · 회색 %d 개는 이 분석에 안 씀) · 염색체 %d" % (ko, ntot, nused, ntot - nused, len(chroms))
    if GRP.get(ko):
        ax.set_title("%s\n%s" % (GRP[ko], ttl), fontsize=9.6, loc="left", color="#222", pad=6)
        ax.title.set_fontweight("bold")
    else:
        ax.set_title(ttl, fontsize=9.6, loc="left", fontweight="bold", color="#222", pad=6)
    if i == len(SP) - 1:                       # 눈금자
        ax.plot([0, 50], [-len(chroms) + 0.1] * 2, "-", color="#333", lw=1.4)
        ax.text(25, -len(chroms) + 0.22, "50 Mb", fontsize=7.6, ha="center", color="#333")

    # 유닛 모식도 — 길이와 유전자만
    ax2 = fig.add_subplot(gs[i, 1])
    gu = g[g.분석에씀 == "Y"]                      # 모식도는 분석에 쓴 array 의 유닛만
    cnt = gu.유닛길이.value_counts()
    dom = int(cnt.index[0])
    ul = [dom]
    far = sorted(set(gu.유닛길이), key=lambda v: -max(v / dom, dom / v))
    if far and max(far[0] / dom, dom / far[0]) >= 1.2:   # 다른 유형은 길이가 1.2 배 넘게 다를 때만 따로 그린다
        ul.append(int(far[0]))
    rng = (int(gu.유닛길이.min()), int(gu.유닛길이.max()))
    for j, u in enumerate(ul[:2]):
        gl = int(g[g.유닛길이 == u].유전자길이.median())
        y = -j * 1.0
        ax2.add_patch(Rectangle((0, y - 0.18), u, 0.36, fc="#FFFFFF", ec="#546E7A", lw=1.0, zorder=2))
        ax2.add_patch(Rectangle((u - gl, y - 0.18), gl, 0.36, fc="#E57373", ec="#B03A2E", lw=1.0, zorder=3))
        narr = int((gu.유닛길이 == u).sum())
        ax2.text(0, y + 0.32, "유닛 %s bp  ·  array %d" % (format(u, ","), narr), fontsize=7.8, ha="left", color="#333")
        ax2.text(u + 0.03 * UBMAX, y, "5S 유전자 %d bp" % gl, fontsize=7.2, ha="left", va="center", color="#B03A2E")
    if rng[0] != rng[1] and len(ul) == 1:
        ax2.text(0, -0.52, "array 마다 %s~%s bp" % (format(rng[0], ","), format(rng[1], ",")),
                 fontsize=7.0, ha="left", color="#777")
    ax2.set_xlim(-0.04 * UBMAX, UBMAX * 1.55)
    ax2.set_ylim(-1.0 * min(len(ul), 2) - 0.35, 0.78)
    ax2.axis("off")
    if i == 0:
        ax2.text(0, 1.02, "반복 유닛 (길이와 5S 유전자만 · 모든 종 같은 bp 눈금)", fontsize=9.4, fontweight="bold",
                 transform=ax2.transAxes, va="bottom")
    if i == len(SP) - 1:
        ax2.plot([0, 500], [-1.0 * min(len(ul), 2) - 0.18] * 2, "-", color="#333", lw=1.4)
        ax2.text(250, -1.0 * min(len(ul), 2) - 0.10, "500 bp", fontsize=7.6, ha="center", color="#333")

fig.suptitle("대주제 A · 종마다 5S array 의 염색체 위치와 반복 유닛 (시아망 제외 · 2026-10-01)", fontsize=13.5, y=0.985)
fig.text(0.085, 0.016,
         "염색체 막대 = 그 염색체(콘티그) 전체 길이 · 세로줄 = 5S array 위치(가운데 좌표) · 색 = 이 분석에 쓴 array · 회색 = 목록에는 있으나 안 쓴 array\n"
         "중심절(센트로미어)은 이 자료에 위치가 없어 그리지 않았다 · 염색체 이름은 제브라피시·금화조만 자료에 있고, 유인원·애기장대는 접근번호 그대로다(추정하지 않음)\n"
         "자료 : 00_catalog/5S_array_structure.tsv · 5S_array_catalog.tsv · ZEBRAFISH_ALL_ARRAYS.tsv · 06_aux/zfinch_5S_array_meth_summary.tsv · 찍은 값 : tables/FIGDATA_A1_chrom_map.tsv · 코드 : code/x29_figA_chrom_map.py",
         fontsize=7.8, color="#444")
base = "FigA1_array_chromosome_map"
for ext in ("png", "svg", "pdf"):
    fig.savefig(TMP + "/%s.%s" % (base, ext))
plt.close(fig)
arch.put(TMP + "/FIGDATA_A1_chrom_map.tsv", "tables/FIGDATA_A1_chrom_map.tsv", REASON)
for ext in ("png", "svg", "pdf"):
    arch.put(TMP + "/%s.%s" % (base, ext), "figures/%s.%s" % (base, ext), REASON)
if os.path.isdir(SVGDIR):
    import shutil
    shutil.copy2(NEW + "/figures/%s.svg" % base, SVGDIR + "/%s.svg" % base)
P("")
P("그림 : latest/figures/%s.{png,svg,pdf}" % base)
LOG.close()
