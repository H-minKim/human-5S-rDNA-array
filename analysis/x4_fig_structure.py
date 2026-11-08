# -*- coding: utf-8 -*-
# 타 생물종 5S array — 유닛 구조 (길이 · 개수 · 변이)
#  메틸 그림과 같은 형식 : 종을 줄로, array 를 점으로.
#  matplotlib 글리프 함정: ASCII 만.
import io, os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["font.family"] = "Malgun Gothic"
rcParams["axes.unicode_minus"] = False
rcParams["figure.dpi"] = 120

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


OUT = _need("RDNA_PREV_OUT_DIR")
LOG = io.open(OUT + "/X4_STRUCT.txt", "w", encoding="utf8", newline="\n")


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


M = pd.read_csv(OUT + "/tables/OTHERSP_UNIT_METH.tsv", sep="\t")
P("=" * 104)
P("유닛 구조 — 유닛 %s · array %d · 종 %d"
  % (format(len(M), ","), M.array_id.nunique(), M["종"].nunique()))
P("=" * 104)
P()

# array 한 줄로 접기 (시료 여럿이면 시료 평균)
AG = M.groupby(["종", "species_key", "array_id"]).agg(
    n_unit=("array_pos", "max"),
    unit_len=("unit_len", "median"),
    unit_len_sd=("unit_len", "std"),
    n_var=("n_var_cons", "median"),
    n_var_strict=("n_var_strict", "median"),
    n_sub=("n_sub", "median"), n_ins=("n_ins", "median"), n_del=("n_del", "median"),
    depth=("depth_median", "median"),
    meth=("meth", "mean"), n_sample=("sample", "nunique")).reset_index()
AG["span_kb"] = AG.n_unit * AG.unit_len / 1e3
AG.to_csv(OUT + "/tables/OTHERSP_ARRAY_STRUCT.tsv", sep="\t", index=False,
          float_format="%.4f")

ORDER = ["제브라피시", "금화조", "애기장대", "시아망", "보르네오오랑",
         "수마트라오랑", "고릴라", "침팬지", "보노보"]
ORDER = [s for s in ORDER if s in set(AG["종"])]
YI = {s: i for i, s in enumerate(ORDER[::-1])}
AG["yi"] = AG["종"].map(YI)

P("  종            array  유닛길이 중앙(범위)       유닛수 중앙(범위)     array kb 중앙   유닛당 변이 중앙")
for s in ORDER:
    g = AG[AG["종"] == s]
    P("  %-12s %4d   %6.0f (%.0f-%.0f)    %5.0f (%d-%d)   %8.1f      %s"
      % (s, len(g), g.unit_len.median(), g.unit_len.min(), g.unit_len.max(),
         g.n_unit.median(), g.n_unit.min(), g.n_unit.max(), g.span_kb.median(),
         ("%.1f" % g.n_var.median()) if g.n_var.notna().any() else " - "))
P()
P("  유닛 길이는 종 안에서 얼마나 일정한가 (array 사이 퍼짐 / 중앙)")
for s in ORDER:
    g = AG[AG["종"] == s]
    if len(g) < 2:
        P("   %-12s array 1 개 — 못 잼" % s)
        continue
    P("   %-12s 중앙 %6.0f bp · array 사이 sd %6.1f · 변동계수 %.3f  (%s)"
      % (s, g.unit_len.median(), g.unit_len.std(),
         g.unit_len.std() / g.unit_len.median(),
         " ".join("%.0f" % v for v in sorted(g.unit_len))))
P()

fig = plt.figure(figsize=(18.0, 10.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.34, wspace=0.26,
                      left=0.085, right=0.978, top=0.905, bottom=0.075)
rng = np.random.RandomState(3)
COL = {"제브라피시": "#2e7d32", "금화조": "#c0392b", "애기장대": "#7e57c2"}


def strip(ax, col, xlab, title, logx=False, cfun=None):
    for r in AG.itertuples(index=False):
        v = getattr(r, col)
        if v != v:
            continue
        y = r.yi + rng.uniform(-0.18, 0.18)
        ax.scatter([v], [y], s=26 + 30 * np.log10(max(r.span_kb, 1)),
                   color=(cfun(r) if cfun else COL.get(r.종, "#1f5fa9")),
                   edgecolor="k", lw=0.45, zorder=3)
    for s, i in YI.items():
        ax.axhline(i, color="#eee", lw=6, zorder=0)
    if logx:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, q: "%g" % v))
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(ORDER[::-1], fontsize=8.5)
    ax.set_xlabel(xlab)
    ax.set_title(title, fontsize=10.5)


ax = fig.add_subplot(gs[0, 0])
strip(ax, "unit_len", "유닛 길이 (bp · 로그)",
      "(a) 유닛 길이 — 금화조/애기장대는 한 값, 제브라피시는 두 유전자족", logx=True)
for v, lab in ((177, "제브라피시 177\n(모계형)"), (1980, "chr18 1,980\n(체세포형)"),
               (368, "금화조 368"), (2240, "사람 2,240")):
    ax.axvline(v, color="#999", ls=":", lw=0.9)
ax.text(2240, len(ORDER) - 0.35, "사람\n2,240", fontsize=7.2, color="#555",
        ha="center")

ax = fig.add_subplot(gs[0, 1])
strip(ax, "n_unit", "array 안 유닛 수 (로그)", "(b) 유닛 개수", logx=True)

ax = fig.add_subplot(gs[0, 2])
strip(ax, "span_kb", "array 길이 (kb · 로그)", "(c) array 총 길이", logx=True)

ax = fig.add_subplot(gs[1, 0])
strip(ax, "n_var", "유닛당 합의 변이 수", "(d) 유닛당 변이 — 종마다 크게 다르다")
P("  [d] 유닛당 변이 (array 중앙값)")
for s in ORDER:
    g = AG[AG["종"] == s].dropna(subset=["n_var"])
    if not len(g):
        continue
    P("   %-12s 변이 %5.1f (엄격판 %5.1f) · 치환 %4.1f 삽입 %4.1f 결실 %4.1f · 깊이 %3.0f"
      % (s, g.n_var.median(), g.n_var_strict.median(), g.n_sub.median(),
         g.n_ins.median(), g.n_del.median(), g.depth.median()))
P()

# (e) 유닛 길이 vs 유닛 수
ax = fig.add_subplot(gs[1, 1])
for s, g in AG.groupby("종"):
    ax.scatter(g.unit_len, g.n_unit, s=34, color=COL.get(s, "#888"),
               edgecolor="k", lw=0.4, label=s)
ax.set_xscale("log"); ax.set_yscale("log")
for a in (ax.xaxis, ax.yaxis):
    a.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, q: "%g" % v))
ax.set_xlabel("유닛 길이 (bp)")
ax.set_ylabel("유닛 수")
ax.legend(fontsize=7.2, frameon=False, ncol=2)
ax.set_title("(e) 짧은 유닛일수록 많이 모이나", fontsize=10.5)
rr = stats.spearmanr(np.log10(AG.unit_len), np.log10(AG.n_unit))
P("  [e] rho(log 유닛길이, log 유닛수) = %+.3f (p=%.3g · array %d)"
  % (rr.statistic, rr.pvalue, len(AG)))
gz = AG[AG["종"] == "제브라피시"]
rz = stats.spearmanr(np.log10(gz.unit_len), np.log10(gz.n_unit))
P("     제브라피시 안에서만 : %+.3f (p=%.3g · array %d)"
  % (rz.statistic, rz.pvalue, len(gz)))

# (f) 변이 vs 메틸
ax = fig.add_subplot(gs[1, 2])
D = AG.dropna(subset=["n_var", "meth"])
for s, g in D.groupby("종"):
    ax.scatter(g.n_var, g.meth, s=34, color=COL.get(s, "#888"), edgecolor="k",
               lw=0.4, label=s)
ax.set_xlabel("유닛당 합의 변이 수")
ax.set_ylabel("array 메틸")
ax.set_ylim(0, 1)
ax.legend(fontsize=7.2, frameon=False, ncol=2)
ax.set_title("(f) 변이가 많은 array 가 메틸이 다른가", fontsize=10.5)
rv = stats.spearmanr(D.n_var, D.meth)
P("  [f] 변이 x 메틸 — 날것으로는 세 보이나 확립되지 않는다")
P("     전체(종 섞음)      rho = %+.3f (p=%.2g · array %d)  <- 종 교락"
  % (rv.statistic, rv.pvalue, len(D)))
D2 = D.copy()
D2["rv"] = D2.n_var - D2.groupby("종").n_var.transform("mean")
D2["rm"] = D2.meth - D2.groupby("종").meth.transform("mean")
rr2 = stats.spearmanr(D2.rv, D2.rm)
P("     종 평균 제거 후     rho = %+.3f (p=%.3g)" % (rr2.statistic, rr2.pvalue))
for s, g in D.groupby("종"):
    if len(g) < 4:
        continue
    r1 = stats.spearmanr(g.n_var, g.meth)
    P("     %-10s n=%2d  rho = %+.3f (p=%.3f)" % (s, len(g), r1.statistic, r1.pvalue))
fz = D[(D.종 == "금화조") & (D.meth > 0.5)]
P("     금화조에서 chr9 저메틸 2 개를 빼면 : n=%d rho=%+.3f p=%.3f — 사라진다"
  % (len(fz), *stats.spearmanr(fz.n_var, fz.meth)[:2]))
gz3 = D[(D.종 == "제브라피시") & (D.n_unit >= 30)]
P("     제브라피시에서 유닛 30 개 이상만 : n=%d rho=%+.3f p=%.3f — 사라진다"
  % (len(gz3), *stats.spearmanr(gz3.n_var, gz3.meth)[:2]))
P("     ⇒ 작은 array 는 변이 중앙값이 불안정하고, 금화조는 chr9 두 개가 끌었다.")
P("       변이-메틸 연관은 확립되지 않았다.")
ax.set_title("(f) 변이 x 메틸 — 종 교락 · 이상치 주도 (확립 안 됨)", fontsize=10.5)
P()

fig.suptitle("타 생물종 5S array 의 유닛 구조 — array %d · 종 %d "
             "(점 크기 = array 총 길이)" % (len(AG), AG["종"].nunique()), fontsize=13)
fig.text(0.085, 0.012,
         "유닛 길이는 종 사이에서 19 배 다르다. 금화조는 array 8 개가 368-372 bp(변동계수 0.005)로 한 값이지만, "
         "제브라피시는 177 bp(모계형 27) · 298 bp(3) · 1,978 bp(체세포형 chr18) 로 갈린다." + chr(10) +
         "시아망 유닛 3,294 bp 는 다른 유인원(2,2xx)과 크게 다르다 — 인계문이 경고한 "
         "'조립본↔read 일치도 0.088' 과 같은 문제일 수 있다. 시아망 유닛 값은 쓰지 말 것.",
         fontsize=8.6, color="#444")
fig.savefig(OUT + "/figures/Fig_othersp_unit_struct.png")
plt.close(fig)
P("  저장 : 10_other_species/figures/Fig_othersp_unit_struct.png · tables/OTHERSP_ARRAY_STRUCT.tsv")
LOG.close()
print("done")
