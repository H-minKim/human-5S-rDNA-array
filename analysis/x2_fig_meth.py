# -*- coding: utf-8 -*-
# 타 생물종 5S array 메틸 — 그림
#
# Reading rules
#   · 종간 절대값 비교 금지 — 시료·호출기 판본이 다르다 (보노보 내부대조 0.13 차)
#   · valid axis 는 'array - 제 flank' 다. 같은 시료·같은 호출기 안의 대비.
#   · 애기장대는 flank 가 bedmethyl 범위 밖이라 Δ 를 못 낸다.
#   · 붉은털은 read 없이 공개 자리별 메틸(피복 1-7%)이라 다른 자료종이다.
#
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
LOG = io.open(OUT + "/X2_FIG.txt", "w", encoding="utf8", newline="\n")


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


T = pd.read_csv(OUT + "/tables/OTHERSP_ARRAY_METH.tsv", sep="\t")
P("=" * 100)
P("타 생물종 5S array 메틸 — array %d · 종 %d" % (len(T), T.species_key.nunique()))
P("=" * 100)
P()

# 종 순서 : 자료질(Tier) -> 계통
ORDER = ["제브라피시", "금화조", "애기장대", "생쥐", "붉은털원숭이",
         "시아망", "보르네오오랑", "수마트라오랑", "고릴라", "침팬지", "보노보"]
ORDER = [s for s in ORDER if s in set(T["종"])]
YI = {s: i for i, s in enumerate(ORDER[::-1])}
T["yi"] = T["종"].map(YI)

# 비교 가능성 딱지
CMP = {"침팬지": "guppy6.3.8", "고릴라": "guppy6.3.8", "보르네오오랑": "guppy6.3.8",
       "수마트라오랑": "guppy6.3.8", "보노보": "6.3.7+6.3.8 혼합", "시아망": "dorado0.9.1",
       "제브라피시": "guppy6.4.2", "금화조": "dorado2.1.2", "애기장대": "dorado2.1.2",
       "생쥐": "dorado v5.2r1", "붉은털원숭이": "read 없음"}

P("  저메틸 array (meth < 0.35) — 전부 나열")
lowT = T[T.meth < 0.35].sort_values("meth")
P("   종            염색체        길이(bp)   유닛     메틸     flank    Δ")
for r in lowT.itertuples(index=False):
    P("   %-12s %-12s %8s   %5s   %.3f   %s   %s"
      % (r.종, str(r.chrom), format(int(r.array_len), ","),
         ("%d" % r.n_unit) if r.n_unit == r.n_unit else "-",
         r.meth,
         ("%.3f" % r.flank) if r.flank == r.flank else "  -  ",
         ("%+.3f" % r.delta) if r.delta == r.delta else "  -   "))
P()

fig = plt.figure(figsize=(18.0, 11.2))
gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1.0], height_ratios=[1.35, 1.0],
                      hspace=0.30, wspace=0.20, left=0.105, right=0.975,
                      top=0.915, bottom=0.075)

# ---------------------------------------------------------------- (a) 종별 array
ax = fig.add_subplot(gs[:, 0])
rng = np.random.RandomState(7)
for r in T.itertuples(index=False):
    y = r.yi + rng.uniform(-0.17, 0.17)
    if r.flank == r.flank:
        ax.plot([r.flank, r.meth], [y, y], "-", color="#bbb", lw=0.9, zorder=1)
        ax.plot([r.flank], [y], "D", color="#888", ms=3.6, zorder=2)
    sz = 22 + 34 * np.log10(max(r.array_len, 1e3) / 1e3)
    c = "#c0392b" if r.meth >= 0.6 else ("#f0a202" if r.meth >= 0.35 else "#1f5fa9")
    ax.scatter([r.meth], [y], s=sz, color=c, edgecolor="k", lw=0.5, zorder=3)
for s, i in YI.items():
    ax.axhline(i, color="#eee", lw=6, zorder=0)
ax.axvline(0.35, color="#1f5fa9", ls=":", lw=1.2)
ax.axvline(0.6, color="#c0392b", ls=":", lw=1.2)
ax.set_yticks(range(len(ORDER)))
ax.set_yticklabels([("%s\n(%s · array %d)"
                     % (s, CMP.get(s, ""), int((T["종"] == s).sum())))
                    for s in ORDER[::-1]], fontsize=8.5)
ax.set_xlim(0, 1.0)
ax.set_xlabel("array 의 CpG 메틸 (5mC)")
ax.set_title("(a) 종별 · array 별 메틸 — 점 = array (크기 = array 길이) · 회색 마름모 = 제 flank\n"
             "종끼리 절대값을 비교하지 말 것 (시료·호출기가 다르다). 같은 줄 안에서만 읽는다.",
             fontsize=10.5)
#  이름표 : 저메틸 array + flank 보다 크게 낮은 array (Delta < -0.10) 도 같이
labT = T[(T.meth < 0.35) | (T.delta < -0.10)].copy()
labT["lab"] = np.where(labT.chrom_name.astype(str) != "nan",
                       labT.chrom_name.astype(str), labT.chrom.astype(str))
used = {}
for r in labT.sort_values("meth").itertuples(index=False):
    k = r.yi
    used[k] = used.get(k, 0) + 1
    dy = 9 if used[k] % 2 else -14
    ax.annotate("%s (%s%s)" % (r.lab, ("%+.2f" % r.delta) if r.delta == r.delta else "Δ?",
                               ""),
                (r.meth, r.yi), textcoords="offset points",
                xytext=(8, dy), fontsize=7.2, color="#12457f",
                arrowprops=dict(arrowstyle="-", color="#9bb8d6", lw=0.6))

# ---------------------------------------------------------------- (b) Δ
ax = fig.add_subplot(gs[0, 1])
D = T.dropna(subset=["delta"])
for r in D.itertuples(index=False):
    y = r.yi + rng.uniform(-0.17, 0.17)
    ax.scatter([r.delta], [y], s=22 + 34 * np.log10(max(r.array_len, 1e3) / 1e3),
               color=("#c0392b" if r.delta > 0 else "#1f5fa9"),
               edgecolor="k", lw=0.5)
ax.axvline(0, color="#333", lw=1.2)
ax.set_yticks(range(len(ORDER)))
ax.set_yticklabels(ORDER[::-1], fontsize=8.5)
ax.set_xlabel("array - 제 flank (100 kb)")
ax.set_title("(b) 진짜 비교축 — 같은 시료·같은 호출기 안의 대비\n"
             "(애기장대는 flank 자료가 없어 빠짐)", fontsize=10.5)
P("  [b] 종별 Δ(array - flank)")
for s in ORDER:
    g = D[D["종"] == s]
    if not len(g):
        P("   %-12s (Δ 없음)" % s)
        continue
    P("   %-12s array %2d · Δ 중앙 %+.3f · 범위 %+.3f ~ %+.3f"
      % (s, len(g), g.delta.median(), g.delta.min(), g.delta.max()))
P()

# ---------------------------------------------------------------- (c) 크기 vs 메틸
ax = fig.add_subplot(gs[1, 1])
P("  [c] array 크기와 메틸 — 종 안에서 상관 (array 3 개 이상인 종만)")
COL = {"제브라피시": "#2e7d32", "금화조": "#c0392b", "애기장대": "#7e57c2",
       "붉은털원숭이": "#ef6c00"}
for s, g in T.groupby("종"):
    if len(g) < 3:
        continue
    ax.scatter(g.array_len, g.meth, s=30, color=COL.get(s, "#888"),
               edgecolor="k", lw=0.4, label="%s (n=%d)" % (s, len(g)))
    rr = stats.spearmanr(np.log10(g.array_len), g.meth, nan_policy="omit")
    P("   %-12s n=%2d · rho(log 길이, 메틸) = %+.3f (p=%.3f)"
      % (s, len(g), rr.statistic, rr.pvalue))
ax.set_xscale("log")
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
    lambda v, q: ("%gk" % (v / 1e3)) if v >= 1e3 else "%g" % v))
ax.set_xlabel("array 길이 (bp · 로그)")
ax.set_ylabel("array 메틸", fontsize=9)
ax.set_ylim(0, 1.0)
ax.legend(fontsize=7.8, frameon=False)
ax.set_title("(c) array 가 클수록 메틸이 높은가 — 종 안에서", fontsize=10.5)
P()

fig.suptitle("타 생물종 5S rDNA array 의 메틸 — array %d 개 · 종 %d "
             "(자료 $ANALYSIS_ROOT/2026-09-11_5S_other_species)"
             % (len(T), T.species_key.nunique()), fontsize=13)
fig.text(0.105, 0.012,
         "파랑 = 메틸 0.35 미만(꺼진 array) · 주황 = 0.35-0.6 · 빨강 = 0.6 이상.  "
         "점 크기는 array 길이(로그).\n"
         "붉은털원숭이는 read 없이 공개 자리별 메틸(자리 피복 1-7%)이라 다른 자료종이다. "
         "생쥐만 뇌 조직, 유인원은 세포주, 제브라피시는 성체 뇌다.", fontsize=8.6, color="#444")
fig.savefig(OUT + "/figures/Fig_othersp_array_meth.png")
plt.close(fig)
P("  저장 : 10_other_species/figures/Fig_othersp_array_meth.png")
LOG.close()
print("done")
