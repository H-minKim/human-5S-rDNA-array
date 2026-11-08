# -*- coding: utf-8 -*-
# read 단위 : 각 read 의 각 유닛 메틸(β = 그 read 안 그 유닛의 n_5mC / n_calls) 분포를
#   종 · array 마다 따로 그리고, 분포 자체에서 봉우리가 하나인지 둘인지 판정한다.
#
#   자료 : $ANALYSIS_ROOT/2026-09-11_5S_other_species/data_v2/03_reads/<species>/*.read_unit.tsv
#          (한 행 = read 하나 x 유닛 하나 · CG · G8 틀)
#          잘 덮인 유닛만 : n_calls >= 그 array 중앙값의 절반
#   사람 참고 : v11 L1 (read x 유닛 · b_call = ML>=128 비율 · 관문 통과)
#
#   판정 (분포에서 직접 · 기준은 미리 고정)
#     밀도 = read x 유닛 β 의 가우스 평활(경계 반사) · 폭 h = max(0.035, 0.5/CpG중앙)
#            (CpG 가 8 개면 β 가 0.125 간격 계단이라 폭을 넓힌다)
#     이봉 = 가장 큰 봉우리와 0.30 이상 떨어진 두 번째 봉우리가 있고
#            골 높이 <= 두 번째 봉우리의 절반, 두 번째 봉우리 쪽 분자 >= 2%
#     단봉 = 그 외 (두 번째 봉우리가 있어도 분자 2% 미만이면 '단봉 + 소수 꼬리')
#   여러 개체 제브라피시는 array 마다 6 개체를 합쳤다 (개체 몫 0.8% · §5-D).
#  matplotlib 글리프 함정: ASCII 만.
import io, os, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["font.family"] = "Malgun Gothic"
rcParams["axes.unicode_minus"] = False
rcParams["figure.dpi"] = 110

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


R = _need("RDNA_DATA_ROOT") + "/03_reads"
OUT = _need("RDNA_PREV_OUT_DIR")
HUM = _need("RDNA_HUMAN_L1_PARQUET")
LOG = io.open(OUT + "/X10_READUNIT_PERARRAY.txt", "w", encoding="utf8", newline="\n")
NM = {"zebrafish_GRCz12tu": "제브라피시", "zfinch_v03_dip": "금화조",
      "arabidopsis_ColCC": "애기장대", "bonobo_mPanPan1": "보노보",
      "chimp_mPanTro3": "침팬지", "gorilla_mGorGor1": "고릴라",
      "orang_bornean_mPonPyg2": "보르네오오랑",
      "orang_sumatran_mPonAbe1": "수마트라오랑", "siamang_mSymSyn1": "시아망"}
SEP, TROUGH, MINOR = 0.30, 0.50, 0.02


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


def density(x, h):
    G = 1001
    grid = np.linspace(0, 1, G)
    cnt = np.bincount(np.clip(np.round(x * (G - 1)).astype(int), 0, G - 1), minlength=G).astype(float)
    pad = np.concatenate([cnt[::-1][:-1], cnt, cnt[::-1][1:]])       # 경계 반사
    sd = h * (G - 1)
    k = np.arange(-int(4 * sd) - 1, int(4 * sd) + 2)
    ker = np.exp(-0.5 * (k / sd) ** 2)
    ker /= ker.sum()
    d = np.convolve(pad, ker, mode="same")[G - 1:2 * G - 1]
    return grid, d / (d.sum() * (grid[1] - grid[0]))


def call(x, ncall_med):
    h = max(0.035, 0.5 / max(ncall_med, 1))
    grid, d = density(x, h)
    pk = [i for i in range(len(d))
          if (i == 0 or d[i] >= d[i - 1]) and (i == len(d) - 1 or d[i] > d[i + 1])]
    i1 = max(pk, key=lambda i: d[i])
    far = [i for i in pk if abs(grid[i] - grid[i1]) >= SEP]
    out = dict(h=h, peak1=grid[i1], peak2=np.nan, trough_ratio=np.nan, minor_mass=0.0)
    if not far:
        out["call"] = "단봉"
        return out, grid, d
    i2 = max(far, key=lambda i: d[i])
    a, b = sorted((i1, i2))
    it = a + int(np.argmin(d[a:b + 1]))
    tr = d[it] / d[i2]
    xt = grid[it]
    minor = float((x > xt).mean()) if grid[i2] > xt else float((x < xt).mean())
    out.update(peak2=grid[i2], trough_ratio=tr, minor_mass=minor)
    if tr <= TROUGH and minor >= MINOR:
        out["call"] = "이봉"
    elif tr <= TROUGH:
        out["call"] = "단봉 + 소수 꼬리"
    else:
        out["call"] = "단봉"
    return out, grid, d


P("=" * 104)
P("read 단위 — 각 read 의 각 유닛 메틸 분포 : 종 · array 마다")
P("=" * 104)
P()
P("  판정 기준(미리 고정) : 두 번째 봉우리가 0.30 이상 떨어져 있고, 골 <= 그 봉우리의 절반,")
P("                        그쪽 분자 >= 2%% 이면 '이봉'")
P()

fs = sorted(glob.glob(R + "/*/*.read_unit.tsv*"))
D = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
D = D[D.context == "CG"].copy()
D = D.drop_duplicates(["species_key", "sample", "read_id", "array_id", "array_pos"])
med = D.groupby("array_id").n_calls.transform("median")
D = D[D.n_calls >= 0.5 * med].copy()
D["b"] = D.n_5mC / D.n_calls
D["종"] = D.species_key.map(NM)
ST = pd.read_csv(OUT + "/tables/OTHERSP_ARRAY_STRUCT.tsv", sep="\t")
UL = dict(zip(ST.array_id, ST.unit_len))
NU = dict(zip(ST.array_id, ST.n_unit))

ROWS, PANEL = [], {}

# ---------------------------------------------------------------- 사람 참고
H = pd.read_parquet(HUM, columns=["gate", "b_call", "n_canon", "sample", "slot"])
H = H[H.gate == 1]
xh = H.b_call.values
r, g, d = call(xh, float(np.median(H.n_canon)))
r.update(group="사람", 종="사람 v11 (참고)", array_id="사람 510 하플로타입 합산", unit_len=2240,
         n_unit=np.nan, n=len(xh), reads=np.nan, calls=float(np.median(H.n_canon)),
         mean=xh.mean(), off=(xh < 0.35).mean(), mid=((xh >= 0.35) & (xh <= 0.65)).mean())
ROWS.append(r); PANEL[("사람", r["array_id"])] = (xh, r["calls"], r)
rs = np.random.RandomState(5)
H["sl"] = H["sample"].astype(str) + "_" + H.slot.astype(str)
hc = []
for sl in rs.choice(H.sl.unique(), 20, replace=False):
    gg = H[H.sl == sl]
    rr, _, _ = call(gg.b_call.values, float(np.median(gg.n_canon)))
    hc.append(rr["call"])
P("  [보정] 사람 하플로타입 20 개를 하나씩 : " +
  " · ".join("%s %d" % (k, v) for k, v in pd.Series(hc).value_counts().items()))
P()

# ---------------------------------------------------------------- 타종
for (sp, aid), gdf in D.groupby(["종", "array_id"]):
    x = gdf.b.values
    cm = float(gdf.n_calls.median())
    r, g, d = call(x, cm)
    grp = "유인원" if sp in ("보노보", "침팬지", "고릴라", "보르네오오랑", "수마트라오랑", "시아망") else sp
    r.update(group=grp, 종=sp, array_id=aid, unit_len=UL.get(aid, np.nan), n_unit=NU.get(aid, np.nan),
             n=len(x), reads=gdf.read_id.nunique(), calls=cm, mean=x.mean(),
             off=(x < 0.35).mean(), mid=((x >= 0.35) & (x <= 0.65)).mean())
    if len(x) < 100:
        r["call"] = "표본 부족 (" + r["call"] + ")"
    ROWS.append(r)
    PANEL[(grp, aid)] = (x, cm, r)

T = pd.DataFrame(ROWS)
T.to_csv(OUT + "/tables/OTHERSP_READUNIT_PERARRAY.tsv", sep="\t", index=False, float_format="%.5f")

P("   종          array                         유닛길이 유닛수  read-유닛   CpG  평균   OFF(<.35) 중간   "
  "봉우리1 봉우리2  골비   소수쪽   판정")
for r in T.sort_values(["group", "종", "mean"]).itertuples(index=False):
    P("   %-10s %-30s %6s %5s %9s  %4.0f  %.3f   %.3f   %.3f   %.2f    %s   %s  %s   %s"
      % (r.종, r.array_id.split("__")[-1][:30],
         ("%.0f" % r.unit_len) if r.unit_len == r.unit_len else "-",
         ("%.0f" % r.n_unit) if r.n_unit == r.n_unit else "-",
         format(r.n, ","), r.calls, getattr(r, "mean"), r.off, r.mid, r.peak1,
         ("%.2f" % r.peak2) if r.peak2 == r.peak2 else "  - ",
         ("%.2f" % r.trough_ratio) if r.trough_ratio == r.trough_ratio else "  - ",
         ("%.3f" % r.minor_mass), r.call))
P()
P("  판정 집계")
for grp, gg in T.groupby("group"):
    P("   %-10s array %2d : %s" % (grp, len(gg), " · ".join("%s %d" % (k, v) for k, v in gg.call.value_counts().items())))
P()

# ---------------------------------------------------------------- 그림 : 종마다 한 장 · array 마다 한 칸
COL = {"사람": "#c0392b", "유인원": "#1f5fa9", "금화조": "#e67e22", "애기장대": "#7e57c2", "제브라피시": "#2e7d32"}
PAGES = [("유인원", ["사람", "유인원"]), ("금화조", ["금화조"]), ("애기장대", ["애기장대"]),
         ("제브라피시", ["제브라피시"])]
for page, grps in PAGES:
    keys = [k for k in PANEL if k[0] in grps]
    keys.sort(key=lambda k: (k[0] != "사람", PANEL[k][2]["종"], -PANEL[k][2]["n_unit"] if PANEL[k][2]["n_unit"] == PANEL[k][2]["n_unit"] else 0))
    if page == "제브라피시":
        keys.sort(key=lambda k: (PANEL[k][2]["unit_len"] if PANEL[k][2]["unit_len"] == PANEL[k][2]["unit_len"] else 0,
                                 -(PANEL[k][2]["n"])))
    ncol = 4 if len(keys) <= 12 else 6
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.9 * nrow + 1.0), squeeze=False)
    for i, k in enumerate(keys):
        ax = axs[i // ncol][i % ncol]
        x, cm, r = PANEL[k]
        if cm < 30:
            nb = int(round(cm))
            edges = np.clip((np.arange(nb + 2) - 0.5) / nb, 0, 1)
        else:
            edges = np.linspace(0, 1, 41)
        w = np.ones(len(x)) / len(x)
        ax.hist(x, bins=edges, weights=w, color=COL[k[0]], alpha=0.85)
        ax.set_xlim(0, 1)
        acc = r["array_id"].split("__")[-1]
        acc = acc.split("_")[0] + ":" + (acc.split("_")[1] if "_" in acc else "")
        ttl = "%s  %s" % (r["종"], acc) if k[0] != "사람" else r["array_id"]
        sub = ("유닛 %s bp x %s · read-유닛 %s · CpG %.0f" %
               (("%.0f" % r["unit_len"]) if r["unit_len"] == r["unit_len"] else "-",
                ("%.0f" % r["n_unit"]) if r["n_unit"] == r["n_unit"] else "-",
                format(r["n"], ","), r["calls"]))
        ax.set_title("%s\n%s\n%s  (OFF %.1f%% · 평균 %.2f)"
                     % (ttl, sub, r["call"], 100 * r["off"], r["mean"]), fontsize=8.2)
        ax.tick_params(labelsize=7)
        if i % ncol == 0:
            ax.set_ylabel("read-유닛 비율", fontsize=8)
        if i // ncol == nrow - 1:
            ax.set_xlabel("read 안 유닛 β", fontsize=8)
    for j in range(len(keys), nrow * ncol):
        axs[j // ncol][j % ncol].axis("off")
    fig.suptitle("%s — 각 read 의 각 유닛 메틸(β) 분포 · array 마다 한 칸" % page, fontsize=13)
    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    fig.savefig(OUT + "/figures/Fig_readunit_perarray_%s.png" % page)
    plt.close(fig)
    P("  저장 : 10_other_species/figures/Fig_readunit_perarray_%s.png (칸 %d)" % (page, len(keys)))
LOG.close()
print("done")
