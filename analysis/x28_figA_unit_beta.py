# -*- coding: utf-8 -*-
# 대주제 A 논문용 그림 ① — array 마다 유닛 β 분포 곡선 (막대그래프가 아니라 분포)
# The layout follows that of human Fig. 3 A4 :
#     · β 0~1 을 50 칸(폭 0.02) 히스토그램 · y = 그 칸에 든 유닛 비율(밀도 x 칸폭) · 선형축
#     · 봉우리·꺼진 몫은 표에서 읽어 그 자리에서 쓴다 — 하드코딩 없음
#   층 : read x 유닛(분자) · β = n_5mC/n_calls (하드콜 ML>=128) · context=CG
#   여과 : 잘 덮인 유닛만(n_calls >= array 중앙값 x 0.5) · 보노보 guppy 6.3.8 만 · 시아망 제외
#   부류 세 가지를 나눠 그린다 :
#     (1) 유인원 5 종 — array 안에서 유닛이 이봉
#     (2) 애기장대 — 전부 고메틸(단봉)
#     (3) 금화조·제브라피시 — array 안은 단봉인데 array 끼리 갈린다(array 수준 이봉)
#   산출 : figures/FigA_unit_beta_distribution.{png,svg,pdf} · tables/FIGDATA_A_unit_beta_hist.tsv · tables/FIGDATA_A_array_summary.tsv
import io, os, sys, glob, collections
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams, gridspec

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
LOG = io.open(NEW + "/logs/X28_FIGA.txt", "w", encoding="utf8", newline="\n")
REASON = "대주제 A 논문용 재작성 — 유닛 β 분포(시아망 제외 · 2026-10-01)"
BW, NB = 0.02, 50
OFFCUT = 0.35 # the OFF definition, the same threshold as elsewhere


def P(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


# ---------------------------------------------------------------- 자료
V2 = pd.read_csv(NEW + "/tables/OTHERSP_READUNIT_PERARRAY_v2_bonobo638.tsv", sep="\t")
V2 = V2[V2.array_id.notna() & (V2.종 != "사람 v11 (참고)")]
BC = pd.read_csv(NEW + "/tables/BONOBO_READ_BASECALLER.tsv", sep="\t")
B637 = set(BC.loc[BC.basecaller.astype(str).str.contains("6.3.7"), "read_id"])
ZF = pd.read_csv(B + "/00_catalog/ZEBRAFISH_ALL_ARRAYS.tsv", sep="\t")
ZNAME = dict(zip(ZF.array_id, ZF.chrom_name))
FN = pd.read_csv(B + "/06_aux/zfinch_5S_array_meth_summary.tsv", sep="\t")
FNAME = {}
for t in FN.itertuples(index=False):
    FNAME[(t.chrom_genbank, int(t.array_start), int(t.array_end))] = t.chrom_bigwig

APES = ["고릴라", "보노보", "침팬지", "보르네오오랑우탄", "수마트라오랑우탄"]
SPKEY = {"고릴라": "gorilla_mGorGor1", "보노보": "bonobo_mPanPan1", "침팬지": "chimp_mPanTro3",
         "보르네오오랑우탄": "orang_bornean_mPonPyg2", "수마트라오랑우탄": "orang_sumatran_mPonAbe1",
         "애기장대": "arabidopsis_ColCC", "금화조": "zfinch_v03_dip", "제브라피시": "zebrafish_GRCz12tu"}


def arr_label(aid):
    """array 이름표 — 염색체 이름이 자료에 있는 종만 이름으로, 없으면 접근번호 그대로"""
    tail0 = aid.split("__")[1]
    acc0, s0, e0 = tail0.rsplit("_", 2)
    if aid in ZNAME and isinstance(ZNAME[aid], str):
        return "%s · %.1f Mb" % (ZNAME[aid], int(s0) / 1e6)
    tail = aid.split("__")[1]
    acc, s, e = tail.rsplit("_", 2)
    k = (acc, int(s), int(e))
    if k in FNAME:
        return FNAME[k]
    return acc                                   # 이름 자료 없음 → 접근번호


# ---------------------------------------------------------------- read x 유닛 적재 → array 마다 히스토그램
HIST, SUMM = {}, []
for sp_ko, sp in SPKEY.items():
    fs = sorted(glob.glob(B + "/03_reads/%s/*.read_unit.tsv" % sp) +
                glob.glob(B + "/03_reads/%s/*.read_unit.tsv.gz" % sp))   # 애기장대는 .gz 로 보관돼 있다
    D = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
    D = D[D.context == "CG"].copy()
    med = D.groupby("array_id").n_calls.transform("median")      # 모든 read 로 중앙값 먼저
    D = D[D.n_calls >= 0.5 * med].copy()
    if sp == "bonobo_mPanPan1":
        D = D[~D.read_id.isin(B637)].copy()
    D["b"] = D.n_5mC / D.n_calls
    for aid, g in D.groupby("array_id"):
        h, _ = np.histogram(g.b.values, bins=NB, range=(0, 1))
        HIST[aid] = h / h.sum()
        v = V2[V2.array_id == aid]
        SUMM.append(dict(종=sp_ko, array_id=aid, 이름표=arr_label(aid), n_read유닛=len(g),
                         평균β=round(float(g.b.mean()), 4), 중앙β=round(float(g.b.median()), 4),
                         꺼진몫=round(float((g.b < OFFCUT).mean()), 4),
                         최빈칸=round(float(np.argmax(h) * BW + BW / 2), 3),
                         판정=(v.call.iloc[0] if len(v) else "-"),
                         봉1=(float(v.peak1.iloc[0]) if len(v) else np.nan),
                         봉2=(float(v.peak2.iloc[0]) if len(v) else np.nan),
                         유닛길이=(float(v.unit_len.iloc[0]) if len(v) else np.nan)))
    P("%-10s array %2d · read-유닛 %s" % (sp_ko, D.array_id.nunique(), format(len(D), ",")))
S = pd.DataFrame(SUMM)
P("합계 array %d (시아망 제외) · read-유닛 %s" % (len(S), format(int(S.n_read유닛.sum()), ",")))

rows = []
ctr = np.arange(NB) * BW + BW / 2
for t in S.itertuples(index=False):
    for c, f in zip(ctr, HIST[t.array_id]):
        rows.append(dict(종=t.종, array_id=t.array_id, 이름표=t.이름표, β칸가운데=round(float(c), 3),
                         유닛비율=round(float(f), 6)))
H = pd.DataFrame(rows)
TMP = os.environ.get("TEMP", ".") + "/x28"
os.makedirs(TMP, exist_ok=True)
H.to_csv(TMP + "/FIGDATA_A_unit_beta_hist.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
S.to_csv(TMP + "/FIGDATA_A_array_summary.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")

# ---------------------------------------------------------------- 그림
fig = plt.figure(figsize=(16.5, 12.4))
gs = gridspec.GridSpec(3, 15, figure=fig, height_ratios=[1.0, 1.12, 0.62], hspace=0.62, wspace=2.1,
                       top=0.885, bottom=0.085, left=0.085, right=0.955)
YMAX = 0.30


def draw_curve(ax, aid, color, lw=1.6, label=None, alpha=1.0):
    ax.plot(ctr, HIST[aid], "-", color=color, lw=lw, label=label, alpha=alpha)


# (1) 유인원 5 종 — array 안 이봉
for i, ko in enumerate(APES):
    ax = fig.add_subplot(gs[0, 3 * i:3 * i + 3])
    r = S[S.종 == ko].iloc[0]
    draw_curve(ax, r.array_id, "#B03A2E", lw=2.0)
    ax.fill_between(ctr, 0, HIST[r.array_id], color="#B03A2E", alpha=0.13, lw=0)
    for pk in (r.봉1, r.봉2):
        if np.isfinite(pk):
            ax.axvline(pk, color="#555", lw=0.9, ls=(0, (3, 3)))
            ax.text(pk, YMAX * 0.97, "%.2f" % pk, fontsize=7.4, ha="center", va="top", color="#333")
    ax.text(0.30, 0.78, "꺼진 유닛\n(β<%.2f)\n%.1f%%" % (OFFCUT, 100 * r.꺼진몫), transform=ax.transAxes,
            fontsize=7.8, ha="left", va="top", color="#B03A2E")
    ax.set_xlim(0, 1); ax.set_ylim(0, YMAX)
    ax.set_title("%s\n%s · 유닛 %d bp · n=%s" % (ko, r.이름표, int(r.유닛길이), format(r.n_read유닛, ",")),
                 fontsize=9.2, loc="left")
    ax.set_xlabel("유닛 β", fontsize=8.4)
    if i == 0:
        ax.set_ylabel("유닛 비율 (β 0.02 칸마다)", fontsize=8.8)
    ax.tick_params(labelsize=7.6)
    ax.grid(alpha=0.22, lw=0.5)
fig.text(0.012, 0.930, "a   유인원 5 종 — array 안에서 유닛이 두 봉우리 (꺼진 분자 + 켜진 분자)",
         fontsize=12.5, fontweight="bold", va="top")

# (2)(3) 애기장대 · 금화조 · 제브라
fig.text(0.012, 0.625, "b   애기장대는 전부 고메틸(단봉) · 금화조·제브라피시는 array 안은 단봉이고 array 끼리 갈린다",
         fontsize=12.5, fontweight="bold", va="top")
panels = [("애기장대", "#1E8449", slice(0, 5)), ("금화조", "#B7950B", slice(5, 10)), ("제브라피시", "#2471A3", slice(10, 15))]
OUTLIER = {}        # 이름표는 자료가 정한다 : 그 종 array 평균 β 중앙값보다 0.08 이상 낮은 array
for ko in ("애기장대", "금화조", "제브라피시"):
    g = S[S.종 == ko]
    OUTLIER[ko] = set(g[g.평균β < g.평균β.median() - 0.08].array_id)
for ko, col, sl in panels:
    ax = fig.add_subplot(gs[1, sl])
    g = S[S.종 == ko].sort_values("평균β", ascending=False)
    ymax = max(max(HIST[a]) for a in g.array_id) * 1.15
    for t in g.itertuples(index=False):
        out = t.array_id in OUTLIER[ko]
        draw_curve(ax, t.array_id, "#B03A2E" if out else col, lw=1.8 if out else 1.1, alpha=1.0 if out else 0.55)
    ax.set_xlim(0, 1); ax.set_ylim(0, ymax)
    lo = g[g.array_id.isin(OUTLIER[ko])]
    hi = g[~g.array_id.isin(OUTLIER[ko])]
    note = ("array %d 개 모두 한 덩어리 · 평균 β %.2f~%.2f" % (len(g), g.평균β.min(), g.평균β.max())) if not len(lo) else \
           ("array %d 개 중 %d 개가 떨어져 나옴 (평균 β %.2f~%.2f) · 나머지 %d 개 %.2f~%.2f"
            % (len(g), len(lo), lo.평균β.min(), lo.평균β.max(), len(hi), hi.평균β.min(), hi.평균β.max()))
    ax.set_title("%s\n%s" % (ko, note), fontsize=9.2, loc="left")
    ax.set_xlabel("유닛 β", fontsize=8.6)
    ax.set_ylabel("유닛 비율 (β 0.02 칸마다)", fontsize=8.8)
    ax.tick_params(labelsize=7.8)
    ax.grid(alpha=0.22, lw=0.5)
    for j, t in enumerate(lo.itertuples(index=False)):
        ax.text(0.03, 0.90 - 0.09 * j, "%s  (평균 β %.2f · n=%s)" % (t.이름표, t.평균β, format(t.n_read유닛, ",")),
                transform=ax.transAxes, fontsize=7.8, color="#B03A2E", va="top")

# (4) array 평균 β — array 수준 갈림
ax = fig.add_subplot(gs[2, :])
order = ["제브라피시", "금화조", "애기장대"] + APES[::-1]
colr = {"애기장대": "#1E8449", "금화조": "#B7950B", "제브라피시": "#2471A3"}
rng = np.random.default_rng(0)
for yi, ko in enumerate(order):
    g = S[S.종 == ko]
    c = colr.get(ko, "#B03A2E")
    ax.scatter(g.평균β, yi + (rng.random(len(g)) - 0.5) * 0.26, s=46, color=c, alpha=0.85,
               edgecolor="white", linewidth=0.6, zorder=3)
    ax.text(1.012, yi, "array %d" % len(g), fontsize=8, va="center", color="#444", transform=ax.get_yaxis_transform())
    for t in g.itertuples(index=False):
        if t.array_id in OUTLIER.get(ko, set()):
            ax.annotate(t.이름표, xy=(t.평균β, yi), xytext=(t.평균β, yi + 0.36), fontsize=7.8, ha="center", color="#B03A2E")
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=9)
ax.set_xlim(0, 1.0); ax.set_ylim(-0.6, len(order) - 0.4)
ax.set_xlabel("그 array 의 유닛 평균 β (read x 유닛 · 호출가중 아님 — 유닛마다 같은 무게)", fontsize=8.8)
ax.grid(axis="x", alpha=0.25, lw=0.5)
ax.set_title("c   array 하나를 점 하나로 : 유인원·애기장대는 한 덩어리, 금화조·제브라피시는 떨어져 나온 array 가 있다",
             fontsize=11.6, loc="left", fontweight="bold")

fig.suptitle("대주제 A · 분자(read x 유닛) 수준 메틸 분포 — 타종 5S (시아망 제외 · 2026-10-01)", fontsize=13.5, y=0.975)
fig.text(0.012, 0.012,
         "층 = read x 유닛(분자) · β = n_5mC / n_calls · 5mC = ML>=128(하드콜) · context=CG · 잘 덮인 유닛만(n_calls >= array 중앙값 x 0.5) · 보노보 = guppy 6.3.8 read 만\n"
         "히스토그램 = β 0~1 을 50 칸(폭 0.02) · y 는 그 칸에 든 유닛 비율 · 봉우리 값은 판정표(x10/x19)에서 읽음 · 찍은 값 : tables/FIGDATA_A_unit_beta_hist.tsv · 코드 : code/x28_figA_unit_beta.py",
         fontsize=7.8, color="#444")

base = "FigA_unit_beta_distribution"
for ext in ("png", "svg", "pdf"):
    fig.savefig(TMP + "/%s.%s" % (base, ext))
plt.close(fig)
for f, rel in ((TMP + "/FIGDATA_A_unit_beta_hist.tsv", "tables/FIGDATA_A_unit_beta_hist.tsv"),
               (TMP + "/FIGDATA_A_array_summary.tsv", "tables/FIGDATA_A_array_summary.tsv")):
    arch.put(f, rel, REASON)
for ext in ("png", "svg", "pdf"):
    arch.put(TMP + "/%s.%s" % (base, ext), "figures/%s.%s" % (base, ext), REASON)
if os.path.isdir(SVGDIR):
    import shutil
    shutil.copy2(NEW + "/figures/%s.svg" % base, SVGDIR + "/%s.svg" % base)
P("")
P("그림 : latest/figures/%s.{png,svg,pdf} · 표 2" % base)
LOG.close()
