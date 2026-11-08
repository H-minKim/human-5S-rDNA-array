# -*- coding: utf-8 -*-
# array 메틸 지형 패널 공용 모듈 — FigE4(보충)와 Fig1 본그림이 같은 코드 로 그린다.
#   (두 군데에 같은 그림을 따로 짜 두면 반드시 갈라진다 — 그래서 한 곳에 둔다)
#
# 두 가지 바꿈
#   ① 두 선을 같은 자로
#      - 옛판 : 주황(유닛 층) = 유닛 β 의 중앙값 · 검정(CpG 창) = 호출 가중 평균 (Σ메틸/Σ호출)
#      - 새판 : 주황도 호출 가중 평균 (Σn_meth / Σn_cpg) — 두 선이 같은 뜻이 된다
#      - 실측 차이(자리 10,499) : |가중평균 − 중앙값| 중앙 0.0385 · 최대 0.625 · 가중평균이 더 큰 자리는 4.5%
#        (유닛 β 가 이진 스위치의 비율 이라 치우친 분포이고, 그런 분포에서는 중앙값과 평균이 갈린다)
#   ② 검정 선을 깔끔하게
#      - 창마다 모으는 최소 호출수를 20 → 120 → 300 으로 올린다(`fine_minc` · 지금은 300). 자료를 버리는 것이 아니라
#        같은 자료를 더 넓은 창으로 모은다 — 점 수가 줄고 선이 덜 튄다. 창 폭은 패널·legend 에 적는다.
#
# 입력(읽기만) : tables/PERARRAY_LANDSCAPE.tsv.gz · tables/PERARRAY_META.tsv · tables/ORIENTATION.tsv ·
#                tables/UNITBETA_PERARRAY_v13.tsv.gz · tables/CHROM_NAMES.tsv
import numpy as np
import pandas as pd

NED, MIDB, FLMAX, OFFCUT = 20, 12, 20000, 0.35
BRK_F, SQ = 2000.0, 4.0
FINE, UCOL = "#4D5656", "#E07B39"
FINE_MINC = 300           # 검정 선 한 점이 모으는 최소 CpG 호출수 (옛판 20 → 120 → 300)


def t(u):
    a = np.asarray(u, float)
    mag = np.where(a < 0, np.abs(a), np.where(np.abs(a) <= BRK_F, np.abs(a), BRK_F + (np.abs(a) - BRK_F) / SQ))
    return np.sign(a) * mag


TOUT = float(t(FLMAX))


def win_edges(inside, imax, unit=None):
    """array 안 창은 유닛 폭의 배수 (2026-10-08) — 창이 유닛보다 좁으면 유닛 안 무늬가 되비쳐
    선이 톱니처럼 보인다(그 무늬는 FigG13 이 보여 준다). 유닛 폭으로 맞추면 선이 깔끔해지고
    같은 자리에서 유닛 층(주황)과 같은 것을 재게 된다."""
    if inside and unit:
        e = list(np.arange(0, max(imax, FLMAX) + unit, float(unit)))
    elif inside:
        e = list(np.arange(0, 2000, 50)) + list(np.arange(2000, 10000, 200)) + list(np.arange(10000, max(imax, 10001), 1000))
    else:
        e = list(np.arange(0, 2000, 100)) + list(np.arange(2000, 10000, 500)) + list(np.arange(10000, 60000 + 1, 5000))
    return np.array([x for x in e if x <= FLMAX + 1e-9], dtype=float)


def adaptive(dd, edges, minc=FINE_MINC):
    """격자 창을 그 자리에서만 넓혀 호출수 minc 를 채운다(자료를 채우는 것이 아니라 창 폭만 바꾼다).
    자리 = 호출수 가중 평균 · β = Σ메틸/Σ호출(하드콜)."""
    got = {round(float(a), 1): (float(b), float(c))
           for a, b, c in zip(dd.d.values, dd.호출수.values, dd.메틸호출수.values)}
    xs, ys, nb = [], [], []
    an = am = 0.0
    bag = []
    for e in edges:
        n, m_ = got.get(round(float(e), 1), (0.0, 0.0))
        an += n
        am += m_
        bag.append((float(e), n))
        if an >= minc:
            xs.append(sum(x * w for x, w in bag) / an)
            ys.append(am / an)
            nb.append(len(bag))
            an = am = 0.0
            bag = []
    return np.array(xs), np.array(ys), nb


def slide_unit(dd, unit, step_frac=0.25, lo=0.0, hi=None):
    """array 안 CpG 창 — 창 폭 = 유닛 1 개 · 걸음 = 유닛의 1/4 (2026-10-10)
    왜 : 창이 유닛보다 좁으면 유닛 안 무늬가 되비쳐 선이 톱니처럼 보인다(그 무늬는 Supp3 가 보여 준다).
         창을 유닛 폭으로 두되 미끄러뜨리면 선이 매끄럽고, 점이 유닛 층보다 촘촘해 두 선이 구분된다.
    자리 = 창 가운데 · β = Σ메틸 / Σ호출(하드콜)."""
    d = dd.sort_values("d")
    x = d.d.values.astype(float)
    n = d.호출수.values.astype(float)
    m = d.메틸호출수.values.astype(float)
    step = max(unit * step_frac, 1.0)
    # array 안에는 flank 상한(FLMAX)을 걸지 않는다 — 걸면 유닛이 긴 종(유인원 20 유닛 = 47 kb)에서
    #   선이 20 kb 에서 끊긴다(2026-10-11 사용자 지적으로 찾은 버그).
    top = (x.max() if len(x) else lo) if hi is None else min(hi, x.max() if len(x) else lo)
    cs = np.arange(lo + unit / 2.0, top + step, step)
    xs, ys = [], []
    for c in cs:
        k = (x >= c - unit / 2.0) & (x < c + unit / 2.0)
        if not k.any():
            continue
        tot = n[k].sum()
        if tot <= 0:
            continue
        xs.append(c)
        ys.append(m[k].sum() / tot)
    return np.array(xs), np.array(ys)


def load(NEW):
    """표를 한 번만 읽어 둔다."""
    LS = pd.read_csv(NEW + "/tables/PERARRAY_LANDSCAPE.tsv.gz", sep="\t")
    M = pd.read_csv(NEW + "/tables/PERARRAY_META.tsv", sep="\t")
    O = pd.read_csv(NEW + "/tables/ORIENTATION.tsv", sep="\t").set_index("array_id")
    UB = pd.read_csv(NEW + "/tables/UNITBETA_PERARRAY_v13.tsv.gz", sep="\t",
                     usecols=["묶음", "array_id", "유닛번호", "n_cpg", "n_meth"])
    GRP = dict(zip(UB.array_id, UB.묶음))
    IDS = {}
    for a, k in dict(zip(UB.array_id, UB.묶음)).items():
        IDS.setdefault(k, []).append(a)
    # 유닛 층 = 호출 가중 평균 : 묶음 x 유닛번호 로 Σ메틸 / Σ호출 (금화조는 mat+pat 을 한 묶음으로 더한다)
    SLOT = (UB.groupby(["묶음", "유닛번호"], as_index=False)[["n_cpg", "n_meth"]].sum())
    SLOT["mean_pooled"] = SLOT.n_meth / SLOT.n_cpg
    # 가로축의 자 — 2026-10-27로 염색체 방향 으로 바꾼다
    #   지시 : "main 의 b 그리고 supple 7-2 의 그림 모두 chr 방향으로 그려주고 거기다가
    #           각각 p, d 를 표시해줘(실측해서). 즉 고릴라도 p-d 방향으로 그리고"
    #
    #   옛 자(2026-10-06 ~ 2026-10-27) : 가닥(전사) — `"높은 좌표" if strand == "+" else "낮은 좌표"`.
    #   그때 지시는 "전사 방향을 모든 종에 통일" 이었고, 그래서 고릴라만 뒤집던 예외를 없앴다
    #   (`x9H_figB_fullunit.py:39-42` 에 그 내력이 적혀 있다).
    #   그 자로는 눈금 글자 `P1`/`D1` 이 거짓이 되었다 — 왼쪽을 가닥으로 정해 놓고 글자는
    #   근위/원위라 적었기 때문이다. 49 칸 중 22 칸 에서 반대 끝을 가리켰다
    #   (Fig 7B 4 칸 : 고릴라·금화조 chr9·애기장대 CP116284.1·제브라피시 chr4·36.22 Mb /
    #    Supp 7-2 18 칸). 수(`tables/FIGDATA_B4_pos_summary.tsv` 의 off_fraction_proximal_end / distal_end)는
    #   처음부터 염색체 기준 이었으므로(`x9H:42` "요약표의 근위/원위 몫은 염색체 기준 그대로 둔다"),
    #   축을 수에 맞추는 쪽 으로 통일한다.
    #
    #   근위/원위는 실측값 : `tables/ORIENTATION.tsv` 의 `근위쪽` — 중심절이 있는 쪽 이 근위다
    #   (2026-10-28 고침 · `tables/CENTROMERE.tsv` 에서 잰다). 53 array 전부에 값이 있다(결측 0).
    #   옛 자 *"가까운 염색체 끝의 반대가 근위"* 는 폐기 — 싼 대리지표였고 6 array 가 틀렸다
    #   (제브라피시 chr4 ×3·chr18 ×1 은 긴 팔 · 금화조 chr9 ×2 는 acrocentric).
    #   ⇒ 이제 왼쪽 = 근위 이므로 `landsc.py` 의 `prox`·`P1`/`D1` 이 이름 그대로 가 된다.
    #
    #   딸려 오는 것 : 전사 방향이 더는 한 방향이 아니다.
    #   고릴라만 strand '-' 이고 근위는 다섯 종 모두 높은 좌표라, 고릴라 칸에서만 전사가 왼→오 다.
    #   Fig 7 범례의 "Right to left is the direction of transcription" 은 더는 참이 아니다.
    # On the horizontal axis the lower coordinate is on the left.
    #   지시 : *"a 의 방향 설정은 chromosome 을 chromosome 번호 1번을 왼쪽 끝으로 가게 해서 한 거야.
    #           b 의 경우도 왼쪽이 더 낮은 number 가 가게 해달라는 거야. proximal 과 distal 은
    #           그렇게 배치한 후에 왼쪽과 오른쪽이 각각 근위인지 원위인지 판단한 다음 적는 거고."*
    #   ⇒ 축은 한 가지 자(좌표 오름차순) 로 고정하고, `P1`/`D1` 은 칸마다 실측해서 붙인다(아래 `panel`).
    #   ⇒ 낮 동안 썼던 "근위를 왼쪽" 은 폐기 한다(그 판은 `archive\2026-10-27_v8_before\` 에 있다).
    LEFT = {a: "낮은 좌표" for a in M.array_id}
    _CN = pd.read_csv(NEW + "/tables/CHROM_NAMES.tsv", sep="\t")
    _CH = dict(zip(_CN.접근번호.astype(str), _CN.염색체번호.astype(str)))
    POS = {}
    for _t in M.itertuples(index=False):
        _c = str(_t.chrom).replace("_mat", "").replace("_pat", "")
        _l = str(_t.이름표).replace("_mat", "").replace("_pat", "")
        _nm = _l if _l.lower().startswith("chr") else _CH.get(_c, _c)
        POS[_t.array_id] = _nm if "Mb" in _nm else "%s : %.1f Mb" % (_nm, float(_t.array시작) / 1e6)
    return dict(LS=LS, M=M, O=O, SLOT=SLOT, GRP=GRP, IDS=IDS, LEFT=LEFT, POS=POS)


def panel(ax, L, mt, fine_minc=FINE_MINC, ends=True, ylab=True, fs_small=6.0,
          lw_fine=0.8, lw_unit=1.6, ms_unit=2.2, marker="o", unit_window=False, col_unit=None,
          scalebar=False, fs_bar=6.0, bar_dy=-0.26):
    """선 굵기·점 크기는 칸 크기에 맞춰 줄일 수 있게 인자로 둔다(Fig1 본그림은 칸이 작다)."""
    """한 칸을 그린다. mt = PERARRAY_META 의 한 행(itertuples). 그린 내용을 dict 로 돌려준다."""
    LS, O, SLOT, GRP, IDS, LEFT = L["LS"], L["O"], L["SLOT"], L["GRP"], L["IDS"], L["LEFT"]
    key = GRP.get(mt.array_id)
    ids = IDS.get(key, [mt.array_id])
    prox = LEFT[mt.array_id]
    dist = "높은 좌표" if prox == "낮은 좌표" else "낮은 좌표"
    UNIT = float(mt.유닛길이)
    TINu = abs(float(t(-(NED + 0.5) * UNIT)))
    MIDW = max(2400.0, 0.30 * TINu)
    X_J1, X_J2 = TOUT, TOUT + TINu + MIDW + TINu
    XMAX = X_J2 + TOUT

    def xleft(u):
        return X_J1 - t(u)

    def xright(u):
        return X_J2 + t(u)

    # 음영은 테두리 없이(lw=0) — `color=` 는 면과 선을 함께 칠한다. 굵기를 안 적으면
    #   matplotlib 이 `stroke-width` 를 SVG 에 안 써서 1.0 pt 테두리가 된다(실측 확인).
    ax.axvspan(0, X_J1, color="#F2F2F2", lw=0, zorder=0)
    ax.axvspan(X_J2, XMAX, color="#F2F2F2", lw=0, zorder=0)
    ax.axvspan(X_J1 + TINu, X_J1 + TINu + MIDW, color="#FFF8E7", lw=0, zorder=0)
    for xf in (xleft, xright):
        x0, x1 = sorted([xf(BRK_F), xf(FLMAX)])
        ax.axvspan(x0, x1, color="#E8EAF6", lw=0, zorder=0)
    d = LS[LS.array_id.isin(ids)]
    if len(ids) > 1:          # 묶음이 두 array 면 호출수를 더해 한 곡선으로
        d = (d.groupby(["쪽", "위치", "d"], as_index=False)[["호출수", "메틸호출수"]].sum())
        d["beta"] = d.메틸호출수 / d.호출수
    npt = 0
    for side, xf in ((prox, xleft), (dist, xright)):
        for loc, sgn in (("array 안", -1), ("flank", +1)):
            dd = d[(d.쪽 == side) & (d.위치 == loc)].sort_values("d")
            if not len(dd):
                continue
            # array 안에도 CpG 창을 그린다 — 유닛 폭 창으로 두면 유닛 층과 겹쳐 보이지 않는다(2026-10-10)
            #   대신 창마다 모으는 호출수를 fine_minc(기본 FINE_MINC = 300 · 부르는 쪽이 덮어쓰지 않는다)로 올려 톱니를 줄인다
            if loc == "array 안":
                dx, dy = slide_unit(dd, UNIT)          # 유닛 폭 창을 1/4 유닛씩 미끄러뜨린다
            else:
                ed = win_edges(False, float(mt.array끝 - mt.array시작))
                dx, dy, _nb = adaptive(dd, ed, minc=fine_minc)
            if not len(dx):
                continue
            x = xf(sgn * dx)
            o2 = np.argsort(x)
            ax.plot(x[o2], dy[o2], "-", color=FINE, lw=lw_fine, zorder=6)
            npt += len(dx)
    offs = {"낮은 좌표": float(mt.부분유닛_낮은쪽), "높은 좌표": float(mt.부분유닛_높은쪽)}
    for side, xf in ((prox, xleft), (dist, xright)):
        if offs[side] > 0:
            x0, x1 = sorted([xf(0.0), xf(-offs[side])])
            ax.axvspan(x0, x1, facecolor="#D5D8DC", edgecolor="#7F8C8D", hatch="///", lw=0, alpha=.85, zorder=1)
    # ── 유닛 층 : 호출 가중 평균 (검정 선과 같은 자)
    s = SLOT[SLOT.묶음 == key].copy()
    N = int(mt.유닛수)
    s["낮은쪽번호"] = s.유닛번호 + 1
    s["높은쪽번호"] = N - s.유닛번호
    pcol = "낮은쪽번호" if prox == "낮은 좌표" else "높은쪽번호"
    dcol = "높은쪽번호" if prox == "낮은 좌표" else "낮은쪽번호"
    ux, uy = [], []
    for k in range(1, NED + 1):
        r_ = s[s[pcol] == k]
        if len(r_):
            ux.append(xleft(-(k - .5) * UNIT))
            uy.append(float(r_.n_meth.sum() / max(r_.n_cpg.sum(), 1)))
    mid = s[(s[pcol] > NED) & (s[dcol] > NED)]
    if len(mid):
        q = np.linspace(0, 1, MIDB + 1)
        pr = (mid[pcol].values - NED) / max(N - 2 * NED, 1)
        xm = np.linspace(X_J1 + TINu + MIDW * .08, X_J1 + TINu + MIDW * .92, MIDB)
        for bi in range(MIDB):
            m_ = (pr >= q[bi]) & ((pr < q[bi + 1]) if bi < MIDB - 1 else (pr <= 1))
            if m_.any():
                ux.append(xm[bi])
                uy.append(float(mid.n_meth.values[m_].sum() / max(mid.n_cpg.values[m_].sum(), 1)))
    for k in range(NED, 0, -1):
        r_ = s[s[dcol] == k]
        if len(r_):
            ux.append(xright(-(k - .5) * UNIT))
            uy.append(float(r_.n_meth.sum() / max(r_.n_cpg.sum(), 1)))
    if ux:
        # 점의 테두리 굵기(mew)를 반드시 적는다 — 안 적으면 기본 1.0 pt 가 숨어 들어가
        #   겉지름이 `ms + 1.0` 이 된다(2026-10-19 실측 : 적힌 ms 0.28 ↔ 실제 겉지름 1.28).
        #   사용자가 보고 고른 크기를 그대로 두려고 규격(0.46)으로 적되 ms 를 그만큼 키운다.
        ax.plot(ux, uy, "-" + (marker or ""), color=(col_unit or UCOL), lw=lw_unit,
                ms=ms_unit + 0.54, markeredgewidth=0.46, zorder=5)
    for xj in (X_J1, X_J2):
        ax.axvline(xj, color="#111", lw=0.46, zorder=4)    # 규격 0.46 (2026-10-18)
    ax.axhline(OFFCUT, color="#888", lw=0.41, ls=(0, (4, 3)), zorder=2)   # 안내선은 데이터 선(0.46)보다 가늘게
    ax.set_xlim(0, XMAX)
    ax.set_ylim(0, 1.06)
    # ⑹ x 축에 첫 유닛 두 개 를 표시한다
    #   P1 = 근위 끝에서 첫 유닛 · D1 = 원위 끝에서 첫 유닛 (유닛 층의 첫 점과 같은 자리다)
    #   저녁 — 글자를 칸마다 실측해서 붙인다 (축을 좌표 오름차순으로
    #   놓은 뒤에 왼쪽·오른쪽이 각각 근위인지 원위인지 판단해서 적는다).
    #   축은 늘 `낮은 좌표`가 왼쪽이므로, 근위가 낮은 좌표면 `P1 … D1`, 높은 좌표면 `D1 … P1` 이다.
    #   박아 놓으면 거짓이 된다 — 실제로 낮 동안 49 칸 중 22 칸에서 반대 끝을 가리키고 있었다.
    #   — 근위/원위가 중심절 좌표 로 다시 정해졌다(`x41_orientation.py` 고침).
    #     그 전에는 `ORIENTATION.tsv` 의 `근위쪽` 이 *"가까운 염색체 끝의 반대"* 였고 6 array 가 틀려
    #     49 칸 중 5 칸(Fig 7B 3 · Supp 7-2 2)에서 `P1`/`D1` 이 반대 끝 을 가리켰다.
    #     ⇒ 표가 고쳐졌으므로 이 줄은 그대로 두고, 겹침만 새로 거른다.
    #   `중심절겹침 == "예"` : array 가 중심절 구간 안 이라 근위/원위가 뜻을 잃는다
    #     ⇒ 글자를 적지 않는다 (애기장대 `CP116282.1_15747509_15754203` 한 칸).
    #     틀린 글자를 적는 것보다 안 적는 것이 낫다 — 다음 사람을 틀린 쪽으로 움직이지 않는다.
    _p_ori = O.근위쪽.get(mt.array_id, "낮은 좌표")
    _ovl = str(O.중심절겹침.get(mt.array_id, "아니오")) == "예" if "중심절겹침" in O.columns else False
    _lab = ["", ""] if _ovl else (["P1", "D1"] if _p_ori == "낮은 좌표" else ["D1", "P1"])
    ax.set_xticks([xleft(-0.5 * UNIT), xright(-0.5 * UNIT)])
    ax.set_xticklabels(_lab, fontsize=fs_small)
    ax.tick_params(axis="x", length=1.5, pad=1, width=0.41)
    c_in_l = mt.array시작 if prox == "낮은 좌표" else mt.array끝
    c_in_r = mt.array끝 if prox == "낮은 좌표" else mt.array시작
    c_out_l = c_in_l - FLMAX if prox == "낮은 좌표" else c_in_l + FLMAX
    c_out_r = c_in_r + FLMAX if prox == "낮은 좌표" else c_in_r - FLMAX
    if ends:
        ax.text(0.0, -0.115, "%.2f Mb" % (c_out_l / 1e6), fontsize=fs_small, ha="left",
                transform=ax.transAxes, color="#555")
        ax.text(1.0, -0.115, "%.2f Mb" % (c_out_r / 1e6), fontsize=fs_small, ha="right",
                transform=ax.transAxes, color="#555")
    if ylab:
        ax.set_ylabel("β", fontsize=8)
    if scalebar:
        # x 축은 구간마다 눌러 그린 축이라 눈금을 못 매긴다 → 실제 bp 자 를 막대로 넣는다.
        #   길이는 종마다 다르다(축이 유닛 길이에 매여 있다) — 끝 20 유닛 폭의 1/4 에 가장 가까운 둥근 값.
        _cand = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]
        _bp = min(_cand, key=lambda v: abs(v - TINu / 4.0))
        # ⑸ 자는 그래프 바로 밑 으로 — 칸 안에 두면 자료를 가린다.
        #   x 는 자료 좌표 로 둔다(축이 구간마다 눌려 있어 자가 가리키는 구간 위에 있어야 한다) ·
        #   y 만 축 아래로 내린다(`get_xaxis_transform`). P1 글자와 겹치지 않게 0.30 유닛폭 오른쪽에서 시작한다.
        _tr = ax.get_xaxis_transform()
        _x0 = X_J1 + TINu * 0.30
        ax.plot([_x0, _x0 + _bp], [bar_dy, bar_dy], "-", color="#222", lw=0.46,
                solid_capstyle="butt", zorder=9, clip_on=False, transform=_tr)
        _lab = ("%g kb" % (_bp / 1000.0)) if _bp >= 1000 else ("%d bp" % _bp)
        # 자 글자가 그림과 겹치지 않게
        #   무엇이 틀렸나 : 글자를 막대 가운데(ha="center") 에 놓아 칸이 좁으면 글자가
        #   array 시작 세로선 왼쪽까지 비어져 나가 "1 kb" 가 "| kb" 로 읽혔다(S12 여러 칸 · Fig 7b 도 같았다).
        #   ⇒ ⒜ 막대 왼끝에 ha="left" 로 붙여 세로선 오른쪽에만 두고
        #      ⒝ 흰 바탕 상자 를 깔아 밑의 선·음영이 글자를 가로지르지 않게 한다.
        ax.text(_x0 + _bp + XMAX * 0.012, bar_dy, _lab, fontsize=fs_bar, ha="left", va="center",
                color="#222", zorder=9, clip_on=False, transform=_tr)
    return dict(묶음=key, 유닛층_점=len(ux), CpG창_점=npt, 왼쪽=prox, 유닛길이=UNIT, 유닛수=N)
