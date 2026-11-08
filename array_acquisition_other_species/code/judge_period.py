# -*- coding: utf-8 -*-
# 옵셋 히스토그램에서 기본 주기와 고차반복(HOR) 주기를 따로 판정한다
#   기본 주기 = 사다리(p·2p·3p…)를 이루는 가장 작은 옵셋
#   HOR = 정렬량이 가장 큰 옵셋 (기본 주기의 배수면 그 배수를 적는다)
#   그리고 내 대표 유닛 길이와 견준다
#   산출 : PERIOD_FINAL.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import glob
import collections

NL = chr(10)
TAB = chr(9)
B = os.path.dirname(os.path.abspath(__file__))
L = os.path.join(B, "lit")
CAN = "$RDNA_OTHERSP_DIR/data_v2/01_canonical"
UL = {}
for p in glob.glob(os.path.join(CAN, "*", "*.unit.fa")):
    aid = os.path.basename(p).replace(".unit.fa", "")
    UL[aid] = len("".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]))
ST = {}
f = io.open("$RDNA_OTHERSP_DIR/data_v2/00_catalog/5S_array_structure.tsv", encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ix = dict((n, i) for i, n in enumerate(hd))
for LL in f:
    r = LL.rstrip(NL).split(TAB)
    ST["%s__%s_%s_%s" % (r[ix["species"]], r[ix["chrom"]], r[ix["array_start"]], r[ix["array_end"]])] = r
TC = {}
if os.path.exists(os.path.join(L, "TILE_CHECK.tsv")):
    g = io.open(os.path.join(L, "TILE_CHECK.tsv"), encoding="utf8")
    h2 = g.readline().rstrip(NL).split(TAB)
    i2 = dict((n, i) for i, n in enumerate(h2))
    for LL in g:
        r = LL.rstrip(NL).split(TAB)
        TC[r[0]] = r
    TCIX = i2
else:
    TCIX = {}

w = io.open(os.path.join(L, "PERIOD_FINAL.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "copies_cms", "array_span", "slice_len", "n_hsp", "my_unit_len",
                  "period_direct", "aln_at_period", "frac_of_slice", "n_rungs",
                  "hor_period", "hor_multiple", "unit_copies_tiled", "unit_cov_tiled", "unit_pident",
                  "verdict_repeat", "verdict_unit"]) + NL)
rows = []
f = io.open(os.path.join(L, "OFFHIST.tsv"), encoding="utf8")
f.readline()
for LL in f:
    r = LL.rstrip(NL).split(TAB)
    aid, slen, nh = r[0], int(r[1]), int(r[2])
    hist = r[3] if len(r) > 3 else ""
    off = {}
    if hist:
        for kv in hist.split(","):
            a, b = kv.split(":")
            off[int(a)] = int(b)
    myu = UL.get(aid, 0)
    cps = ST[aid][ix["copies"]] if aid in ST else "NA"
    span = ST[aid][ix["array_end"]] if aid in ST else "NA"
    span = (int(ST[aid][ix["array_end"]]) - int(ST[aid][ix["array_start"]]) + 1) if aid in ST else 0
    fund, af, nr, hor, horm = 0, 0, 0, 0, 0
    if off:
        # 20 bp 바구니 경계 때문에 대각선이 두 칸에 걸린다 -> 인접 칸을 합쳐서 본다
        def wsum(p):
            return off.get(p - 20, 0) + off.get(p, 0) + off.get(p + 20, 0)
        keys = sorted(off)
        mw = max(wsum(k) for k in keys)
        cand = [k for k in keys if k >= 130 and wsum(k) >= 0.35 * mw]
        # 같은 봉우리의 옆칸 중복 제거 (더 큰 쪽만)
        cand = [k for k in cand if wsum(k) >= wsum(k - 20) and wsum(k) >= wsum(k + 20)]
        for q in cand:
            tol = max(40, int(0.04 * q))
            rungs = 0
            for k in range(1, 13):
                d = q * k
                if d > slen * 0.8 or d > 30000:
                    break
                hit = [o for o in keys if abs(o - d) <= tol and wsum(o) >= 0.20 * mw]
                if not hit:
                    break
                rungs += 1
            if rungs >= 3:
                fund, af, nr = q, wsum(q), rungs
                break
        if not fund and cand:
            fund, af, nr = cand[0], wsum(cand[0]), 1
        if fund:
            # 진짜 고차반복 = 기본 주기의 배수인데 합친 정렬량이 기본보다 확실히 큰 것
            for k in keys:
                if k < fund * 1.8:
                    continue
                q = k / float(fund)
                if abs(q - round(q)) <= 0.05 and wsum(k) > 1.25 * af:
                    hor, horm = k, int(round(q))
    frac = af / float(slen) if slen else 0.0
    t = TC.get(aid)
    ucp = t[TCIX["U_copies"]] if t else "NA"
    ucov = t[TCIX["U_cov"]] if t else "NA"
    upid = t[TCIX["U_pident_med"]] if t else "NA"
    if not fund:
        vr = "반복 신호 없음"
    elif nr >= 3 and frac >= 0.50:
        vr = "탠덤 확인"
    elif nr >= 3 and frac >= 0.15:
        vr = "일부만 탠덤"
    elif nr >= 3:
        vr = "약한 탠덤"
    else:
        vr = "사다리 없음"
    if not fund or not myu:
        vu = "판정 불가"
    elif abs(fund - myu) <= max(0.10 * myu, 25):
        vu = "유닛 = 주기 (맞음)"
    elif myu > fund * 1.6:
        vu = "유닛이 주기의 %.1f 배 (덩어리)" % (myu / float(fund))
    elif fund > myu * 1.6:
        vu = "주기가 유닛의 %.1f 배" % (fund / float(myu))
    else:
        vu = "어긋남 (주기 %d · 유닛 %d)" % (fund, myu)
    rows.append((aid, cps, span, slen, nh, myu, fund, af, frac, nr, hor, horm, ucp, ucov, upid, vr, vu))
    w.write(TAB.join([aid, str(cps), str(span), str(slen), str(nh), str(myu), str(fund), str(af),
                      "%.3f" % frac, str(nr), str(hor), str(horm), ucp, ucov, upid, vr, vu]) + NL)
w.close()
print("%-50s %5s %6s %6s %6s %5s %5s %5s  %-12s %s" % ("array_id", "copy", "unit", "주기", "p정렬", "비율", "단", "HOR", "반복", "유닛"))
for t in rows:
    aid, cps, span, slen, nh, myu, fund, af, frac, nr, hor, horm, ucp, ucov, upid, vr, vu = t
    hs = ("%d(x%d)" % (hor, horm)) if horm > 1 else "-"
    print("%-50s %5s %6d %6d %6d %5.2f %5d %5s  %-12s %s" % (aid, cps, myu, fund, af, frac, nr, hs, vr, vu))
print()
print("반복 판정:", dict(collections.Counter(t[-2] for t in rows)))
print("유닛 판정:", dict(collections.Counter(t[-1].split(" (")[0] for t in rows)))
