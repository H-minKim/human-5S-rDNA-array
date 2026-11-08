# -*- coding: utf-8 -*-
# 최종 판정 : array 마다 (1) 직접 측정 주기 (2) 유닛 하나에 든 유전자 수 k (3) 지금 G8 유닛이 맞는가
#   앞 판들의 결함 정정
#     - 기본 주기를 "가장 작은 사다리 옵셋" 으로 고르면 부차 봉우리(240 등)에 끌려간다
#       -> 합친 정렬량이 가장 큰 옵셋을 잡고, 그것의 약수 가운데 정렬량이 0.9 배 이상인 것이 있으면 그것으로 내린다
#     - 간격을 부류로 나누는 방식은 넓은 분포를 잘게 쪼갠다
#       -> k = round(직접주기 / 유전자 간격 중앙값) 으로 센다 (간격은 중앙값 2 배 이하만)
#   산출 : UNIT_STATUS.tsv (77 행) · FIXLIST.tsv (고칠 array 만)
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import glob
import statistics
import collections

NL = chr(10)
TAB = chr(9)
BASE = "$RDNA_OTHERSP_DIR/data_v2"
LOC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lit")
ST, P, OFF = {}, {}, {}
f = io.open(os.path.join(BASE, "00_catalog", "5S_array_structure.tsv"), encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ix = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    ST["%s__%s_%s_%s" % (r[ix["species"]], r[ix["chrom"]], r[ix["array_start"]], r[ix["array_end"]])] = r
f = io.open(os.path.join(BASE, "00_catalog", "REFUNIT_PARAMS.tsv"), encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ipr = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    P[r[1]] = r
f = io.open(os.path.join(LOC, "OFFHIST.tsv"), encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    d = {}
    if len(r) > 3 and r[3]:
        for kv in r[3].split(","):
            a, b = kv.split(":")
            d[int(a)] = int(b)
    OFF[r[0]] = (int(r[1]), int(r[2]), d)
HIT = {}
for p in glob.glob(os.path.join(BASE, "cms_tbl", "*.tbl")):
    sp = os.path.basename(p).replace(".tbl", "")
    H = collections.defaultdict(list)
    for L in io.open(p, encoding="utf8"):
        if L.startswith("#"):
            continue
        ff = L.split()
        if len(ff) < 9:
            continue
        a, b = int(ff[7]), int(ff[8])
        H[ff[0]].append((min(a, b), max(a, b)))
    for c in H:
        H[c].sort()
    HIT[sp] = H


def fundamental(off, slen):
    if not off:
        return 0, 0
    def ws(p):
        return off.get(p - 20, 0) + off.get(p, 0) + off.get(p + 20, 0)
    keys = [k for k in sorted(off) if k >= 130]
    if not keys:
        return 0, 0
    top = max(keys, key=ws)
    best = top
    for q in keys:
        if q >= top:
            break
        m = top / float(q)
        if abs(m - round(m)) <= 0.06 and round(m) >= 2 and ws(q) >= 0.90 * ws(top):
            best = q
            break
    return best, ws(best)


w = io.open(os.path.join(LOC, "UNIT_STATUS.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "species_key", "copies", "array_span", "g8_refunit_len", "g8_usable",
                  "g8_units_in_range", "g8_breaks", "g8_median_identity", "period_direct",
                  "repeat_frac", "gap_med", "genes_per_unit_k", "unit_len_expected",
                  "status", "action", "note"]) + NL)
fx = io.open(os.path.join(LOC, "FIXLIST.tsv"), "w", encoding="utf8", newline=NL)
rows = []
for aid in sorted(OFF):
    if aid not in ST:
        continue
    r = ST[aid]
    sp = aid.split("__")[0]
    slen, nh, off = OFF[aid]
    pdv, apd = fundamental(off, slen)
    frac = apd / float(slen) if slen else 0.0
    c, s0, e1 = r[ix["chrom"]], int(r[ix["array_start"]]), int(r[ix["array_end"]])
    hs = [h for h in HIT.get(sp, {}).get(c, []) if s0 <= h[0] and h[1] <= e1]
    ends = sorted(h[1] for h in hs)
    gaps = [ends[i + 1] - ends[i] for i in range(len(ends) - 1)]
    gm = 0
    if gaps:
        m0 = statistics.median(gaps)
        keep = [g for g in gaps if g <= 2 * m0]
        gm = int(statistics.median(keep)) if keep else int(m0)
    pp = P.get(aid)
    g8 = pp[ipr["refunit_len"]] if pp and pp[ipr["refunit_len"]] != "NA" else "NA"
    us = pp[ipr["refunit_usable"]] if pp else "N"
    uir = pp[ipr["n_units_in_range"]] if pp else "0"
    brk = pp[ipr["n_break"]] if pp else "NA"
    mid = pp[ipr["units_median_identity"]] if pp else "NA"
    # k = 유닛 하나에 든 유전자 수 : 간격 차례의 자기상관으로 센다
    #   단량체면 gap[i] ~ gap[i+1] · 이량체면 gap[i] != gap[i+1] 이지만 gap[i] ~ gap[i+2]
    k, ule, agr = 0, 0, []
    for j in (1, 2, 3, 4):
        if len(gaps) <= j + 2:
            break
        ok = sum(1 for i in range(len(gaps) - j)
                 if abs(gaps[i] - gaps[i + j]) <= 0.08 * max(gaps[i], gaps[i + j]))
        fr = ok / float(len(gaps) - j)
        agr.append("%d:%.2f" % (j, fr))
        if fr >= 0.80 and k == 0:
            k = j
    if k:
        sums = [sum(gaps[i:i + k]) for i in range(0, len(gaps) - k + 1)]
        sums.sort()
        ule = int(sums[len(sums) // 2])
    cp = int(r[ix["copies"]])
    if us != "Y":
        status, action = "제외 (G8 대표유닛 부적격)", "없음"
        note = pp[ipr["note"]] if pp else "params 없음"
    elif not pdv or frac < 0.15:
        status, action = "유닛 없음 (일정한 주기 없음)", "대표유닛·지도 제외"
        note = "직접 주기 정렬량이 조각의 %.2f 뿐" % frac
    elif k == 0:
        status, action = "간격 불규칙 (유닛 경계 흔들림)", "검토 필요"
        note = "간격 차례 자기상관 %s · 간격 중앙 %d · 직접 주기 %d" % ("/".join(agr), gm, pdv)
    elif k >= 2:
        status, action = "%d 량체 (유전자 %d 개/유닛)" % (k, k), "유전자 %d 개 틀로 재산출" % k
        note = "간격 차례 자기상관 %s · 유닛 길이 %d · 직접 주기 %d" % ("/".join(agr), ule, pdv)
    elif g8 != "NA" and (abs(int(g8) - ule) <= max(0.10 * int(g8), 30) or abs(int(g8) - pdv) <= max(0.10 * int(g8), 30)):
        status, action = "확인 (유닛 = 직접 주기)", "없음"
        note = ""
    else:
        status, action = "G8 유닛과 직접 주기 어긋남", "직접 주기로 재산출"
        note = "G8 %s · 유닛(간격) %d · 직접 %d" % (g8, ule, pdv)
    if mid != "NA" and mid != "" and us == "Y":
        try:
            if float(mid) < 90 and "확인" in status:
                status = "유닛 동질성 낮음 (%s%%)" % mid
                action = "검토 필요"
        except ValueError:
            pass
    rows.append((aid, sp, cp, e1 - s0 + 1, g8, us, uir, brk, mid, pdv, frac, gm, k, ule, status, action, note))
    w.write(TAB.join([aid, sp, str(cp), str(e1 - s0 + 1), g8, us, uir, brk, mid, str(pdv),
                      "%.3f" % frac, str(gm), str(k), str(ule), status, action, note]) + NL)
    if action not in ("없음",):
        fx.write(TAB.join([aid, status, str(k if k >= 2 else 1), action]) + NL)
w.close()
fx.close()
print("%-50s %5s %6s %6s %5s %4s %6s  %s" % ("array", "copy", "G8유닛", "직접", "비율", "k", "간격중앙", "판정"))
for t in rows:
    if t[14].startswith("확인"):
        continue
    print("%-50s %5d %6s %6d %5.2f %4d %6d  %s | %s" % (t[0].split("__")[1], t[2], t[4], t[9], t[10], t[12], t[11], t[14], t[15]))
print()
print("집계:", dict(collections.Counter(t[14].split(" (")[0] for t in rows)))
