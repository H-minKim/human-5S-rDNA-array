# -*- coding: utf-8 -*-
# 00_catalog/5S_unit_status.tsv : array 77 개의 대표 유닛 상태를 한 표로 (2026-09-13 직접 검산 결과 반영)
#   status : CONFIRMED(직접 주기 = 유닛) · FIXED(재산출) · MIXED(두 구역) · SPARSE(유닛 맞으나 반복 드묾)
#            LOW_QUALITY(유닛 동질성 낮음) · NO_UNIT(일정한 주기 없음) · NON_5S(5S 아닌 반복) · EXCLUDED(G8 부적격)
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import glob

NL = chr(10)
TAB = chr(9)
V = "$RDNA_OTHERSP_DIR/data_v2"
LOC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lit")
US = {}
f = io.open(os.path.join(LOC, "UNIT_STATUS.tsv"), encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
iu = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    US[r[0]] = r
P = {}
f = io.open(os.path.join(V, "00_catalog", "REFUNIT_PARAMS.tsv"), encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ip = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    P[r[1]] = r
FIX = {
    "rhesus_T2T_MMU8v2__CM111659.2_20358081_20566882": (
        "FIXED", "유전자 2 개짜리 7,523 bp 이중 유닛으로 재산출 (옛 4,503)",
        "유전자 간격이 2,9xx / 4,4xx 로 교대 · 끊김 28 -> 3 · 유닛 26 · 일치도 중앙 99.79%"),
    "arabidopsis_ColCC__CP116284.1_15916357_15935906": (
        "FIXED", "749 bp 로 재산출 (옛 500)",
        "5S 3' 조각(덮힘 0.58) 25 개를 경계에서 뺐다 · 주기 746 · 유닛 20 · 일치도 중앙 93.19%"),
    "zebrafish_GRCz12tu__CM111182.1_56528778_56649899": (
        "MIXED", "지도 두 개 : 177 단량체 구역 + 477 이중 구역(.k2)",
        "간격이 164/313 교대인 구역(이중 유닛 160 개 · 일치도 100%)과 164 만 이어지는 구역(188 유닛)으로 갈린다"),
    "rhesus_T2T_MMU8v2__CM111659.2_86347364_86428918": (
        "LOW_QUALITY", "지도 유지 · 유닛별 결론 금지",
        "직접 반복 비율 0.14 · 유닛 동질성 중앙 78.8% · G8 끊김 22/34"),
    "rhesus_T2T_MMU8v2__CM111659.2_86448017_86604726": (
        "LOW_QUALITY", "지도 유지 · 유닛별 결론 금지",
        "직접 반복 비율 0.11 · 유닛 동질성 중앙 77.6% · G8 끊김 45/63"),
    "rhesus_T2T_MMU8v2__CM111677.2_38968264_39012812": (
        "NON_5S", "5S array 아님 (이미 기각) · 275 bp 위성 반복",
        "옵셋 사다리 5 단·비율 0.96 로 진짜 탠덤이나 유닛 동질성 77.6% · 5S 유전자와 78-82 bp 만 84-88% 일치"),
    "rhesus_T2T_MMU8v2__CM111659.2_103988607_104060482": (
        "NO_UNIT", "대표 유닛 쓰지 말 것",
        "옵셋 칸 867 개에 최상위도 조각의 6% · G8 도 '유전자가 끝 아님' 으로 부적격"),
    "zebrafish_GRCz12tu__CM111182.1_50334795_50356311": (
        "SPARSE", "유닛(177)은 맞음 · array 의 9% 만 반복",
        "유전자 15 개가 21.5 kb 에 흩어져 있다"),
    "zebrafish_GRCz12tu__CM111182.1_64539016_64551899": (
        "SPARSE", "유닛(177)은 맞음 · array 의 10% 만 반복",
        "유전자 12 개가 12.9 kb 에 흩어져 있다"),
}
w = io.open(os.path.join(V, "00_catalog", "5S_unit_status.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "species_key", "copies", "array_span", "refunit_len", "genes_per_unit",
                  "period_direct", "repeat_frac", "units_median_identity", "g8_units", "g8_breaks",
                  "status", "action", "evidence"]) + NL)
import collections
cnt = collections.Counter()
rows = []
for aid in sorted(US):
    r = US[aid]
    pp = P.get(aid)
    gpu = "2" if aid in ("rhesus_T2T_MMU8v2__CM111659.2_20358081_20566882",) else "1"
    if aid == "zebrafish_GRCz12tu__CM111182.1_56528778_56649899":
        gpu = "1+2"
    rl = pp[ip["refunit_len"]] if pp else "NA"
    mid = pp[ip["units_median_identity"]] if pp else "NA"
    gu = pp[ip["n_units_in_range"]] if pp else "NA"
    gb = pp[ip["n_break"]] if pp else "NA"
    if aid in FIX:
        st, act, ev = FIX[aid]
    elif r[iu["g8_usable"]] != "Y":
        ev0 = "G8 부적격 : " + (pp[ip["note"]] if pp else "")
        pdv = int(r[iu["period_direct"]]) if r[iu["period_direct"]].isdigit() else 0
        rlv = int(rl) if rl.isdigit() else 0
        if pdv and rlv and abs(pdv - rlv) <= max(0.10 * rlv, 30) and float(r[iu["repeat_frac"]]) >= 0.30:
            ev0 += " | 단 직접 측정은 유닛 %d 를 지지한다 (옵셋 %d 에 정렬 비율 %s) — copy 가 적어 제외된 것" % (rlv, pdv, r[iu["repeat_frac"]])
        st, act, ev = "EXCLUDED", "대표 유닛·지도 없음", ev0
    elif r[iu["status"]].startswith("확인"):
        st, act, ev = "CONFIRMED", "없음", "직접 주기 %s = 유닛 %s · 반복 비율 %s" % (r[iu["period_direct"]], rl, r[iu["repeat_frac"]])
    elif r[iu["status"]].startswith("유닛 동질성"):
        st, act, ev = "LOW_QUALITY", "유닛별 결론 조심", "유닛 동질성 중앙 %s%% (직접 주기 %s = 유닛 %s)" % (mid, r[iu["period_direct"]], rl)
    else:
        pdv = int(r[iu["period_direct"]]) if r[iu["period_direct"]].isdigit() else 0
        rlv = int(rl) if rl.isdigit() else 0
        if pdv and rlv and abs(pdv - rlv) <= max(0.10 * rlv, 30):
            st, act = "CONFIRMED", "없음"
            ev = "직접 주기 %d = 유닛 %d · 반복 비율 %s · (유전자 간격 차례는 불규칙)" % (pdv, rlv, r[iu["repeat_frac"]])
        else:
            st, act, ev = "REVIEW", "검토", r[iu["status"]] + " | " + r[iu["note"]]
    cnt[st] += 1
    rows.append(aid)
    w.write(TAB.join([aid, r[iu["species_key"]], r[iu["copies"]], r[iu["array_span"]], rl, gpu,
                      r[iu["period_direct"]], r[iu["repeat_frac"]], mid, gu, gb, st, act, ev]) + NL)
w.close()
print("5S_unit_status.tsv : %d 행" % len(rows))
for k, v in cnt.most_common():
    print("   %-12s %d" % (k, v))
