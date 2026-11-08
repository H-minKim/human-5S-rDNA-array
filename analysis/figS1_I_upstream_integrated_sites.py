# -*- coding: utf-8 -*-
# → INTEGRATED_13SITES.tsv (기탁본은 14 행 뒤판 · 그 소스는 미보존)
#   단계 : S1 I 상류
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 통합 : 삽입/결실(앵커·DP 축) + 구조 이상(배치자리 축 · 역위 양성대조 포함) 을 한 그림으로
# Output : INTEGRATED_STRUCTURAL.md and INTEGRATED_13SITES.tsv
import io, os, collections, statistics
K = os.environ["RDNA_HANDOFF_AD6468DD"]; B = chr(92)
D = os.path.join(K, "result_placement_slot_anchor_20260930")
A = os.path.join(os.environ["RDNA_V11"], "09_tests", "unit_indel_all_v2_20260916",
             "anchor_check_20260918", "results", "ANCHOR_SITES.tsv")
TAB = chr(9); NL = chr(10)
def rows(p):
    f = io.open(p, encoding="utf8"); h = f.readline().rstrip(NL).split(TAB)
    for L in f:
        if L.strip(): yield dict(zip(h, L.rstrip(NL).split(TAB)))
rate = {};
for r in rows(os.path.join(D, "ALIGNED_POS_RATE.tsv")):
    rate[(r["sample"], r["slot"], int(r["array_pos"]))] = float(r["local_best_rate"])
byhap = collections.defaultdict(list)
for (s, sl, ap), v in rate.items(): byhap[(s, sl)].append(v)
hub = {("HG01258","1","D1"):("모자이크","확정",10,6,"0.9865"), ("HG03516","2","D4"):("모자이크","확정",17,5,"0.9950"),
 ("HG04115","2","D2"):("모자이크","보류",13,10,"0.9897"), ("HG00609","2","I1"):("모자이크","보류",9,8,"0.9919"),
 ("HG01150","1","I2"):("모자이크","보류",3,5,"0.9932"), ("HG00609","2","I13"):("반박","크기 섞임",10,10,"0.9948"),
 ("NA20827","2","I8"):("반박","정본형 다수",9,15,"0.9850"), ("HG01106","1","D8"):("read 동의","§3",14,6,"0.9567"),
 ("NA20850","1","D6"):("read 동의","§3",12,5,"0.9364"), ("HG01081","2","D9"):("read 동의","§3",4,4,"0.9578"),
 ("HG02559","2","D6"):("read 동의","§3 · 정본형 0",6,0,"–"), ("HG03270","1","D1"):("read 동의","수리됨",7,2,"0.9850"),
 ("NA21133","2","D22"):("재현 안 됨","–",None,None,"–")}
dec = {("HG03270","1"):"수리함 v11-수정1(09-17)", ("HG01081","2"):"수리 안 함 09-22 '냅두자'", ("HG01106","1"):"수리 안 함 09-22 '냅두자'",
 ("HG02559","2"):"수리 안 함 09-22 '냅두자'", ("NA20850","1"):"수리 안 함 09-22 '냅두자'", ("NA21133","2"):"수리 안 함(재현안됨)"}
out = []
for r in rows(A):
    if not (int(r["n_reads"]) >= 3 and float(r["support_frac"]) >= 0.7) or r["sample"] == "HG00272": continue
    s, sl, t = r["sample"], r["slot"], r["type"] + r["size_units"]; lo, hi = int(r["ap_lo"]), int(r["ap_hi"])
    h = hub[(s, sl, t)]
    aps = [ap for ap in range(lo - 2, hi + 3) if (s, sl, ap) in rate]
    mn = min(rate[(s, sl, ap)] for ap in aps) if aps else None
    med = statistics.median(byhap.get((s, sl), [0]))
    out.append((h[0], s, sl, t, lo, hi, int(r["n_reads"]), int(r["n_spanning"]), h[2], h[3], h[4],
                dec.get((s, sl), "수리 안 함 09-22(모자이크·조립본 형 유지)"), mn, med))
order = {"read 동의": 0, "모자이크": 1, "반박": 2, "재현 안 됨": 3}
out.sort(key=lambda x: (order[x[0]], x[1]))
g = io.open(os.path.join(D, "INTEGRATED_13SITES.tsv"), "w", encoding="utf8", newline=NL)
g.write(TAB.join(["class", "sample", "slot", "event", "ap_lo", "ap_hi", "our_n_reads", "our_n_spanning",
                  "hub_event_reads", "hub_canon_reads", "hub_canon_identity", "decision",
                  "newscan_min_local_best", "hap_median_local_best", "newscan_recovered"]) + NL)
for o in out:
    g.write(TAB.join([o[0], o[1], o[2], o[3], str(o[4]), str(o[5]), str(o[6]), str(o[7]),
                      "-" if o[8] is None else str(o[8]), "-" if o[9] is None else str(o[9]), o[10], o[11],
                      "-" if o[12] is None else "%.2f" % o[12], "%.2f" % o[13],
                      "-" if o[12] is None else ("Y" if o[12] <= 0.5 else "N")]) + NL)
g.close()
T = []
T.append("# 통합 그림 — 구조적 불일치 하나로 (2026-09-30)")
T.append("")
T.append("이 세션이 찾던 것은 read 와 정본 사이의 **구조적 불일치**다. 두 축으로 봤고, 각 축의 **검출력을 따로 실측**했다.")
T.append("")
T.append("| | 축 A — 분자 사건 (앵커 · DP) | 축 B — 정본 구조 (배치 자리) |")
T.append("|---|---|---|")
T.append("| 무엇을 보나 | read 한 분자가 정본보다 유닛이 **몇 개 많고 적은가** | 정본 자리 하나하나가 **read 와 맞는가·차례가 맞는가** |")
T.append("| 어떻게 | 고유 앵커 두 개 사이 **bp 거리**가 유닛 정수배만큼 바뀌나 | 배치된 자리의 자카드가 **±5 이웃 중 최고인가**(국소 최고율) |")
T.append("| 검출력(실측) | 길이 분포의 **81.6%가 유닛 정수배**(균등 10.7%) · 밀림 인공산물 상한 ~11% | **심어 넣은 역위 검출** : 3 유닛 **91.8%** · 6 유닛 **97.4%** · 12 유닛 **98.6%**(거짓양성 1.0~1.3% · 구간 밖 0.97→0.97 불변) |")
T.append("| 사각지대 | 동점 53.7% 라 개별 유닛 자리는 특정 못 함(사건 판정은 앵커 사이 거리로만) | **1 유닛 사건**(문턱 0.50 미달) · 깊이 0~4 자리 7.8% · 22 하플로 · **소수 read 분자 수준** |")
T.append("| 결과 | 자리 단위 **ONT 17.81** [13.96–22.39] · **HiFi 4.21** [1.93–8.00] /10 만 read×유닛 · 정본 오류 의심 **13 자리** | 저율 구간 56 → 5 축 판정 → read 합의 대조 → **구조(역위·블록 이동·중복) 0 건** |")
T.append("")
T.append("## 두 축을 맞물린 교차검증")
T.append("- 축 B 로 축 A 의 결과를 다시 보면 **정본 오류 의심 13 중 11(85%)** 이 국소최고율 ≤0.5 로 재검출된다(최소값 중앙 **0.19** · 대조 31%·0.80).")
T.append("- 못 잡은 2 곳은 **둘 다 1 유닛 결실**(HG03270 h1 0.55 · HG01258 h1 0.70). ⇒ 축 B 의 민감도 경계 = 1 유닛.")
T.append("- 반대로 축 A 의 **주 값 73 자리는 축 B 로 7%** 만 보인다(대조 8%) — 축 B 는 자리별 read **평균**이라 소수 read 사건은 묻힌다. **두 축은 상보적이고, 어느 하나로 다른 하나를 대신할 수 없다.**")
T.append("- ⇒ 축 B 의 '역위 0 건' 은 **정본 수준** 구조에 대한 음성이다. 소수 read 에만 있는 **분자 수준 역위는 미측정**이다.")
T.append("")
T.append("## 정본 오류 의심 13 자리 — 성격 · 결정 · 축 B 재검출 (전건)")
T.append("")
T.append("| 성격 | 자리 | 사건 | 우리 지지 | 사건형/정본형 read (정본형 일치도) | 수리 결정 | 축 B 국소최고율 최소 (하플로 중앙) |")
T.append("|---|---|---|---|---|---|---|")
for o in out:
    T.append("| %s | %s h%s | %s %d–%d | %d/%d | %s / %s (%s) | %s | %s (%.2f) %s |" % (
        o[0], o[1], o[2], o[3], o[4], o[5], o[6], o[7],
        "-" if o[8] is None else o[8], "-" if o[9] is None else o[9], o[10], o[11],
        "-" if o[12] is None else "**%.2f**" % o[12], o[13], "" if o[12] is None else ("재검출" if o[12] <= 0.5 else "못잡음")))
T.append("")
T.append("## 한 문단으로")
T.append("구조적 불일치는 이 세션이 다룬 범위에서 **유닛 단위 삽입/결실 한 가지 모양**으로만 나타났다. "
         "그 사건들은 실재하고(길이의 81.6%가 유닛 정수배 · HiFi 교차 확인), 그중 **read 전원이 정본에 반대하는 5 자리**가 정본 오류에 해당하며 **2 자리(HG00133 h2 · HG03270 h1)가 실제로 수리**됐고 나머지는 **2026-09-22 에 '수리 안 함' 으로 결정**됐다. "
         "**역위 · 블록 재배치 · 중복은 전용 검정으로 0 건**이다 — 그 검정은 정본에 6 유닛 역위를 심으면 **97.4%** 로 잡고(거짓양성 1.2%) 기존 정본 오류 13 중 11 을 다시 찾으므로, 이 음성은 **근거 있는 음성**이다. "
         "다만 음성의 범위는 **정본 자리의 92.2% · 594/616 하플로 · 정본 수준**이며, **1 유닛 사건과 소수 read 분자 수준 역위는 이 검정의 밖**이다.")
T.append("")
T.append("## 파일")
T.append("`INTEGRATED_13SITES.tsv`(위 표 원자료) · 근거 : `ALIGNED_SCAN_ALL.txt`(역위 검출력) · `RECOVER_INDELS.txt`(교차검증) · `ALIGNED_COVERAGE.txt`(분모) · `RUN_CONSENSUS.tsv` · `STRUCTURAL_CHECK.txt`")
p = os.path.join(D, "INTEGRATED_STRUCTURAL.md")
io.open(p, "w", encoding="utf8", newline=NL).write(NL.join(T) + NL)
print("ctrl", sum(1 for c in io.open(p, encoding="utf8").read() if ord(c) < 32 and c not in "\n\r\t"))
print(NL.join(T[:10]))
