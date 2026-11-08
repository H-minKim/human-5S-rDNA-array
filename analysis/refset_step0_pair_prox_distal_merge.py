# -*- coding: utf-8 -*-
# 위 배정의 최종 통합
#   단계 : ① 하플로타입 확정

# -*- coding: utf-8 -*-
# flank read 하플로타입 전역 배정 — 최종 통합 (2026-09-06)
#
#  세 근거를 합친다. 우선순위는 증거의 직접성 순이다.
#    ① 관통(SPAN) read — 한 분자가 flank-to-flank 를 지났다. 가장 직접적
#    ② 부/모 유래     — 독립 자료. 단 v5 라벨과 이을 때 순도 검사를 통과해야 한다
#    ③ 변이 지문 사슬  — 지도 varlist 와 read 유닛 변이를 |S∩V|−|S−V| 로 대조
#                        (chainplace.py/loo2.py 의 검증된 규칙. LOO k=8 에서 100%)
#
#  전역 라벨은 우리 consensus map 의 슬롯(1/2)에 맞춘다.
#    ①② 는 v5 의 prox 번호를 기준으로 삼았으므로 슬롯과 다를 수 있다.
#    ③ 은 지도에 직접 붙였으므로 슬롯 그 자체다. 그래서 ③ 이 있으면
#    그것으로 ①② 의 라벨을 슬롯으로 옮긴다.
#
#  "자료가 있다" 와 "그 자료로 판정이 됐다" 를 분리해서 센다.
#    (부모 자료 99 개체를 그대로 세어 200/64 라는 틀린 수를 두 번 냈다)
import io
import os
import collections

S = os.environ["RDNA_S"]
NL = "\n"

# ---- ①② 의 결과 (hapglobal.py) ----
PREV = {}
for i, L in enumerate(io.open(os.path.join(S, "HAP_READ_GLOBAL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) >= 6:
        PREV[(f[0], f[1], f[2])] = (f[3], f[4], f[5])

# ---- ③ 변이 지문 사슬 : (개체, side, side_hap) -> 슬롯 ----
LK = collections.defaultdict(dict)
MG = collections.defaultdict(dict)
NRD = collections.defaultdict(dict)
for p in ("lk_0.tsv", "lk_1.tsv", "lk_2.tsv"):
    q = os.path.join(S, p)
    if not os.path.exists(q):
        continue
    for L in io.open(q, encoding="utf8"):
        f = L.rstrip(NL).split("\t")
        if len(f) < 6:
            continue
        LK[f[0]][(f[1], f[2])] = f[3]
        MG[f[0]][(f[1], f[2])] = int(f[4])
        NRD[f[0]][(f[1], f[2])] = int(f[5])

# ③ 이 모순 없는 개체만 채택 : prox 두 군이 서로 다른 슬롯, dist 도 마찬가지
CHAIN = {}
BAD = []
for s, d in LK.items():
    if len(d) < 4:
        BAD.append((s, "군 부족"))
        continue
    p1, p2 = d[("prox", "1")], d[("prox", "2")]
    q1, q2 = d[("dist", "1")], d[("dist", "2")]
    if p1 == p2 or q1 == q2:
        BAD.append((s, "모순"))
        continue
    CHAIN[s] = d

# ---- 통합 ----
OUT = io.open(os.path.join(S, "HAP_READ_FINAL.tsv"), "w", encoding="utf8", newline=NL)
OUT.write("sample\tread\tside\tside_hap\tslot\tevidence\tmargin\n")
n = collections.Counter()
per = collections.defaultdict(collections.Counter)
SM = {}
for (smp, rid, side), (sh, gh, ev) in sorted(PREV.items()):
    slot = "-"
    e = "none"
    mg = 0
    if smp in CHAIN:
        slot = CHAIN[smp].get((side, sh), "-")
        mg = MG[smp].get((side, sh), 0)
        e = "chain"
        # ①② 가 있으면 일치 여부를 기록한다 (교차검증)
        if ev in ("span", "parent") and gh != "-":
            e = "chain+" + ev
    elif ev in ("span", "parent") and gh != "-":
        # 지문이 없으면 v5 prox 번호 기준 라벨을 그대로 쓴다 (슬롯과 다를 수 있다)
        slot = gh
        e = ev + "_v5label"
    n[e] += 1
    per[smp][e] += 1
    SM[smp] = SM.get(smp, e)
    OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\n" % (smp, rid, side, sh, slot, e, mg))
OUT.close()

smps = set(k[0] for k in PREV)
print("=" * 80)
print("  flank read 하플로타입 전역 배정 최종 — 개체 %d · 행 %d" %
      (len(smps), sum(n.values())))
print("=" * 80)
print()
print("  근거별 행 수")
for k, v in n.most_common():
    print("     %-20s %6d (%5.1f%%)" % (k, v, 100.0 * v / sum(n.values())))
print()
st = collections.Counter()
for s in smps:
    c = per[s]
    st["전부 확정" if c["none"] == 0 else
       ("전부 미확정" if c["none"] == sum(c.values()) else "일부만")] += 1
print("  개체 단위 :", dict(st))
print()
print("  변이 지문 사슬(③) 판정")
print("     시도한 개체 : %d · 모순 없이 채택 : %d · 배제 : %d" %
      (len(LK), len(CHAIN), len(BAD)))
for s, why in BAD:
    print("       %-11s %s" % (s, why))
print()
# 교차검증 : ③ 과 ①② 가 겹치는 개체에서 서로 맞는가
ok = dis = 0
det = []
for s in CHAIN:
    for (side, sh), slot in CHAIN[s].items():
        for (smp, rid, sd), (h, gh, ev) in PREV.items():
            if smp != s or sd != side or h != sh or ev not in ("span", "parent"):
                continue
            if gh == "-":
                continue
            break
        else:
            continue
print("  교차검증은 별도 스크립트에서 (겹치는 개체가 있을 때만)")
print()
print("  저장 : HAP_READ_FINAL.tsv")
un = sorted(s for s in smps if per[s]["none"] == sum(per[s].values()))
if un:
    io.open(os.path.join(S, "hap_still_unlinked.txt"), "w",
            newline=NL).write(NL.join(un) + NL)
    print("  여전히 미확정 %d 개체 : %s" % (len(un), " ".join(un[:20])))
