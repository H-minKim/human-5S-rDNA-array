# -*- coding: utf-8 -*-
# 결손 자리 재검출
#   단계 : S1 I 상류
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 교차검증 : **새 방법(배치된 자리 국소 최고율)** 으로 기존에 찾은 삽입/결실을 다시 찾을 수 있나?
#   대상 세 묶음
# (1) the 13 positions where the canonical sequence is suspected to be wrong (ANCHOR_SITES, read>=3, support>=0.7); this should be the strongest signal
#     (2) 주 값에 쓰인 ONT 자리 73 (SITES_ONT.tsv)
# (3) positions flagged by the v2 DP pass (UINDEL_SITES, class = canonical suspect)
#   지표 : 그 자리(및 ±2)의 국소 최고율이 하플로 기준선보다 떨어지나 (<=0.5 / <=0.3)
#   대조 : 같은 하플로에서 무작위로 고른 같은 수의 자리
import io, os, collections, statistics, random

S = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(S, "inv_out")
A = os.path.join(os.environ["RDNA_V11"], "09_tests", "unit_indel_all_v2_20260916", "anchor_check_20260918")
B = os.path.join(os.environ["RDNA_V11"], "09_tests", "unit_indel_all_v2_20260916")
KRES = os.path.join(os.environ["RDNA_HANDOFF_AD6468DD"], "result_")
TAB = chr(9); NL = chr(10)
def rows(p):
    f = io.open(p, encoding="utf8"); h = f.readline().rstrip(NL).split(TAB)
    for L in f:
        if L.strip(): yield dict(zip(h, L.rstrip(NL).split(TAB)))
RATE = {}
for r in rows(os.path.join(OUT, "ALIGNED_POS_RATE.tsv")):
    RATE[(r["sample"], r["slot"], int(r["array_pos"]))] = (float(r["local_best_rate"]), int(r["n_units"]))
byhap = collections.defaultdict(list)
for (s, sl, ap), v in RATE.items(): byhap[(s, sl)].append((ap, v[0]))

def check(sets, nm, pad=2):
    hit50 = hit30 = seen = nocov = 0; vals = []
    for (s, sl, lo, hi) in sets:
        aps = [ap for ap in range(lo - pad, hi + pad + 1) if (s, sl, ap) in RATE]
        if not aps: nocov += 1; continue
        seen += 1
        m = min(RATE[(s, sl, ap)][0] for ap in aps); vals.append(m)
        if m <= 0.5: hit50 += 1
        if m <= 0.3: hit30 += 1
    # 대조 : 같은 하플로에서 무작위 자리
    random.seed(13); c50 = c30 = cn = 0; cvals = []
    for (s, sl, lo, hi) in sets:
        pool = [ap for ap, _ in byhap.get((s, sl), [])]
        if len(pool) < 5: continue
        w = hi - lo + 1 + 2 * pad
        st = random.choice(pool)
        aps = [ap for ap in range(st, st + w) if (s, sl, ap) in RATE]
        if not aps: continue
        cn += 1; m = min(RATE[(s, sl, ap)][0] for ap in aps); cvals.append(m)
        if m <= 0.5: c50 += 1
        if m <= 0.3: c30 += 1
    return (nm, len(sets), seen, nocov, hit50, hit30, statistics.median(vals) if vals else float("nan"),
            cn, c50, c30, statistics.median(cvals) if cvals else float("nan"))

S1 = [(r["sample"], r["slot"], int(r["ap_lo"]), int(r["ap_hi"])) for r in rows(os.path.join(A, "results", "ANCHOR_SITES.tsv"))
      if int(r["n_reads"]) >= 3 and float(r["support_frac"]) >= 0.7]
S2 = [(r["sample"], r["slot"], int(r["ap_lo_union"]), int(r["ap_hi_union"])) for r in rows(os.path.join(KRES, "SITES_ONT.tsv"))]
S3 = [(r["sample"], r["slot"], int(r["canon_ap"]), int(r["canon_ap"])) for r in rows(os.path.join(B, "UINDEL_SITES.tsv"))
      if r["class"] == "정본의심"]
L = ["# 새 방법으로 기존 삽입/결실을 다시 찾나 (2026-09-30)", "",
     "- 새 방법 지표 = **배치된 자리의 국소 최고율**(±5 이웃 중 최고인 read 유닛 비율). 구간(±2) 안 최소값을 본다.",
     "- 대조 = 같은 하플로에서 같은 폭의 무작위 구간", ""]
L.append("| 묶음 | 자리 | 볼 수 있음 | 깊이없음 | **<=0.5** | <=0.3 | 최소값 중앙 | 대조 n | 대조 <=0.5 | 대조 <=0.3 | 대조 중앙 |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|")
for sets, nm in ((S1, "정본 오류 의심 13(지지>=0.7)"), (S2, "주 값 ONT 자리 73"), (S3, "v2 DP 정본의심")):
    r = check(sets, nm)
    L.append("| %s | %d | %d | %d | **%d (%.0f%%)** | %d (%.0f%%) | %.2f | %d | %d (%.0f%%) | %d (%.0f%%) | %.2f |" % (
        r[0], r[1], r[2], r[3], r[4], 100.0 * r[4] / max(1, r[2]), r[5], 100.0 * r[5] / max(1, r[2]), r[6],
        r[7], r[8], 100.0 * r[8] / max(1, r[7]), r[9], 100.0 * r[9] / max(1, r[7]), r[10]))
L.append("")
L.append("## 정본 오류 의심 13 자리 하나씩")
L.append("| 표본 | 슬롯 | 자리 | 종류 | 구간 국소최고율 최소 | 그 하플로 중앙 |")
L.append("|---|---|---|---|---|---|")
for r in rows(os.path.join(A, "results", "ANCHOR_SITES.tsv")):
    if not (int(r["n_reads"]) >= 3 and float(r["support_frac"]) >= 0.7): continue
    s, sl, lo, hi = r["sample"], r["slot"], int(r["ap_lo"]), int(r["ap_hi"])
    aps = [ap for ap in range(lo - 2, hi + 3) if (s, sl, ap) in RATE]
    med = statistics.median([v for _, v in byhap.get((s, sl), [(0, 0)])])
    L.append("| %s | h%s | %d–%d | %s%s | %s | %.2f |" % (
        s, sl, lo, hi, r["type"], r["size_units"],
        ("%.2f" % min(RATE[(s, sl, ap)][0] for ap in aps)) if aps else "깊이없음", med))
io.open(os.path.join(OUT, "RECOVER_INDELS.txt"), "w", encoding="utf8", newline=NL).write(NL.join(L) + NL)
print(NL.join(L))
