# -*- coding: utf-8 -*-
# 단계 ③ 슬롯 배정 (하플로타입마다 파일 1 개)
#   단계 : ③ — hap_final2.tsv + seqidx_all.tsv → hap_seq.tsv → hap_seq_node_v8.tsv
#   이 글이 슬롯 번호를 정한다 — 방법 문서가 *"가장 위험한 곳"* 이라 적은 단계다.

# -*- coding: utf-8 -*-
# 2단계. 확정 하플로타입에 실제 서열을 배정한다 (2026-09-05)
#
#  입력
#    hap_final2.tsv  — 314 개체 × 슬롯 (유닛 수 · 채택 근거)   프레임 UNITS
#    seqidx_all.tsv  — cmp4d 의 .fa 1,516 개 전수 색인          프레임 UNITS
#
#  배정 규칙
#    ① 유닛 수가 확정값과 ±1 안에 드는 .fa 만 후보로 본다.
#       (하플로타입을 고른 근거는 이미 1단계에서 정했다. 여기서는 그 유닛 수를 가진
#        서열 중 가장 정확한 것을 고르는 문제만 남는다.)
#    ② 서열 정확도 순으로 고른다:  공개조립본 > 우리 hifiasm > 관통 read
#       조립본은 다듬어진 서열(~99.9%+)이고 관통 read 는 ONT 원시 read(~99%)다.
#       관통이 유닛 수 싸움에서 이겼더라도, 같은 유닛 수를 내는 조립 서열이 있으면
#         그 서열을 쓴다 (유닛 수는 관통이, 염기는 조립본이 정확하다).
#    ③ 같은 등급 안에서는 사슬연속 % 높은 것 -> bp 큰 것.
#    ④ 한 개체의 두 슬롯에 같은 파일을 배정하지 않는다.
#
#  슬롯은 실제 서열로 뒷받침돼야 한다 (2026-09-05 정정).
#    1 단계가 1·2 위로 뽑은 군집에 대응하는 .fa 가 없으면 다음 순위 군집으로 내려간다.
#    실측으로 이것이 두 건을 바로잡았다:
#      NA19120  91(관통 read 사슬계수) -> 171 : 91 유닛짜리 서열이 아예 없고
#               t1 조립 서열은 96, hifiasm 은 171 이다. 91 은 96 의 사슬 저계수였다.
#      HG01028  45 -> 40 : 같은 관통 read 로 만든 t1 참조 서열이 40 유닛이다.
#    끝까지 서열이 없으면 그 슬롯은 유닛 수만 남기고 서열 없음으로 표시한다
#    (재조합관통 단독 군집 — 합성 서열이라 쓸 수 없다).
#
#  재조합관통(옛 refseq_*.fa)은 합성 서열(정확히 N x 2240)이라 서열로 쓰지 않는다.
import io
import os
import collections

S = os.environ["RDNA_S"]
NL = "\n"
TOL = 1
KPRI = {"pub": 0, "hifi": 1, "span": 2}
KNAME = {"pub": "공개조립본", "hifi": "우리hifiasm", "span": "관통"}

IDX = collections.defaultdict(list)
for i, L in enumerate(io.open(os.path.join(S, "seqidx_all.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 8:
        continue
    IDX[f[1]].append({"file": f[0], "kind": f[2], "dir": f[3],
                      "bp": int(f[4]), "u": int(f[5]), "cp": float(f[6])})

HAP = collections.defaultdict(list)
for i, L in enumerate(io.open(os.path.join(S, "hap_clusters.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    HAP[f[0]].append({"rank": int(f[1]), "u": int(f[2]), "n": int(f[3]),
                      "strong": f[4], "nread": int(f[5]), "tier": f[6],
                      "m": f[7], "twin": f[8]})

OUT = io.open(os.path.join(S, "hap_seq.tsv"), "w", encoding="utf8", newline=NL)
OUT.write("sample\tslot\tunits_final\trank\ttier\tmethods\t"
          "seq_file\tseq_kind\tseq_dir\tseq_bp\tseq_units\tchain_pct\n")
ok = collections.Counter()
miss = []
kcnt = collections.Counter()
demoted = []
for s in sorted(HAP):
    used = set()
    got = []
    left = sorted(HAP[s], key=lambda x: x["rank"])
    # ① 서열이 있는 군집을 순위대로 채운다
    for h in left:
        if len(got) >= 2:
            break
        cand = [c for c in IDX.get(s, [])
                if abs(c["u"] - h["u"]) <= TOL and c["file"] not in used]
        if not cand:
            continue
        cand.sort(key=lambda c: (KPRI[c["kind"]], abs(c["u"] - h["u"]), -c["cp"], -c["bp"]))
        b = cand[0]
        used.add(b["file"])
        got.append((h, b))
        if h["rank"] > 2:
            demoted.append((s, h, b))
    # ② 그래도 2 개가 안 차면 서열 없는 상위 군집을 유닛 수만으로 남긴다
    if len(got) < 2:
        for h in left:
            if len(got) >= 2:
                break
            if any(g[0]["u"] == h["u"] for g in got):
                continue
            got.append((h, None))
    for i, (h, b) in enumerate(got):
        if b is None:
            ok["서열없음"] += 1
            miss.append((s, h))
            OUT.write("%s\t%d\t%d\t%d\t%s\t%s\t-\t-\t-\t-\t-\t-%s"
                      % (s, i + 1, h["u"], h["rank"], h["tier"], h["m"], NL))
        else:
            ok["배정"] += 1
            kcnt[b["kind"]] += 1
            OUT.write("%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.1f%s"
                      % (s, i + 1, b["u"], h["rank"], h["tier"], h["m"], b["file"],
                         KNAME[b["kind"]], b["dir"], b["bp"], b["u"], b["cp"], NL))
OUT.close()

tot = ok["배정"] + ok["서열없음"]
print("=" * 80)
print("  2단계. 서열 배정 — 확정 하플로타입 %d" % tot)
print("=" * 80)
print("  서열 배정됨 : %d (%.1f%%) · 서열 없음 %d" %
      (ok["배정"], 100.0 * ok["배정"] / tot, ok["서열없음"]))
print()
print("  배정된 서열의 출처")
for k, v in kcnt.most_common():
    print("     %-12s %4d (%5.1f%%)" % (KNAME[k], v, 100.0 * v / ok["배정"]))
print()
print("  === 서열이 없어 다음 순위로 내려간 건 ===")
for s, h, b in demoted:
    print("     %-11s 순위%d 채택 %3d 유닛 (%s) — 상위 군집엔 서열이 없었다" %
          (s, h["rank"], b["u"], KNAME[b["kind"]]))
print("  내려간 건 : %d" % len(demoted))
print()
if miss:
    print("  서열 없는 슬롯 %d 건 (유닛 수만 남김 — 재조합관통 단독 등)" % len(miss))
    for s, h in miss:
        print("     %-12s %3d 유닛  [%s]  그 개체의 서열 후보: %s" %
              (s, h["u"], h["m"],
               ",".join("%s%d" % (c["kind"], c["u"]) for c in sorted(
                   IDX.get(s, []), key=lambda c: c["u"])) or "없음"))
