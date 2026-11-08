# -*- coding: utf-8 -*-
# 개체당 2 개 확정 → hap_final2.tsv (③ 단계가 읽는다)
#   단계 : ① 하플로타입 확정

# -*- coding: utf-8 -*-
# 1단계 v2. 하플로타입 확정 — 개체당 정확히 2 개 (2026-09-05)
#
#  v1 에서 실측으로 찾아낸 결함 세 가지
#
#  ① 신규추출 44 개체가 통째로 빠졌다 (270 개체만 나왔다).
#     원인은 "1kGP 필터" 가 아니다 — 실측하니 ex44 ∩ L_kgp = NA18879 단 1 명이다.
#     진짜 원인은 이들의 공개조립본 표 파일명이 `asm44_all.tsv` 로 달라
#     v1 이 그 파일을 아예 읽지 않은 것이다.  -> 270 + 44 = 314
#     ※ 이 44 개체는 공개조립본만 있다 (hifiasm·관통·재조합관통 자료 0).
#
#  ② 유닛 수 프레임이 두 개였다 (인접 군집이 사방에서 갈라진 원인).
#     원표(pub_all·asm44_all·asm_units_all·hisum_v2·refseq_headers)의 유닛 컬럼은
#     "2,100-2,400 bp 정상 간격의 개수" 다.  cmp4 는 "간격 + 1 = 진짜 유닛 수" 다.
#     실측 대조 : pub 595/610 · hifi 708/728 · 관통 178/210 이 정확히 +1 차이.
#     ⇒ 원표 값은 전부 GAPS 프레임이므로 +1 해서 UNITS 프레임으로 바꾼다.
#       변환은 이 파일 안 add_units() 한 곳에서만 한다 (프레임 사고 재발 방지).
#
#  ③ 관통 read 계수 잡음을 두 번째 하플로타입으로 오인했다.
#     한 개체 안 강한 관통 read 들의 인접 간격 분포를 재니 뚜렷한 이봉이다:
#       간격 0:194  1:21  2:8  ||  3:3  4:2  5:1  6:2 ... (평평한 꼬리)
#     간격 <=2 는 같은 분자의 계수 잡음, >=3 부터가 진짜 다른 하플로타입이다.
#     ⇒ 강한 관통 read 는 간격 2 이하 단일연결로 먼저 묶고 중앙값 하나로 쓴다.
#
#  정정된 우선순위 ("관통 read 가 있다면 관통 read 가 정답이야")
#    관통 read 는 한 분자가 flank-to-flank 를 물리적으로 지난 직접 증거다.
#    다만 관통 단독도 틀리는 일이 있으므로, 다른 방법의 지지를 받을 때 최우선으로 둔다.
#      1) 강한 관통(양쪽 flank 확인) + 다른 방법 지지   -> 최우선
#      2) 공개조립본
#      3) 강한 관통 단독
#      4) 그 외 — 교차지지 많은 것 -> hifiasm -> 관통 -> 재조합관통
#
#    공개조립본을 절대 우선으로 두지 않는 근거: 공개조립본이 밀린 12 건을 실측하니
#    관통 read 는 전부 반대편을 지지했고(공개조립본 지지 0 건),
#    HG00128 의 공개조립본 81 유닛 array 에는 44,784 bp 구멍이 있었다.
#
#  강한 관통의 정의 (verify3_all.tsv)
#    kind=="SPAN" · PROX50k>0 · DIST50k>0 (양쪽 flank) · 사슬연속 >=95% · 정상간격 >=19
import io
import os
import re
import glob
import collections

S = os.environ["RDNA_S"]
D = os.environ["RDNA_ASM_COMPARE"]
U = os.environ["RDNA_UNIFIED_REFMAP"]
NL = "\n"
TOL = 1          # 같은 하플로타입으로 볼 유닛 차이 (군집)
SNOISE = 2       # 관통 read 계수 잡음 폭 (실측 이봉 경계)
PRI = {"공개조립본": 0, "우리hifiasm": 1, "관통": 2, "재조합관통": 3}

UNITS = collections.defaultdict(lambda: collections.defaultdict(set))
PUBCTG = collections.defaultdict(dict)      # smp -> contig -> 유닛수 (조립본 contig 단위)


def add_units(smp, method, value, frame):
    """유닛 수 프레임 변환 — 이 함수 한 곳에서만 한다.
       frame='GAPS'  : 정상 간격 수      -> 유닛 = 간격 + 1
       frame='UNITS' : 이미 유닛 수      -> 그대로"""
    UNITS[smp][method].add(value + 1 if frame == "GAPS" else value)


# ㉮ 공개조립본 — asm44_all.tsv 포함 (v1 누락분) · 프레임 GAPS
for p in (os.path.join(D, "asm_units_all.tsv"),
          os.path.join(S, "pub_all.tsv"),
          os.path.join(S, "asm44_all.tsv")):
    if not os.path.exists(p):
        continue
    for L in io.open(p, encoding="utf8"):
        f = L.rstrip(NL).split("\t")
        if len(f) >= 9 and f[8] == "OK":
            try:
                add_units(f[0], "공개조립본", int(f[7]), "GAPS")
                PUBCTG[f[0]][f[3]] = int(f[7]) + 1      # contig 별로도 남긴다
            except ValueError:
                pass

# ㉯ 우리 hifiasm — 프레임 GAPS
for L in io.open(os.path.join(S, "hisum_v2.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) >= 6 and f[1] == "OK":
        s = re.sub(r"_\d+$", "", f[0])                    # _N 을 먼저 뗀다
        for pat in (r"_(sw|f)[2-6]$", r"_l0$", r"_sub250$", r"_noSPAN$"):
            s = re.sub(pat, "", s)
        try:
            add_units(s, "우리hifiasm", int(f[5]), "GAPS")
        except ValueError:
            pass

# ㉰ 관통 + ㉱ 재조합관통 (옛 refseq 헤더) — 프레임 GAPS
for L in io.open(os.path.join(S, "refseq_headers.txt"), encoding="utf8"):
    f = L.strip().split("|")
    if len(f) < 4:
        continue
    m = re.match(r"^N(\d+)$", f[3])
    if not m:
        continue
    k = "관통" if f[2] == "tier1" else ("재조합관통" if f[2] == "tier2" else None)
    if k:
        add_units(f[0], k, int(m.group(1)), "GAPS")

# cmp4 가 서열에서 직접 잰 값 — 프레임 UNITS
KN = {"pub": "공개조립본", "hifi": "우리hifiasm", "span": "관통"}
for L in io.open(os.path.join(S, "cmp4_all.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) < 9:
        continue
    for col, kk in ((3, f[1]), (4, f[2])):
        if kk not in KN or f[col] in ("-", ""):
            continue
        try:
            add_units(f[0], KN[kk], int(f[col]), "UNITS")
        except ValueError:
            pass

# ---- 강한 관통 (양쪽 flank 확인) — 프레임 GAPS, 잡음 병합 ----
RAW = collections.defaultdict(list)
for L in io.open(os.path.join(S, "verify3_all.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) < 10 or f[1] != "SPAN":
        continue
    try:
        px, dx, frac, gaps = int(f[4]), int(f[5]), float(f[8]), int(f[9])
    except ValueError:
        continue
    if px > 0 and dx > 0 and frac >= 95.0 and gaps >= 19:
        RAW[f[0]].append(gaps + 1)

STRONG = {}
NREAD = {}
for s, vs in RAW.items():
    vs.sort()
    grp = [[vs[0]]]
    for v in vs[1:]:
        if v - grp[-1][-1] <= SNOISE:      # 간격 <=2 는 같은 분자의 계수 잡음
            grp[-1].append(v)
        else:
            grp.append([v])
    reps = [g[len(g) // 2] for g in grp]
    STRONG[s] = set(reps)
    NREAD[s] = dict((g[len(g) // 2], len(g)) for g in grp)
    for r in reps:
        add_units(s, "관통", r, "UNITS")

# ---- 코호트 ----
hprc = set()
for p in glob.glob(os.path.join(U, "refmap_hap_ip-*.tsv")):
    for i, L in enumerate(io.open(p, encoding="utf8")):
        if i == 0:
            continue
        f = L.rstrip(NL).split("\t")
        if len(f) >= 6:
            hprc.add(f[0])
EX44 = set(x.strip() for x in io.open(os.path.join(S, "ex44.txt"), encoding="utf8") if x.strip())
kgp = set(x.strip() for x in io.open(os.path.join(S, "L_kgp.txt"), encoding="utf8") if x.strip())
kgp -= (hprc | EX44)          # 두 코호트에 걸친 개체는 제외 대상이 아니다
TEST = ("PARAM_", "D5_", "D50_", "D200_")
smps = sorted(s for s in UNITS if s not in kgp and not s.startswith(TEST))

TN = {0: "강한관통+교차지지", 1: "공개조립본", 2: "강한관통 단독", 3: "그 외(교차지지順)"}


def tier_of(c):
    if c["strong"] and c["n"] >= 2:
        return 0
    if "공개조립본" in c["m"]:
        return 1
    if c["strong"]:
        return 2
    return 3


def cluster(smp):
    pts = [(v, k) for k, vs in UNITS[smp].items() for v in vs]
    pts.sort(key=lambda z: -z[0])
    cl = []
    for v, k in pts:
        for c in cl:
            if abs(v - c["u"]) <= TOL:
                c["m"].add(k)
                c["v"].append(v)
                c["u"] = sorted(c["v"])[len(c["v"]) // 2]
                break
        else:
            cl.append({"u": v, "v": [v], "m": set([k])})
    for c in cl:
        c["n"] = len(c["m"])
        c["pri"] = min(PRI[x] for x in c["m"])
        c["strong"] = False
    # 강한 관통 대표값 하나는 가장 가까운 군집 하나에만 붙인다.
    #   허용오차(2)만으로 표식하면 2 유닛 떨어진 두 군집이 둘 다 강한관통이 되어
    #   같은 분자가 하플로타입 두 개를 차지한다 (HG03193 69/67 사례).
    for x in STRONG.get(smp, set()):
        near = min(cl, key=lambda c: abs(c["u"] - x))
        if abs(near["u"] - x) <= SNOISE:
            near["strong"] = True
    cl.sort(key=lambda c: (tier_of(c), -c["n"], c["pri"], -c["u"]))
    return cl


def pick(smp):
    """두 슬롯을 고른다.

    쌍둥이 예외 — 정답지 실측상 진짜 하플로타입 쌍의 4.4% 가 1 유닛 이내다.
    이런 쌍은 군집(TOL=1)이 하나로 합쳐버려 두 번째 슬롯을 엉뚱한 값이 채웠다
    (실측 3 건: HG02040 정답 70,71 -> 71,83 / NA20805 75,76 -> 45,76 /
     HG01258 74,75 -> 72,74).
    공개 조립본이 서로 다른 contig 2 개를 냈다는 것은 분자가 2 개라는 직접 증거이므로,
    두 contig 가 같은 군집에 들어갈 만큼 가까우면 그 두 값으로 두 슬롯을 채운다."""
    cl = cluster(smp)
    v = sorted(PUBCTG.get(smp, {}).values())
    if len(v) == 2 and v[1] - v[0] <= TOL and v[1] != v[0]:
        src = dict((c["u"], c) for c in cl)
        base = min(cl, key=lambda c: abs(c["u"] - v[1])) if cl else None
        if base is not None and abs(base["u"] - v[1]) <= TOL:
            out = []
            for x in (v[1], v[0]):
                d = dict(base)
                d["u"] = x
                d["twin"] = True
                out.append(d)
            return out, cl
    return cl[:2], cl


print("=" * 86)
print("  1단계 v2. 하플로타입 확정 — 개체 %d  (v1 은 270)" % len(smps))
print("=" * 86)
print("  코호트  HPRC ONT %d · 신규추출44 %d · 그 외 %d" %
      (len(set(smps) & hprc), len(set(smps) & EX44), len(set(smps) - hprc - EX44)))
print("  강한 관통(양쪽 flank)을 가진 개체 %d · 그중 서로 다른 관통군 2 개 이상 %d" %
      (len([s for s in smps if STRONG.get(s)]),
       len([s for s in smps if len(STRONG.get(s, ())) >= 2])))
print()

nc = collections.Counter()
sup = collections.Counter()
tiers = collections.Counter()
final = {}
ntwin = 0
for s in smps:
    got, cl = pick(s)
    nc[min(len(cl), 5)] += 1
    final[s] = got
    if any(c.get("twin") for c in got):
        ntwin += 1
    for c in got:
        sup[c["n"]] += 1
        tiers[tier_of(c)] += 1

print("  후보 군집 수 분포 (추리기 전)")
for k in sorted(nc):
    print("     %-3s개 : %3d 개체 (%5.1f%%)" %
          ("5+" if k == 5 else k, nc[k], 100.0 * nc[k] / len(smps)))
print()
n2 = sum(1 for s in smps if len(final[s]) >= 2)
tot = sum(tiers.values())
print("  2 슬롯 채운 개체 : %d / %d (%.1f%%) · 1 슬롯 %d · 0 슬롯 %d" %
      (n2, len(smps), 100.0 * n2 / len(smps),
       sum(1 for s in smps if len(final[s]) == 1),
       sum(1 for s in smps if not final[s])))
print("  쌍둥이 예외로 구제한 개체 : %d" % ntwin)
print("  확정 하플로타입 : %d / %d 슬롯 = %.1f%%" %
      (tot, 2 * len(smps), 100.0 * tot / (2 * len(smps))))
print()
print("  채택 근거")
for k in sorted(tiers):
    print("     %-20s %4d (%5.1f%%)" % (TN[k], tiers[k], 100.0 * tiers[k] / tot))
print()
print("  교차지지 (몇 개 방법이 같은 유닛 수를 지목했나)")
for k in sorted(sup, reverse=True):
    print("     %d 개 방법 : %4d (%5.1f%%)" % (k, sup[k], 100.0 * sup[k] / tot))
print()
print("  단독 지목인 하플로타입의 출처")
one = collections.Counter()
for s in smps:
    for c in final[s]:
        if c["n"] == 1:
            one[list(c["m"])[0]] += 1
for k, v in one.most_common():
    print("     %-12s %d" % (k, v))
print()

# ---- 공개조립본이 밀린 건 ----
print("  === 공개조립본이 밀린 건 (관통이 이긴 경우) ===")
print("  %-11s %-14s %-14s %s" % ("개체", "채택", "버린 공개조립본", "강한관통(read수)"))
print("  " + "-" * 76)
lost = 0
for s in smps:
    got, cl = pick(s)
    if any(c.get("twin") for c in got):
        continue
    if len(cl) <= 2:
        continue
    drop = [c for c in cl[2:] if "공개조립본" in c["m"]]
    if not drop:
        continue
    lost += 1
    print("  %-11s %-14s %-14s %s" %
          (s, ",".join("%d" % c["u"] for c in got),
           ",".join("%d" % c["u"] for c in drop),
           ",".join("%d(%d)" % (x, NREAD[s][x]) for x in sorted(STRONG.get(s, ())))
           or "없음"))
print("  공개조립본이 밀린 개체 : %d" % lost)
print()

# ---- 관통과 hifiasm 이 정면충돌한 건 (사람이 볼 것) ----
print("  === 관통 단독이 hifiasm+재조합관통 을 이긴 건 (검토 필요) ===")
flag = 0
for s in smps:
    got, cl = pick(s)
    if not any(tier_of(c) == 2 for c in got):
        continue
    drop = [c for c in cl[2:] if c["n"] >= 2]
    if not drop:
        continue
    flag += 1
    print("  %-11s 채택 %-12s 버림 %s" %
          (s, ",".join("%d(%s)" % (c["u"], TN[tier_of(c)][:6]) for c in got),
           " · ".join("%d(%s)" % (c["u"], ",".join(sorted(c["m"]))) for c in drop)))
print("  검토 대상 : %d" % flag)

OUT = io.open(os.path.join(S, "hap_final2.tsv"), "w", encoding="utf8", newline=NL)
OUT.write("sample\tslot\tunits\tn_method\tstrong_span\tn_span_read\ttier\tmethods\n")
for s in smps:
    for i, c in enumerate(final[s]):
        nr = 0
        for x, k in NREAD.get(s, {}).items():
            if abs(c["u"] - x) <= SNOISE:
                nr = max(nr, k)
        OUT.write("%s\t%d\t%d\t%d\t%s\t%d\t%s\t%s\n" %
                  (s, i + 1, c["u"], c["n"], "Y" if c["strong"] else "N", nr,
                   TN[tier_of(c)], ",".join(sorted(c["m"]))))
OUT.close()
print()
print("  저장 : hap_final2.tsv")
