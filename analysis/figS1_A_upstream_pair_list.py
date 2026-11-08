# -*- coding: utf-8 -*-
# → pubspan2_pairs.tsv
#   단계 : S1 A 상류
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 짝 목록 v3 — 사용자 규칙 두 개를 넣는다
#   ① 관통 read 는 유닛을 덜 셀 수 있다 (blastn 이 유전자 HSP 를 놓치면 그 유닛이 사라진다)
#      그래서 유닛 수만 보지 말고 길이 ÷ 유닛길이(2,240) 로도 대본다
#   ② 관통 후보가 여럿이면 가장 잘 맞는 것을 쓴다 (±1 안에서만 고르지 않는다)
#   짝짓기 점수 = |길이차(유닛 환산)| 을 1순위, |유닛수차| 를 2순위
#   정본 유닛 수도 조립본 fasta 를 units_of() 로 직접 센 값을 쓴다 (관통과 같은 방법)
import io, os, collections

V = os.environ["RDNA_V11"]
S = os.path.dirname(os.path.abspath(__file__))
UNIT = 2240.0

# 지도 UNIT 개수 (비교용)
cnt = collections.Counter()
f = io.open(os.path.join(V, r"03_maps\g8canonv10\CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
for x in f:
    r = x.rstrip("\n").split("\t")
    if r[ix["seg_type"]] == "UNIT":
        cnt[(r[0], r[1])] += 1

# 정본 슬롯 -> 파일
src = {}
for x in io.open(os.path.join(V, r"02_canonical\hap_seq_node_v8.tsv"), encoding="utf8"):
    r = x.rstrip("\n").split("\t")
    if len(r) >= 3:
        src[(r[0], r[1])] = r[2]


def rd(fn):
    d = {}
    for x in io.open(os.path.join(S, fn), encoding="utf8"):
        a, n, l = x.rstrip("\n").split("\t")
        d[a] = (int(n), int(l))
    return d


measured = rd("canon_units.tsv")        # 정본 조립본 : 유닛수, bp
cand = rd("cmp4d_units.tsv")            # 관통 후보  : 유닛수, bp

body = lambda fn: fn.split("__")[-1]


def kind(fn):
    b = body(fn)
    for t in ("hiarr_", "arr_", "t1_", "c2seq_", "pat"):
        if b.startswith(t):
            return t.rstrip("_")
    return "?"


span = collections.defaultdict(list)
for fn, (n, l) in cand.items():
    if not body(fn).startswith("t1_"):
        continue
    span[body(fn).split("_")[1]].append((fn, n, l))

# 지도 vs 조립본 실측 차이
dif = collections.Counter()
for (smp, slot), fn in src.items():
    a = cnt.get((smp, slot), 0); b = measured.get(fn, (0, 0))[0]
    if a and b:
        dif[b - a] += 1
print("지도 UNIT 수 vs 조립본 실측(units_of) 차이 : " +
      " · ".join("%+d : %d 슬롯" % (k, v) for k, v in sorted(dif.items())))
print()

K = collections.Counter(); pairs = []; drop = []
for (smp, slot), fn in sorted(src.items()):
    kd = kind(fn)
    K["정본_" + kd] += 1
    if kd not in ("arr", "hiarr"):
        drop.append((smp, slot, kd, "정본출처_관통계열", fn)); continue
    if fn not in measured:
        drop.append((smp, slot, kd, "정본파일_노드에없음", fn)); continue
    nu, lu = measured[fn]
    K["대상슬롯"] += 1
    if smp not in span:
        drop.append((smp, slot, kd, "관통후보없음", "")); continue
    K["관통있는개체의슬롯"] += 1
    pairs.append((smp, slot, kd, nu, lu, fn))

by = collections.defaultdict(list)
for p in pairs:
    by[p[0]].append(p)

out = []
# 개체 안에서 (슬롯 x 후보) 를 모두 만들고 잘 맞는 것부터 1:1 로 붙인다.
#   강제로 붙이지 않는다 — 문턱(길이 1 유닛 이내 또는 유닛수 +-1) 을 넘으면 안 붙이고 사유를 남긴다
#   (후보 1 개인 개체의 두 번째 슬롯이 엉뚱한 후보를 먹는 것을 막는다)
for smp, ps in by.items():
    combos = []
    for smp2, slot, kd, nu, lu, fn in ps:
        for c in span[smp]:
            combos.append((abs((c[2] - lu) / UNIT), abs(c[1] - nu), slot, kd, nu, lu, fn, c))
    combos.sort()
    usedS, usedC = set(), set()
    for dl, du, slot, kd, nu, lu, fn, c in combos:
        if slot in usedS or c[0] in usedC:
            continue
        if dl > 1.0 and du > 1:
            continue
        usedS.add(slot); usedC.add(c[0])
        out.append((smp, slot, kd, nu, lu, fn, c[0], c[1], c[2],
                    c[1] - nu, round((c[2] - lu) / UNIT, 3)))
    for smp2, slot, kd, nu, lu, fn in ps:
        if slot in usedS:
            continue
        rest = [c for c in span[smp] if c[0] not in usedC]
        if not rest:
            drop.append((smp, slot, kd, "관통후보_이미다씀", "정본%d유닛" % nu))
        else:
            b = min(rest, key=lambda c: abs((c[2] - lu) / UNIT))
            drop.append((smp, slot, kd, "안맞음(문턱밖)",
                         "정본%d유닛 %dbp / 가장가까운후보 %d유닛 %dbp (길이차 %.2f유닛)"
                         % (nu, lu, b[1], b[2], (b[2] - lu) / UNIT)))

with io.open(os.path.join(S, "pubspan2_pairs.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("sample\tslot\tsrc\tunits\tasm_bp\tpubfile\tspanfile\tspan_units\tspan_bp\td_unit\td_len_unit\n")
    for r in sorted(out):
        g.write("\t".join(map(str, r)) + "\n")
with io.open(os.path.join(S, "pubspan2_drop.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("sample\tslot\tsrc\treason\tdetail\n")
    for r in sorted(drop):
        g.write("\t".join(map(str, r)) + "\n")

print("정본 슬롯 출처 : " + " · ".join("%s %d" % (k, v) for k, v in sorted(K.items())))
print("관통 후보 %d 개 · 개체 %d" % (sum(len(v) for v in span.values()), len(span)))
print()
print("짝 %d 슬롯 · 개체 %d   (출처별 : %s)"
      % (len(out), len(set(x[0] for x in out)),
         " · ".join("%s %d" % (k, v) for k, v in sorted(collections.Counter(x[2] for x in out).items()))))
print()
print("## 검정 ① 관통 read 와 array 유닛 수가 맞는가")
print("| 지표 | 정확 0 | ±1 | ±2~3 | 4 이상 |")
print("|---|---:|---:|---:|---:|")
for nm, gi in (("유닛 수 차이", 9), ("길이 차이(유닛 환산)", 10)):
    b = collections.Counter()
    for r in out:
        d = abs(r[gi])
        b["0" if d < 0.5 else "1" if d < 1.5 else "3" if d < 3.5 else "4"] += 1
    t = float(len(out))
    print("| %s | %d (%.1f%%) | %d (%.1f%%) | %d (%.1f%%) | %d (%.1f%%) |"
          % (nm, b["0"], 100 * b["0"] / t, b["1"], 100 * b["1"] / t,
             b["3"], 100 * b["3"] / t, b["4"], 100 * b["4"] / t))
print()
print("탈락 사유")
for k, v in collections.Counter(x[3] for x in drop).most_common():
    print("   %-24s %d" % (k, v))
