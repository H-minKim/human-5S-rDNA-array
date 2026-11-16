# -*- coding: utf-8 -*-
# → pubspan2_pairs.tsv
#   Step: S1 A upstream.
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# The pair list, with two further rules:
#   1) An array-spanning read can undercount units: if blastn misses a gene HSP that unit disappears,
#      so the pairing also compares length divided by the unit length of 2,240, not the unit count alone.
#   2) Where several spanning reads are candidates, the best-matching one is used, rather than only those within +-1.
#   The pairing score ranks on the absolute length difference in units first, and the absolute unit-count difference second.
#   The canonical unit count is also counted directly from the assembly fasta with units_of(), the same way as for the spanning read.
import io, os, collections

V = os.environ["RDNA_V11"]
S = os.path.dirname(os.path.abspath(__file__))
UNIT = 2240.0

# the unit count from the map, for comparison
cnt = collections.Counter()
f = io.open(os.path.join(V, r"03_maps/g8canonv10/CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
h = f.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
for x in f:
    r = x.rstrip("\n").split("\t")
    if r[ix["seg_type"]] == "UNIT":
        cnt[(r[0], r[1])] += 1

# canonical slot -> file
src = {}
for x in io.open(os.path.join(V, r"02_canonical/hap_seq_node_v8.tsv"), encoding="utf8"):
    r = x.rstrip("\n").split("\t")
    if len(r) >= 3:
        src[(r[0], r[1])] = r[2]


def rd(fn):
    d = {}
    for x in io.open(os.path.join(S, fn), encoding="utf8"):
        a, n, l = x.rstrip("\n").split("\t")
        d[a] = (int(n), int(l))
    return d


measured = rd("canon_units.tsv")        # canonical assembly: unit count, bp
cand = rd("cmp4d_units.tsv")            # spanning candidate: unit count, bp

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

# the measured difference between the map and the assembly
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
# Within an individual every (slot, candidate) pair is formed and matched one to one, best match first.
#   No match is forced: beyond the threshold (within one unit of length, or +-1 in unit count) the pair is left unmatched and the reason recorded,
#   which stops the second slot of an individual with one candidate from taking the wrong one.
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
