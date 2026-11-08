# -*- coding: utf-8 -*-
# 정정판 : d 는 read 가 덮은 서로 다른 array_pos 의 개수로 센다
#   (앞 판은 AO/JX·cfg 가 겹쳐 같은 read 를 여러 번 세어 d/c 가 2~4 까지 나왔다)
import io, os, gzip, glob, collections, statistics

# -- (repo) paths come from the environment --------------------------------
#   기본값이 분석 때 쓴 경로와 같다 ⇒ 작업 폴더에서는 전과 똑같이 돈다.
#   그림 내용은 한 글자도 바뀜이 없다 — png 화소로 맞대어 확인한다.
V = os.environ["V_ROOT"]

MAP = os.path.join(V, r"03_maps\g8canonv10\CONSENSUS_MAP_G8_full.tsv")
HSN = os.path.join(V, r"02_canonical\hap_seq_node_v8.tsv")
RAW = os.path.join(V, r"08_methylation\raw")
# 원본은 자기 스크립트 폴더에 썼다. 이 사본은 `analysis/` 에 있으므로 그대로 두면
#   코드 저장소 안에 자료 파일이 생긴다(이 저장소는 코드만 담는다).
#   ⇒ `DC_OUT` 으로 자리를 가리킬 수 있게 한다. 기본값은 전과 같다(스크립트 옆).
OUT = os.environ.get("DC_OUT",
                     os.path.join(os.path.dirname(os.path.abspath(__file__)), "DC_TABLE.tsv"))

c = collections.Counter()
f = io.open(MAP, encoding="utf8")
h = f.readline().rstrip("\n").split("\t")
ix = {k: i for i, k in enumerate(h)}
for L in f:
    r = L.rstrip("\n").split("\t")
    if r[ix["seg_type"]] == "UNIT":
        c[(r[0], r[1])] += 1

BODY = {"arr": "HPRC", "hiarr": "hifiasm", "t1": "span", "c2seq": "cons2", "pat": "T2T"}
src = {}
for L in io.open(HSN, encoding="utf8"):
    r = L.rstrip("\n").split("\t")
    if len(r) >= 3:
        src[(r[0], r[1])] = BODY.get(r[2].split("__")[-1].split("_")[0], "?")

pos = collections.defaultdict(lambda: collections.defaultdict(set))   # (smp,slot) -> read -> {array_pos}
kinds = collections.defaultdict(lambda: collections.defaultdict(set))  # read -> {AO/JX}
for path in sorted(glob.glob(os.path.join(RAW, "METH_READUNIT_*.tsv.gz"))):
    tag = "JX" if "_jx" in os.path.basename(path) else "AO"
    with gzip.open(path, "rt", encoding="utf8") as g:
        hh = g.readline().rstrip("\n").split("\t")
        j = {k: i for i, k in enumerate(hh)}
        for L in g:
            r = L.rstrip("\n").split("\t")
            k = (r[0], r[j["slot"]])
            pos[k][r[j["read_id"]]].add(int(r[j["array_pos"]]))
            kinds[k][r[j["read_id"]]].add(tag)
    print("  읽음 %s" % os.path.basename(path), flush=True)

rows = []
for k in sorted(c):
    if k not in pos or not pos[k]:
        continue
    best = max(pos[k].items(), key=lambda kv: len(kv[1]))
    d = len(best[1])
    rows.append((k[0], k[1], src.get(k, "?"), c[k], d, d / float(c[k]), len(pos[k]),
                 "+".join(sorted(kinds[k][best[0]]))))

with io.open(OUT, "w", encoding="utf8", newline="\n") as g:
    g.write("sample\tslot\tsource\tc_units\td_longest_read_units\td_over_c\tn_reads\tlongest_read_kind\n")
    for r in rows:
        g.write("%s\t%s\t%s\t%d\t%d\t%.4f\t%d\t%s\n" % r)
print("표 %d 행 → %s" % (len(rows), OUT))

print()
print("## 출처별 d/c  (d = 가장 긴 read 가 덮은 서로 다른 유닛 수)")
print("| 출처 | 하플로타입 | c 중앙 | d 중앙 | d/c 중앙 | d/c 최소 | d/c 최대 | d/c>1 |")
print("|---|---:|---:|---:|---:|---:|---:|---:|")
by = collections.defaultdict(list)
for smp, sl, s, cc, dd, q, nr, kd in rows:
    by[s].append((cc, dd, q, smp, sl))
for s in sorted(by, key=lambda x: -len(by[x])):
    v = by[s]
    qs = sorted(x[2] for x in v)
    print("| %s | %d | %d | %d | %.3f | **%.3f** | %.3f | %d |"
          % (s, len(v), statistics.median([x[0] for x in v]), statistics.median([x[1] for x in v]),
             statistics.median(qs), qs[0], qs[-1], sum(1 for x in qs if x > 1.0)))

hi = by.get("hifiasm", [])
mn = min(hi, key=lambda x: x[2])
thr = mn[2]
print()
print("hifiasm 하플로타입 %d 개 · 최소 d/c = %.4f  (%s 슬롯 %s · c=%d · d=%d)"
      % (len(hi), thr, mn[3], mn[4], mn[0], mn[1]))
qs = sorted(x[2] for x in hi)
print("   hifiasm d/c 분포 : q10 %.3f · q25 %.3f · 중앙 %.3f · q75 %.3f · 최대 %.3f"
      % (qs[int(len(qs) * .1)], qs[int(len(qs) * .25)], statistics.median(qs), qs[int(len(qs) * .75)], qs[-1]))

print()
print("## 문턱별 제거 대상 read (전 하플로타입)")
tot = sum(len(pos[k]) for k in pos)
print("| 문턱 d/c | 제거 read | 비율 | 영향 하플로타입 |")
print("|---|---:|---:|---:|")
for t2 in (thr, 0.20, 0.30, 0.50, 0.80, 1.00):
    rem = 0; aff = 0
    for k in pos:
        if k not in c or not c[k]:
            continue
        r2 = sum(1 for u in pos[k].values() if len(u) / float(c[k]) > t2)
        rem += r2
        if r2:
            aff += 1
    print("| %.4f%s | %s | %.2f%% | %d / %d |"
          % (t2, " 최소" if abs(t2 - thr) < 1e-9 else "", format(rem, ","), 100.0 * rem / tot, aff, len(pos)))
print("read 합계 %s" % format(tot, ","))
