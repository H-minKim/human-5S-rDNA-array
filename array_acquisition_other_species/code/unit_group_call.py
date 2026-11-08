# -*- coding: utf-8 -*-
# array 를 유닛으로 나눌 때 유닛마다 먼저 그룹을 따진다 — 유닛 한 줄짜리 판정표
#
# 유닛마다 서로 독립인 두 증거를 나란히 놓는다
#   증거 A (변이 사건) : 그룹을 가르는 사건들을 이 유닛이 가졌는가 → 비트 패턴 · ev_group
#   증거 B (정렬 점수) : 그룹 참조 둘에 각각 정렬한 점수 차 → aln_group · margin_frac
#   두 증거가 같으면 확정, 다르거나 약하면 미정 으로 두고 array 참조(공통 프레임)를 쓴다
#
# 두 그룹이 진짜인지도 이 표로 검사한다
#   진짜 두 종류면 점수 차 d = score(grp1) - score(grp2) 가 이봉이고 0 근처가 비어 있다.
#   연속체(사실은 한 종류)면 d 가 0 을 중심으로 한 덩어리다.
#   그래서 |d| < 0.2*Dmax 인 '가운데 구역' 유닛 비율을 센다 — 높으면 두 그룹 주장을 접어야 한다.
#   덤 : 사건 패턴이 두 그룹 패턴 중 어느 것도 아닌(섞인) 유닛 수도 센다.
#
# 인자 : array_id 지도(g8mr.tsv) 그룹멤버(members.tsv) 종키 출력폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import re
import sys
import csv
import collections

NL = chr(10)
TAB = chr(9)
AID, MAPF, MEMF, SP, OUT = sys.argv[1:6]
os.makedirs(OUT, exist_ok=True)
NUM = re.compile(r"^([0-9]+)")
pos_of = lambda v: int(NUM.match(v).group(1)) if NUM.match(v) else -1
EVGAP = 20


def events(vs):
    ps = sorted(set(pos_of(v) for v in vs if pos_of(v) >= 0))
    ev, cur = [], []
    for p in ps:
        if cur and p - cur[-1] > EVGAP:
            ev.append(cur)
            cur = []
        cur.append(p)
    if cur:
        ev.append(cur)
    return ev


def comps(cols, n, thr=0.8):
    m = len(cols)
    adj = collections.defaultdict(set)
    for i in range(m):
        a, na = cols[i], len(cols[i])
        for j in range(i + 1, m):
            b, nb = cols[j], len(cols[j])
            num = n * len(a & b) - na * nb
            den = (na * (n - na) * nb * (n - nb)) ** 0.5
            if den > 0 and abs(num / den) >= thr:
                adj[i].add(j)
                adj[j].add(i)
    seen, out = set(), []
    for i in range(m):
        if i in seen or i not in adj:
            continue
        st, c = [i], []
        seen.add(i)
        while st:
            x = st.pop()
            c.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        out.append(sorted(c))
    return out


rows = list(csv.DictReader(io.open(MAPF, encoding="utf8"), delimiter=TAB))
ok = [r for r in rows if r["verify"] == "OK"]
n = len(ok)
# 사건은 공통 프레임(array 참조) 변이로 뽑는다 — 그래야 유닛끼리 자리가 같다
V = [set() if r["varlist_arrayref"] in ("-", "NOALN", "") else set(r["varlist_arrayref"].split(","))
     for r in ok]
cnt = collections.Counter()
for s in V:
    cnt.update(s)
poly = sorted([v for v, c in cnt.items() if 0.10 * n <= c <= 0.90 * n], key=pos_of)
cols = [set(i for i in range(n) if v in V[i]) for v in poly]
bl = [c for c in comps(cols, n) if len(events([poly[k] for k in c])) >= 3]
bl.sort(key=lambda c: -len(events([poly[k] for k in c])))
evv = []
if bl:
    bv = [poly[k] for k in bl[0]]
    evv = [[v for v in bv if pos_of(v) in set(e)] for e in events(bv)]
key = [tuple(1 if sum(1 for v in e if v in V[i]) * 2 >= len(e) else 0 for e in evv) for i in range(n)]
pat = collections.Counter(key)
top = [p for p, c in pat.most_common() if c >= 5][:2]
MEM = dict((r["array_pos"], r["group"]) for r in
           csv.DictReader(io.open(MEMF, encoding="utf8"), delimiter=TAB))

HDR = ["species_key", "array_id", "array_pos", "chrom_start", "chrom_end", "unit_len",
       "ev_pattern", "ev_group", "ev_mismatch_bits", "score_grp1", "score_grp2", "score_diff",
       "margin_frac", "aln_group", "aln_confidence", "final_group", "agree", "ref_to_use"]
o = io.open(os.path.join(OUT, "UNIT_GROUP_CALL_%s.tsv" % AID), "w", encoding="utf8", newline=NL)
o.write(TAB.join(HDR) + NL)
mid = collections.Counter()
fin = collections.Counter()
diffs = []
for i, r in enumerate(ok):
    k = key[i]
    pstr = "".join(str(x) for x in k)
    if top and k == top[0]:
        evg, mmb = "grp1", 0
    elif len(top) > 1 and k == top[1]:
        evg, mmb = "grp2", 0
    elif top:
        d0 = sum(1 for a, b in zip(k, top[0]) if a != b)
        d1 = sum(1 for a, b in zip(k, top[1]) if a != b) if len(top) > 1 else 10 ** 6
        evg = "grp1" if d0 < d1 else ("grp2" if d1 < d0 else "-")
        mmb = min(d0, d1)
    else:
        evg, mmb = "-", -1
    s1 = int(r.get("score_grp1", 0) or 0)
    s2 = int(r.get("score_grp2", 0) or 0)
    d = s1 - s2
    diffs.append(d)
    alng = r["ref_best_raw"]
    conf = r["assign"]
    mid[("가운데" if conf == "AMBIG" else "바깥")] += 1
    agree = (evg == alng and evg in ("grp1", "grp2"))
    if agree and conf in ("OK", "WEAK"):
        fg, use = evg, evg
        fin["확정 " + evg] += 1
    elif conf == "AMBIG" or evg == "-":
        fg, use = "미정", "arrayref"
        fin["미정"] += 1
    else:
        fg, use = "불일치", "arrayref"
        fin["불일치"] += 1
    o.write(TAB.join([SP, AID, r["array_pos"], r["chrom_start"], r["chrom_end"], r["unit_len"],
                      pstr, evg, str(mmb), str(s1), str(s2), str(d), r["margin_frac"], alng, conf,
                      fg, "예" if agree else "아니오", use]) + NL)
o.close()
ds = sorted(diffs)
mm = max(abs(x) for x in ds) if ds else 1
nmid = sum(1 for x in ds if abs(x) < 0.2 * mm)
purep = sum(1 for i in range(n) if top and (key[i] == top[0] or (len(top) > 1 and key[i] == top[1])))
print("%-46s 유닛 %3d · 사건 %2d · 패턴순수 %3d(%4.1f%%) · 점수차 가운데 %3d(%4.1f%%) · %s"
      % (AID[:46], n, len(evv), purep, 100.0 * purep / n, nmid, 100.0 * nmid / n,
         " ".join("%s=%d" % (k, v) for k, v in sorted(fin.items()))))
hist = collections.Counter()
for x in ds:
    hist[int(round(10.0 * x / mm))] += 1
print("      점수차 분포(-10=grp2 쪽 ... +10=grp1 쪽) : "
      + " ".join("%d:%d" % (b, hist[b]) for b in sorted(hist)))
