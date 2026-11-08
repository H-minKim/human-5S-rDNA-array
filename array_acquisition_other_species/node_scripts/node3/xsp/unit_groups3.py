# -*- coding: utf-8 -*-
# array 안 유닛이 한 종류인가 여러 종류인가 — 3 판
#   판정 축 : "서로 떨어진 변이 사건이 블록으로 같이 다니는가"
#     · 한 종류   = 변이가 유닛마다 독립적으로 흩어진다 (연관 없음)
#     · 여러 종류 = 떨어진 여러 사건이 같은 유닛들에 함께 나타난다 (강한 연관)
#   2 판에서 고친 것
#     · 2 판은 좁은 구간에 몰린 변이(예 1397~1527 안 27 개 = 발산 구간 하나)를 '변이 27 개'로 세어
#       단일 사건을 블록으로 승격시켰다. 3 판은 연속 자리를 사건 하나로 접는다(간격 <= EVGAP).
#   방법
#     1) 다형성 변이 : 유닛의 10~90% 에 있는 변이
#     2) 변이쌍 연관 phi · |phi| >= 0.8 이면 연결 → 연결 성분
#     3) 성분의 변이 자리를 사건으로 접고(간격 <= 20 bp 는 한 사건) 사건 >= 3 인 성분만 '블록'
#     4) 귀무 대조 : 변이마다 같은 빈도로 유닛을 무작위 재배치해 같은 절차를 NPERM 회
#     5) 위치 응집 : 그룹이 array 위에서 이어져 있는가(구역 분화) vs 섞여 있는가(유닛 종류)
#   유닛 20 개 미만이면 판정하지 않는다
#   인자 : 종키 fasta 지도폴더 출력폴더 [셔플 횟수=50] [EVGAP=20]
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import re
import sys
import csv
import random
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
SP, FA, MAPD, OUT = sys.argv[1:5]
NPERM = int(sys.argv[5]) if len(sys.argv) > 5 else 50
EVGAP = int(sys.argv[6]) if len(sys.argv) > 6 else 20
random.seed(20260913)
os.makedirs(os.path.join(OUT, "groups"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))
NUM = re.compile(r"^([0-9]+)")


def pos_of(v):
    m = NUM.match(v)
    return int(m.group(1)) if m else -1


def n_events(vs):
    # 변이 자리를 정렬해 간격 <= EVGAP 이면 한 사건으로 접는다
    ps = sorted(set(pos_of(v) for v in vs if pos_of(v) >= 0))
    if not ps:
        return 0
    k = 1
    for a, b in zip(ps, ps[1:]):
        if b - a > EVGAP:
            k += 1
    return k


def blocks(cols, poly, n, thr=0.8, minev=3):
    m = len(cols)
    adj = collections.defaultdict(set)
    for i in range(m):
        a = cols[i]
        na = len(a)
        for j in range(i + 1, m):
            b = cols[j]
            nb = len(b)
            num = n * len(a & b) - na * nb
            den = (na * (n - na) * nb * (n - nb)) ** 0.5
            if den > 0 and abs(num / den) >= thr:
                adj[i].add(j)
                adj[j].add(i)
    seen, out = set(), []
    for i in range(m):
        if i in seen or i not in adj:
            continue
        st, comp = [i], []
        seen.add(i)
        while st:
            x = st.pop()
            comp.append(x)
            for y in adj[x]:
                if y not in seen:
                    seen.add(y)
                    st.append(y)
        comp.sort()
        if n_events([poly[k] for k in comp]) >= minev:
            out.append(comp)
    return out


def adjacency_same(lab, order):
    return sum(1 for a, b in zip(order, order[1:]) if lab[a] is not None and lab[a] == lab[b])


HDR = ["species", "array_id", "n_units_OK", "n_var_sites", "n_poly_var", "n_poly_event", "n_blocks",
       "block_var_sizes", "block_event_sizes", "null_blocks_mean", "null_blocks_max", "verdict",
       "n_groups", "group_sizes", "largest_group_pct", "pos_clustering_ratio", "pos_clustering_p",
       "pos_pattern", "block1_vars", "note"]
SUM = io.open(os.path.join(OUT, "UNIT_GROUPS3_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
SUM.write(TAB.join(HDR) + NL)


def w(d):
    SUM.write(TAB.join(str(d.get(k, "-")) for k in HDR) + NL)


for fn in sorted(os.listdir(MAPD)):
    if not fn.endswith(".g8map.tsv"):
        continue
    aid = fn[:-len(".g8map.tsv")]
    rows = list(csv.DictReader(io.open(os.path.join(MAPD, fn), encoding="utf8"), delimiter=TAB))
    ok = [r for r in rows if r["verify"] == "OK"]
    n = len(ok)
    base = {"species": SP, "array_id": aid, "n_units_OK": n}
    if n < 20:
        base["verdict"] = "판정불가(OK 유닛<20)"
        w(base)
        continue
    V = [set() if r["varlist"] in ("-", "NOALN", "REV", "") else set(r["varlist"].split(",")) for r in ok]
    cnt = collections.Counter()
    for s in V:
        cnt.update(s)
    poly = sorted([v for v, c in cnt.items() if 0.10 * n <= c <= 0.90 * n], key=pos_of)
    base.update({"n_var_sites": len(cnt), "n_poly_var": len(poly), "n_poly_event": n_events(poly)})
    if n_events(poly) < 3:
        base.update({"n_blocks": 0, "verdict": "한 종류", "n_groups": 1, "group_sizes": n,
                     "largest_group_pct": "100.0",
                     "note": "다형성 사건 %d 개뿐(변이 %d)" % (n_events(poly), len(poly))})
        w(base)
        print("%s : 유닛 %d · 다형성 변이 %d/사건 %d -> 한 종류" % (aid, n, len(poly), n_events(poly)))
        continue
    cols = [set(i for i in range(n) if v in V[i]) for v in poly]
    bl = blocks(cols, poly, n)
    nullc = []
    allidx = list(range(n))
    for _ in range(NPERM):
        nullc.append(len(blocks([set(random.sample(allidx, len(c))) for c in cols], poly, n)))
    nm = sum(nullc) / float(len(nullc))
    base.update({"n_blocks": len(bl),
                 "block_var_sizes": ",".join(str(len(b)) for b in bl) if bl else "-",
                 "block_event_sizes": ",".join(str(n_events([poly[k] for k in b])) for b in bl) if bl else "-",
                 "null_blocks_mean": "%.2f" % nm, "null_blocks_max": max(nullc)})
    if not bl:
        base.update({"verdict": "한 종류(연관 블록 없음)", "n_groups": 1, "group_sizes": n,
                     "largest_group_pct": "100.0", "note": "변이는 있으나 떨어진 사건끼리 묶이지 않는다"})
        w(base)
        print("%s : 유닛 %d · 다형성 변이 %d/사건 %d · 블록 0 (귀무 %.1f) -> 한 종류"
              % (aid, n, len(poly), n_events(poly), nm))
        continue
    bl.sort(key=lambda b: -n_events([poly[k] for k in b]))
    b1 = bl[0]
    base["block1_vars"] = ";".join(poly[k] for k in b1[:10])
    key = [tuple(1 if poly[k] in V[i] else 0 for k in b1) for i in range(n)]
    pat = collections.Counter(key)
    grp = [(p, [i for i in range(n) if key[i] == p]) for p, c in pat.most_common() if c >= 5]
    sizes = [len(g[1]) for g in grp]
    if len(grp) <= 1:
        base.update({"verdict": "한 종류(블록 있으나 패턴 1)", "n_groups": 1, "group_sizes": n,
                     "largest_group_pct": "100.0"})
        w(base)
        continue
    lab = [None] * n
    for k, (p, mem) in enumerate(grp):
        for i in mem:
            lab[i] = k
    order = sorted(range(n), key=lambda i: int(ok[i]["array_pos"]))
    obs = adjacency_same(lab, order)
    nulls = []
    for _ in range(NPERM):
        sh = lab[:]
        random.shuffle(sh)
        nulls.append(adjacency_same(sh, order))
    em = sum(nulls) / float(len(nulls))
    ratio = obs / em if em > 0 else float("nan")
    pv = (sum(1 for x in nulls if x >= obs) + 1.0) / (len(nulls) + 1.0)
    base.update({"pos_clustering_ratio": "%.2f" % ratio, "pos_clustering_p": "%.3f" % pv,
                 "pos_pattern": "구역 분화(이어져 있음)" if (ratio >= 1.5 and pv <= 0.05) else "섞여 있음"})
    regs, meta = [], []
    for k, (p, mem) in enumerate(grp):
        for i in mem:
            r = ok[i]
            regs.append("%s:%s-%s" % (r["chrom"], r["chrom_start"], r["chrom_end"]))
            meta.append((k, r["strand"]))
    rf = os.path.join(OUT, "groups", aid + ".regions.txt")
    io.open(rf, "w", encoding="utf8", newline=NL).write(NL.join(regs) + NL)
    rr = subprocess.run(["samtools", "faidx", "-r", rf, FA], capture_output=True, text=True)
    seqs, cur = [], []
    for L in rr.stdout.splitlines():
        if L.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(L.strip())
    if cur:
        seqs.append("".join(cur))
    os.remove(rf)
    bg = collections.defaultdict(list)
    for (k, st), s in zip(meta, seqs):
        bg[k].append(rc(s.upper()) if st == "-" else s.upper())
    files = []
    for k in sorted(bg):
        us = bg[k]
        ml = collections.Counter(len(u) for u in us).most_common(1)[0][0]
        pool = [u for u in us if len(u) == ml] or us
        con = "".join(collections.Counter(u[i] for u in pool).most_common(1)[0][0] for i in range(ml))
        p = os.path.join(OUT, "groups", "%s.grp%d.unit.fa" % (aid, k + 1))
        io.open(p, "w", encoding="utf8", newline=NL).write(
            ">%s_grp%d n=%d len=%d%s%s%s" % (aid, k + 1, len(us), ml, NL, con, NL))
        files.append(os.path.basename(p))
    base.update({"verdict": "여러 종류", "n_groups": len(grp),
                 "group_sizes": ",".join(str(x) for x in sizes),
                 "largest_group_pct": "%.1f" % (100.0 * sizes[0] / n),
                 "note": "유닛파일 " + ",".join(files)})
    w(base)
    print("%s : 유닛 %d · 다형성 변이 %d/사건 %d · 블록 사건수 %s (귀무 %.1f) -> 그룹 %d %s · 위치응집 %.2f p=%.3f"
          % (aid, n, len(poly), n_events(poly), [n_events([poly[k] for k in b]) for b in bl], nm,
             len(grp), sizes, ratio, pv))
SUM.close()
