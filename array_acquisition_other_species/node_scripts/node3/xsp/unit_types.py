# -*- coding: utf-8 -*-
# array 안 유닛이 한 종류인가 여러 종류인가 — 최종판 (unit_types)
#
# 무엇을 판정하는가
#   사용자 요구 : array 마다 참조 유닛을 두되, array 안이 같은 유닛의 반복인지 먼저 확인하고
#   아니면 보고하고 유닛을 그룹으로 나눈다. 그룹마다 참조 유닛을 따로 만든다.
#
# 판정 축 : "서로 떨어진 변이 사건이 블록으로 같이 다니는가"
#   · 한 종류   = 변이가 유닛마다 독립적으로 흩어진다 (연관 없음)
#   · 여러 종류 = 떨어진 여러 사건이 같은 유닛들에 함께 나타난다 (강한 연관 = 하플로타입)
#
# 이전 판에서 고친 것 (모두 실측으로 드러난 함정)
#   1 판 : 변이 패턴이 다르기만 하면 나눠 사람 array 를 6 그룹으로 쪼갰다 → 연관(phi) 기준으로 교체
#   2 판 : 좁은 구간에 몰린 변이(예 1397~1527 안 27 개)를 '변이 27 개'로 세어 단일 사건을 블록으로
#          승격시켰다 → 연속 자리를 사건 하나로 접는다
#   3 판 : 분할 열쇠가 변이 하나하나라 한 그룹이 잘게 쪼개졌다(수마트라오랑우탄 2·3·5 가 0.3~0.5% 차이)
#          → 분할 열쇠를 사건 단위로
#
# 절차
#   0) 상태표에서 NON_5S / NO_UNIT 인 array 는 판정 대상에서 뺀다
#   1) 다형성 변이 = 유닛의 10~90% 에 있는 변이 · 연속 자리(간격 <= EVGAP)는 사건 하나로 접는다
#   2) 변이쌍 연관 |phi| >= 0.8 로 이어 연결 성분 · 사건 >= 3 인 성분만 '블록'
#   3) 가장 큰 블록의 사건 보유 패턴으로 유닛을 나눈다 (유닛 >= 5 인 패턴만 그룹)
#   4) 귀무 대조 : 변이마다 같은 빈도로 유닛을 무작위 재배치해 같은 절차 NPERM 회
#   5) 홀/짝 교차검증 : 블록 사건을 홀·짝 두 묶음으로 갈라 각각 분할 → ARI 로 재현성
#   6) 위치 응집 : 그룹이 array 위에서 이어져 있는가(구역 분화) vs 섞여 있는가(유닛 종류)
#   7) 그룹마다 자리별 다수결로 그룹 참조 유닛 생성
#   OK 유닛 20 개 미만이면 판정하지 않는다
#
# 인자 : 종키 fasta 지도폴더 출력폴더 상태표.tsv [NPERM=50] [EVGAP=20]
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
SP, FA, MAPD, OUT, STF = sys.argv[1:6]
NPERM = int(sys.argv[6]) if len(sys.argv) > 6 else 50
EVGAP = int(sys.argv[7]) if len(sys.argv) > 7 else 20
random.seed(20260913)
os.makedirs(os.path.join(OUT, "groups"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))
NUM = re.compile(r"^([0-9]+)")
pos_of = lambda v: int(NUM.match(v).group(1)) if NUM.match(v) else -1
STATUS = {}
if os.path.exists(STF):
    for r in csv.DictReader(io.open(STF, encoding="utf8"), delimiter=TAB):
        STATUS[r["array_id"]] = r


def events(vs):
    # 변이 자리를 정렬해 간격 <= EVGAP 이면 한 사건으로 접는다 → [[자리,...], ...]
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


def blocks(cols, poly, n, minev=3):
    return [c for c in comps(cols, n) if len(events([poly[k] for k in c])) >= minev]


def ev_key(V, n, evlist):
    # 유닛이 사건 j 를 가졌는가 = 그 사건의 변이 중 절반 이상을 가졌는가
    return [tuple(1 if sum(1 for v in e if v in V[i]) * 2 >= len(e) else 0 for e in evlist) for i in range(n)]


def partition(V, n, evlist, minsize=5):
    key = ev_key(V, n, evlist)
    lab = [None] * n
    k = 0
    for p, c in collections.Counter(key).most_common():
        if c < minsize:
            continue
        for i in range(n):
            if key[i] == p:
                lab[i] = k
        k += 1
    return lab


def ari(a, b):
    A = [x if x is not None else ("a", i) for i, x in enumerate(a)]
    B = [x if x is not None else ("b", i) for i, x in enumerate(b)]
    ct = collections.Counter(zip(A, B))
    ra, rb = collections.Counter(A), collections.Counter(B)
    ch2 = lambda x: x * (x - 1) / 2.0
    s = sum(ch2(v) for v in ct.values())
    sa = sum(ch2(v) for v in ra.values())
    sb = sum(ch2(v) for v in rb.values())
    exp = sa * sb / ch2(len(A)) if len(A) > 1 else 0.0
    mx = (sa + sb) / 2.0
    return (s - exp) / (mx - exp) if mx != exp else 1.0


def adj_same(lab, order):
    return sum(1 for a, b in zip(order, order[1:]) if lab[a] is not None and lab[a] == lab[b])


def group_ref(us, tmpdir, tag):
    """그룹 참조 유닛 하나를 만든다.
    자리별 다수결은 유닛 길이가 제각각이면 삽입/결실 뒤쪽을 뭉갠다(실측 : orang_bornean 그룹1
      대표가 실제 유닛 어느 것과도 88% 밖에 안 맞았다). 그래서
        · 최빈 길이 유닛이 그룹의 30% 이상이면 그 무리로 자리별 다수결(기존 방식)
        · 아니면 medoid — 다른 유닛들과 가장 잘 맞는 실제 유닛 하나를 그대로 쓴다
      어느 쪽이든 만든 뒤 그룹 전원을 대표에 정렬해 일치도 중앙값을 재고, 낮으면 표에 적는다."""
    cnt = collections.Counter(len(u) for u in us)
    ml, mc = cnt.most_common(1)[0]
    if mc >= max(5, 0.30 * len(us)):
        pool = [u for u in us if len(u) == ml]
        return ("".join(collections.Counter(u[i] for u in pool).most_common(1)[0][0]
                        for i in range(ml)), "modal_consensus", mc)
    smp = us if len(us) <= 60 else random.sample(us, 60)
    qp = os.path.join(tmpdir, tag + ".q.fa")
    tp = os.path.join(tmpdir, tag + ".t.fa")
    io.open(qp, "w", encoding="utf8", newline=NL).write(
        "".join(">q%d%s%s%s" % (i, NL, u, NL) for i, u in enumerate(smp)))
    io.open(tp, "w", encoding="utf8", newline=NL).write(
        "".join(">t%d%s%s%s" % (i, NL, u, NL) for i, u in enumerate(smp)))
    rr = subprocess.run("minimap2 -cx asm20 -t 4 %s %s 2>/dev/null" % (tp, qp),
                        shell=True, capture_output=True, text=True)
    tot = collections.Counter()
    for L in rr.stdout.splitlines():
        f = L.split(TAB)
        if len(f) < 12 or f[0][1:] == f[5][1:]:
            continue
        for x in f[12:]:
            if x.startswith("AS:i:"):
                tot[int(f[5][1:])] += int(x[5:])
    os.remove(qp)
    os.remove(tp)
    if not tot:
        return (smp[0], "medoid_fallback", 1)
    return (smp[tot.most_common(1)[0][0]], "medoid", len(smp))


def ref_quality(us, ref, tmpdir, tag):
    """그룹 전원을 대표에 정렬해 일치도 중앙값(%) 을 낸다."""
    qp = os.path.join(tmpdir, tag + ".mq.fa")
    tp = os.path.join(tmpdir, tag + ".mt.fa")
    io.open(qp, "w", encoding="utf8", newline=NL).write(
        "".join(">q%d%s%s%s" % (i, NL, u, NL) for i, u in enumerate(us)))
    io.open(tp, "w", encoding="utf8", newline=NL).write(">ref%s%s%s" % (NL, ref, NL))
    rr = subprocess.run("minimap2 -cx asm20 -t 4 %s %s 2>/dev/null" % (tp, qp),
                        shell=True, capture_output=True, text=True)
    best = {}
    for L in rr.stdout.splitlines():
        f = L.split(TAB)
        if len(f) < 11:
            continue
        idt = 100.0 * int(f[9]) / max(1, int(f[10]))
        if f[0] not in best or idt > best[f[0]]:
            best[f[0]] = idt
    v = sorted(best.values())
    return (v[len(v) // 2] if v else 0.0), len(v)


HDR = ["species_key", "array_id", "array_status", "refunit_len", "n_units_OK", "n_var_sites",
       "n_poly_var", "poly_per_kb", "n_poly_event", "n_blocks", "block_event_sizes",
       "null_blocks_mean", "null_blocks_max", "holdout_ARI", "unit_type_verdict", "n_groups",
       "group_sizes", "largest_group_pct", "pos_clustering_ratio", "pos_clustering_p",
       "pos_pattern", "block1_events", "group_unit_files", "group_ref_quality", "note"]
SUM = io.open(os.path.join(OUT, "UNIT_TYPES_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
SUM.write(TAB.join(HDR) + NL)


def w(d):
    SUM.write(TAB.join(str(d.get(k, "-")) for k in HDR) + NL)


for fn in sorted(os.listdir(MAPD)):
    if not fn.endswith(".g8map.tsv"):
        continue
    aid = fn[:-len(".g8map.tsv")]
    stt = STATUS.get(aid, {})
    rl = stt.get("refunit_len", "")
    base = {"species_key": SP, "array_id": aid, "array_status": stt.get("status", "UNLISTED"),
            "refunit_len": rl or "-"}
    if base["array_status"] in ("NON_5S", "NO_UNIT"):
        base["unit_type_verdict"] = "제외(%s)" % base["array_status"]
        base["note"] = stt.get("action", "-")
        w(base)
        continue
    rows = [r for r in csv.DictReader(io.open(os.path.join(MAPD, fn), encoding="utf8"), delimiter=TAB)
            if r["verify"] == "OK"]
    n = len(rows)
    base["n_units_OK"] = n
    if n < 20:
        base["unit_type_verdict"] = "판정불가(OK 유닛<20)"
        w(base)
        continue
    V = [set() if r["varlist"] in ("-", "NOALN", "REV", "") else set(r["varlist"].split(",")) for r in rows]
    cnt = collections.Counter()
    for s in V:
        cnt.update(s)
    poly = sorted([v for v, c in cnt.items() if 0.10 * n <= c <= 0.90 * n], key=pos_of)
    ul = collections.Counter(int(r["unit_len"]) for r in rows).most_common(1)[0][0]
    base.update({"n_var_sites": len(cnt), "n_poly_var": len(poly), "n_poly_event": len(events(poly)),
                 "poly_per_kb": "%.1f" % (1000.0 * len(poly) / ul)})
    if len(events(poly)) < 3:
        base.update({"n_blocks": 0, "unit_type_verdict": "한 종류", "n_groups": 1, "group_sizes": n,
                     "largest_group_pct": "100.0",
                     "note": "다형성 사건 %d 개뿐(변이 %d)" % (len(events(poly)), len(poly))})
        w(base)
        continue
    cols = [set(i for i in range(n) if v in V[i]) for v in poly]
    bl = blocks(cols, poly, n)
    nullc = [len(blocks([set(random.sample(range(n), len(c))) for c in cols], poly, n)) for _ in range(NPERM)]
    nm = sum(nullc) / float(len(nullc))
    base.update({"n_blocks": len(bl), "null_blocks_mean": "%.2f" % nm, "null_blocks_max": max(nullc),
                 "block_event_sizes": ",".join(str(len(events([poly[k] for k in b]))) for b in bl) or "-"})
    if not bl:
        base.update({"unit_type_verdict": "한 종류(연관 블록 없음)", "n_groups": 1, "group_sizes": n,
                     "largest_group_pct": "100.0", "note": "변이는 있으나 떨어진 사건끼리 묶이지 않는다"})
        w(base)
        continue
    bl.sort(key=lambda b: -len(events([poly[k] for k in b])))
    bv = [poly[k] for k in bl[0]]
    ev = events(bv)
    evv = [[v for v in bv if pos_of(v) in set(e)] for e in ev]
    base["block1_events"] = ";".join("%d-%d" % (e[0], e[-1]) if len(e) > 1 else str(e[0]) for e in ev)
    lab = partition(V, n, evv)
    ng = len(set(x for x in lab if x is not None))
    # 홀/짝 교차검증
    od, evn = evv[0::2], evv[1::2]
    base["holdout_ARI"] = "%.3f" % ari(partition(V, n, od), partition(V, n, evn)) if od and evn else "-"
    if ng <= 1:
        base.update({"unit_type_verdict": "한 종류(블록 있으나 패턴 1)", "n_groups": 1, "group_sizes": n,
                     "largest_group_pct": "100.0"})
        w(base)
        continue
    order = sorted(range(n), key=lambda i: int(rows[i]["array_pos"]))
    obs = adj_same(lab, order)
    nulls = []
    for _ in range(NPERM):
        sh = lab[:]
        random.shuffle(sh)
        nulls.append(adj_same(sh, order))
    em = sum(nulls) / float(len(nulls))
    ratio = obs / em if em > 0 else float("nan")
    pv = (sum(1 for x in nulls if x >= obs) + 1.0) / (len(nulls) + 1.0)
    base.update({"pos_clustering_ratio": "%.2f" % ratio, "pos_clustering_p": "%.3f" % pv,
                 "pos_pattern": "구역 분화(이어져 있음)" if (ratio >= 1.5 and pv <= 0.05) else "섞여 있음"})
    mem = collections.defaultdict(list)
    for i, k in enumerate(lab):
        if k is not None:
            mem[k].append(i)
    ks = sorted(mem, key=lambda k: -len(mem[k]))
    regs, meta = [], []
    for nk, k in enumerate(ks):
        for i in mem[k]:
            r = rows[i]
            regs.append("%s:%s-%s" % (r["chrom"], r["chrom_start"], r["chrom_end"]))
            meta.append((nk, r["strand"]))
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
    for (nk, stnd), s in zip(meta, seqs):
        bg[nk].append(rc(s.upper()) if stnd == "-" else s.upper())
    files, sizes, qual = [], [], []
    for nk in sorted(bg):
        us = bg[nk]
        sizes.append(len(us))
        con, how, npool = group_ref(us, os.path.join(OUT, "groups"), "%s.g%d" % (aid, nk + 1))
        md, nq = ref_quality(us, con, os.path.join(OUT, "groups"), "%s.g%d" % (aid, nk + 1))
        qual.append("grp%d:%s(%d):%.2f%%" % (nk + 1, how, npool, md))
        p = os.path.join(OUT, "groups", "%s.grp%d.unit.fa" % (aid, nk + 1))
        io.open(p, "w", encoding="utf8", newline=NL).write(
            ">%s_grp%d n=%d len=%d method=%s member_median_identity=%.2f%s%s%s"
            % (aid, nk + 1, len(us), len(con), how, md, NL, con, NL))
        files.append(os.path.basename(p))
    mf = io.open(os.path.join(OUT, "groups", "%s.members.tsv" % aid), "w", encoding="utf8", newline=NL)
    mf.write(TAB.join(["array_id", "array_pos", "chrom", "chrom_start", "chrom_end", "strand", "group"]) + NL)
    for i, r in enumerate(rows):
        g = lab[i]
        mf.write(TAB.join([aid, r["array_pos"], r["chrom"], r["chrom_start"], r["chrom_end"], r["strand"],
                           ("grp%d" % (ks.index(g) + 1)) if g is not None else "-"]) + NL)
    mf.close()
    base.update({"unit_type_verdict": "여러 종류", "n_groups": len(sizes),
                 "group_sizes": ",".join(str(x) for x in sizes),
                 "largest_group_pct": "%.1f" % (100.0 * sizes[0] / n),
                 "group_unit_files": ",".join(files),
                 "group_ref_quality": " · ".join(qual),
                 "note": "그룹 밖 유닛 %d 개(어느 패턴도 5 개 미만)" % sum(1 for x in lab if x is None)})
    w(base)
    print("%s : 유닛 %d · 사건 %d · 블록 %s · ARI %s -> 그룹 %d %s · 위치응집 %.2f p=%.3f"
          % (aid, n, len(events(poly)), base["block_event_sizes"], base["holdout_ARI"],
             len(sizes), sizes, ratio, pv))
SUM.close()
