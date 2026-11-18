# -*- coding: utf-8 -*-
# An independent check on the call that several types are present: the block variants are split as events into an odd and an even half,
#   the units are partitioned separately from each half, and the two partitions are compared.
#   - If the subtypes (haplotypes) are real, half the variants still produce the same partition.
#   - If it is noise, the two partitions disagree.
#   Agreement is the Adjusted Rand Index (0 = chance, 1 = identical) plus the fraction of unit pairs that agree.
#   Also reported: how many positions the group representative sequences actually differ at, and the variant frequency within and outside each group.
#   Arguments: the map folder, [EVGAP = 20].
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import re
import sys
import csv
import collections

TAB = chr(9)
MAPD = sys.argv[1]
EVGAP = int(sys.argv[2]) if len(sys.argv) > 2 else 20
NUM = re.compile(r"^([0-9]+)")
pos_of = lambda v: int(NUM.match(v).group(1)) if NUM.match(v) else -1


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


def ari(a, b):
    # Adjusted Rand Index; a None label is treated as its own unique label.
    A = [x if x is not None else ("a", i) for i, x in enumerate(a)]
    B = [x if x is not None else ("b", i) for i, x in enumerate(b)]
    ct = collections.Counter(zip(A, B))
    ra = collections.Counter(A)
    rb = collections.Counter(B)
    ch2 = lambda x: x * (x - 1) / 2.0
    s = sum(ch2(v) for v in ct.values())
    sa = sum(ch2(v) for v in ra.values())
    sb = sum(ch2(v) for v in rb.values())
    n = len(A)
    exp = sa * sb / ch2(n) if n > 1 else 0.0
    mx = (sa + sb) / 2.0
    return (s - exp) / (mx - exp) if mx != exp else 1.0


def partition(V, n, sel):
    key = [tuple(1 if v in V[i] else 0 for v in sel) for i in range(n)]
    pat = collections.Counter(key)
    lab = [None] * n
    for k, (p, c) in enumerate(pat.most_common()):
        if c < 5:
            continue
        for i in range(n):
            if key[i] == p:
                lab[i] = k
    return lab


print("%-46s %5s %6s %6s %7s %7s %s" % ("array", "유닛", "블록변이", "사건", "ARI", "pair match", "홀/짝 그룹크기"))
for fn in sorted(os.listdir(MAPD)):
    if not fn.endswith(".g8map.tsv"):
        continue
    aid = fn[:-len(".g8map.tsv")]
    rows = [r for r in csv.DictReader(io.open(os.path.join(MAPD, fn), encoding="utf8"), delimiter=TAB)
            if r["verify"] == "OK"]
    n = len(rows)
    if n < 20:
        continue
    V = [set() if r["varlist"] in ("-", "NOALN", "REV", "") else set(r["varlist"].split(",")) for r in rows]
    cnt = collections.Counter()
    for s in V:
        cnt.update(s)
    poly = sorted([v for v, c in cnt.items() if 0.10 * n <= c <= 0.90 * n], key=pos_of)
    if len(events(poly)) < 3:
        continue
    cols = [set(i for i in range(n) if v in V[i]) for v in poly]
    cs = [c for c in comps(cols, n) if len(events([poly[k] for k in c])) >= 3]
    if not cs:
        continue
    cs.sort(key=lambda c: -len(events([poly[k] for k in c])))
    b1 = cs[0]
    bv = [poly[k] for k in b1]
    ev = events(bv)
    byev = collections.defaultdict(list)
    for v in bv:
        for j, e in enumerate(ev):
            if pos_of(v) in e:
                byev[j].append(v)
                break
    odd = [v for j in byev if j % 2 == 0 for v in byev[j]]
    even = [v for j in byev if j % 2 == 1 for v in byev[j]]
    if len(odd) < 1 or len(even) < 1:
        continue
    la, lb = partition(V, n, odd), partition(V, n, even)
    same = tot = 0
    for i in range(n):
        for j in range(i + 1, n):
            tot += 1
            if (la[i] == la[j] and la[i] is not None) == (lb[i] == lb[j] and lb[i] is not None):
                same += 1
    sa = sorted(collections.Counter(x for x in la if x is not None).values(), reverse=True)
    sb = sorted(collections.Counter(x for x in lb if x is not None).values(), reverse=True)
    print("%-46s %5d %6d %6d %7.3f %7.3f %s / %s"
          % (aid[:46], n, len(bv), len(ev), ari(la, lb), same / float(tot), sa, sb))
