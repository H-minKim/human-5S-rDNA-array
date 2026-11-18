# -*- coding: utf-8 -*-
# Are the units within an array of one type or several? Version 2, based on association, with a null control and positional clustering reported alongside.
#   Splitting whenever the variant patterns merely differ is not used: it cuts a human array into 6 groups, since random variants also give different patterns.
#   The criterion: do the variants travel together as a block?
#     - In a single-type array the variants are scattered independently across units, with no association.
#     - The presence of subtypes means several variants appear together in the same units, a strong association.
#   method
#     1) A polymorphic variant is one present in 10-90% of units.
#     2) the association phi between variant pairs; |phi| >= 0.8 links them
#     3) a connected component of at least 3 counts as a block, and the units are divided by the pattern over the largest block
#     4) null control: the units are reassigned at random at the same frequency per variant and the whole procedure repeated NPERM times, recording the number of blocks arising by chance
#     5) Positional clustering: are the groups contiguous along the array (zonation) or interleaved (unit types)?
#        plus the observed over expected ratio of neighbouring unit pairs in the same group, with a permutation p value
#   No call is made with fewer than 20 units, the frequency estimates being too coarse.
#   Arguments: species key, fasta, map folder, output folder, [shuffles = 50].
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import csv
import random
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
SP, FA, MAPD, OUT = sys.argv[1:5]
NPERM = int(sys.argv[5]) if len(sys.argv) > 5 else 50
random.seed(20260913)
os.makedirs(os.path.join(OUT, "groups"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))


def blocks(cols, n, thr=0.8, minsize=3):
    # cols: per variant, the set of unit indices carrying it. They are linked at an association |phi| >= thr and the components built.
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
        if len(comp) >= minsize:
            out.append(sorted(comp))
    return out


def adjacency_same(lab, order):
    # the number of times two neighbouring units in array order fall in the same group
    s = 0
    for a, b in zip(order, order[1:]):
        if lab[a] is not None and lab[a] == lab[b]:
            s += 1
    return s


HDR = ["species", "array_id", "n_units_OK", "n_var_sites", "n_poly_var", "n_blocks", "block_sizes",
       "null_blocks_mean", "null_blocks_max", "verdict", "n_groups", "group_sizes", "largest_group_pct",
       "pos_clustering_ratio", "pos_clustering_p", "pos_pattern", "block1_vars", "note"]
SUM = io.open(os.path.join(OUT, "UNIT_GROUPS2_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
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
        base.update({"verdict": "판정불가(OK 유닛<20)"})
        w(base)
        continue
    V = [set() if r["varlist"] in ("-", "NOALN", "REV", "") else set(r["varlist"].split(",")) for r in ok]
    cnt = collections.Counter()
    for s in V:
        cnt.update(s)
    poly = sorted([v for v, c in cnt.items() if 0.10 * n <= c <= 0.90 * n])
    base.update({"n_var_sites": len(cnt), "n_poly_var": len(poly)})
    if len(poly) < 3:
        base.update({"n_blocks": 0, "verdict": "한 종류", "n_groups": 1, "group_sizes": n,
                     "largest_group_pct": "100.0", "note": "다형성 변이 %d 개뿐" % len(poly)})
        w(base)
        print("%s : units %d, polymorphic %d -> one kind" % (aid, n, len(poly)))
        continue
    cols = [set(i for i in range(n) if v in V[i]) for v in poly]
    bl = blocks(cols, n)
    nullc = []
    allidx = list(range(n))
    for _ in range(NPERM):
        nullc.append(len(blocks([set(random.sample(allidx, len(c))) for c in cols], n)))
    nm = sum(nullc) / float(len(nullc))
    base.update({"n_blocks": len(bl), "block_sizes": ",".join(str(len(b)) for b in bl) if bl else "-",
                 "null_blocks_mean": "%.2f" % nm, "null_blocks_max": max(nullc)})
    if not bl:
        base.update({"verdict": "한 종류(연관 블록 없음)", "n_groups": 1, "group_sizes": n,
                     "largest_group_pct": "100.0", "note": "변이는 있으나 서로 묶이지 않는다"})
        w(base)
        print("%s : units %d, polymorphic %d, blocks 0 (null %.1f) -> one kind" % (aid, n, len(poly), nm))
        continue
    bl.sort(key=lambda b: -len(b))
    b1 = bl[0]
    base["block1_vars"] = ";".join(poly[k] for k in b1[:8])
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
    # the reference unit per group
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
    print("%s : units %d, polymorphic %d, blocks %s (null %.1f) -> groups %d %s, positional clustering %.2f p=%.3f"
          % (aid, n, len(poly), [len(b) for b in bl], nm, len(grp), sizes, ratio, pv))
SUM.close()
