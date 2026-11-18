# -*- coding: utf-8 -*-
# Are the units within an array of one type or several (unit_types).
#
# what is being decided
#   Each array gets a reference unit, but only after confirming that the array really is a repeat of one unit;
#   otherwise this is reported and the units are divided into groups, each with its own reference unit.
#
# The axis of the decision: do variant events that lie far apart travel together as a block?
#   - One type: the variants are scattered independently across units, with no association.
#   - Several types: distant events appear together in the same units, a strong association, that is a haplotype.
#
# Three rules, and what each avoids:
#   The criterion is association (phi): splitting whenever the variant patterns merely differ cuts a human array into 6 groups.
#   Adjacent positions are collapsed into one event: counting variants clustered in a narrow interval (27 within 1397-1527)
#          as 27 variants promotes a single event to a block.
#   The partition is keyed on events: keying it on individual variants fragments one group
#          (Sumatran orangutan groups 2, 3 and 5 differ by 0.3-0.5%).
#
# procedure
#   0) Arrays marked NON_5S or NO_UNIT in the status table are excluded.
#   1) A polymorphic variant is one present in 10-90% of units; adjacent positions (gap <= EVGAP) collapse into one event.
#   2) Variant pairs with |phi| >= 0.8 are linked into connected components, and only components of at least 3 events count as a block.
#   3) Units are partitioned by their event pattern over the largest block; only patterns with at least 5 units form a group.
#   4) Null control: units are reassigned at random at the same frequency per variant and the whole procedure repeated NPERM times.
#   5) Odd/even cross-validation: the block events are split into an odd and an even half, each partitioning the units, and reproducibility measured by ARI.
#   6) Positional clustering: are the groups contiguous along the array (zonation) or interleaved (unit types)?
#   7) A group reference unit is built per group by a per-position majority.
#   No call is made with fewer than 20 OK units.
#
# Arguments: species key, fasta, map folder, output folder, status table, [NPERM = 50], [EVGAP = 20].
# MIT License. Copyright (c) 2026 Hyungmin Kim
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
    # The variant positions are sorted and collapsed into one event where the gap is at most EVGAP, giving [[position, ...], ...].
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
    # A unit carries event j if it carries at least half of that event's variants.
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
    """Build one reference unit for a group.
    Per-position majority vote blurs everything after an insertion or deletion when unit lengths vary (measured : the representative
      of orang_bornean group 1 matched no real unit better than 88%). So
        - if units of the modal length are 30% or more of the group, per-position majority vote over them
        - otherwise the medoid: the one real unit that matches the others best is used as it is
      Either way every member of the group is then aligned to the representative, the median identity is measured, and a low value is written in the table."""
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
    """Align every member of the group to the representative and give the median identity (%)."""
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
    os.remove(qp)
    os.remove(tp)
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
    # odd/even cross-validation
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
    print("%s : units %d, events %d, blocks %s, ARI %s -> groups %d %s, positional clustering %.2f p=%.3f"
          % (aid, n, len(events(poly)), base["block_event_sizes"], base["holdout_ARI"],
             len(sizes), sizes, ratio, pv))
SUM.close()
