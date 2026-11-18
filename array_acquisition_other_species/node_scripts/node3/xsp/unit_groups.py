# -*- coding: utf-8 -*-
# Are the units within an array of one type or several? They are grouped, and a reference unit built per group.
#   Even a species with only one array, as the great apes have, is checked for whether its array is a repeat of one unit, and divided into groups if not.
#   The method uses the per-unit variant lists of the G8 map:
#     - a polymorphic variant is one present in 10-90% of units (a rarer one is noise in an individual unit, and above 90% the reference unit is wrong);
#     - with fewer than 3 polymorphic variants the array counts as a single type;
#     - otherwise the units are clustered on their polymorphic variant vectors, grouping identical patterns and then merging patterns within 15% of each other;
#     - only a cluster with at least 5 units and at least 3 distinguishing variants counts as a group; the rest are other.
#   Per group: the unit sequences are pooled into a group reference unit by a per-position majority, and the within-group identity measured by blastn.
#   Arguments: species key, fasta, map folder, output folder.
#   Outputs: groups/<array_id>.groups.tsv, groups/<array_id>.g<k>.unit.fa and UNIT_GROUPS_<species>.tsv.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import sys
import csv
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
SP, FA, MAPD, OUT = sys.argv[1:5]
os.makedirs(os.path.join(OUT, "groups"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))
SUM = io.open(os.path.join(OUT, "UNIT_GROUPS_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
SUM.write(TAB.join(["species", "array_id", "n_units_OK", "n_var_sites", "n_poly_var", "verdict", "n_groups",
                    "group_sizes", "largest_group_pct", "group_defining_vars", "note"]) + NL)
for fn in sorted(os.listdir(MAPD)):
    if not fn.endswith(".g8map.tsv"):
        continue
    aid = fn[:-len(".g8map.tsv")]
    rows = [r for r in csv.DictReader(io.open(os.path.join(MAPD, fn), encoding="utf8"), delimiter=TAB)]
    ok = [r for r in rows if r["verify"] == "OK"]
    if len(ok) < 5:
        SUM.write(TAB.join([SP, aid, str(len(ok)), "NA", "NA", "판정불가(OK 유닛<5)", "NA", "NA", "NA", "NA", ""]) + NL)
        continue
    V = []
    for r in ok:
        s = set() if r["varlist"] in ("-", "NOALN", "REV", "") else set(r["varlist"].split(","))
        V.append(s)
    cnt = collections.Counter()
    for s in V:
        cnt.update(s)
    n = len(V)
    poly = sorted([v for v, c in cnt.items() if 0.10 * n <= c <= 0.90 * n])
    if len(poly) < 3:
        SUM.write(TAB.join([SP, aid, str(n), str(len(cnt)), str(len(poly)), "한 종류", "1", str(n), "100.0", "-",
                            "다형성 변이 %d 개뿐" % len(poly)]) + NL)
        print("%s : units %d, polymorphic variants %d -> one kind" % (aid, n, len(poly)))
        continue
    idx = dict((v, i) for i, v in enumerate(poly))
    vec = [tuple(1 if v in s else 0 for v in poly) for s in V]
    pat = collections.Counter(vec)
    reps = [p for p, c in pat.most_common()]
    thr = max(1, int(0.15 * len(poly)))
    groups = []          # [(the representative pattern, [unit indices])]
    for p in reps:
        placed = False
        for g in groups:
            if sum(1 for a, b in zip(p, g[0]) if a != b) <= thr:
                g[1].append(p)
                placed = True
                break
        if not placed:
            groups.append([p, [p]])
    gmem = collections.defaultdict(list)
    for i, v in enumerate(vec):
        for gi, g in enumerate(groups):
            if v in g[1]:
                gmem[gi].append(i)
                break
    big = [(gi, m) for gi, m in gmem.items() if len(m) >= 5]
    big.sort(key=lambda t: -len(t[1]))
    if len(big) <= 1:
        SUM.write(TAB.join([SP, aid, str(n), str(len(cnt)), str(len(poly)), "한 종류(다형성 있으나 묶임 1)", "1",
                            str(n), "100.0", "-", "다형성 변이 %d 개" % len(poly)]) + NL)
        print("%s : units %d, polymorphic %d, clusters 1 -> one kind" % (aid, n, len(poly)))
        continue
    # the distinguishing variants and the reference unit, per group
    sizes, defs, paths = [], [], []
    regs, meta = [], []
    for k, (gi, mem) in enumerate(big):
        sizes.append(len(mem))
        inn = collections.Counter()
        for i in mem:
            inn.update(V[i])
        out = collections.Counter()
        for i in range(n):
            if i not in mem:
                out.update(V[i])
        nout = n - len(mem)
        d = [v for v in poly if inn[v] >= 0.8 * len(mem) and out[v] <= 0.2 * nout]
        defs.append(";".join(d[:8]) if d else "-")
        for i in mem:
            r = ok[i]
            regs.append("%s:%s-%s" % (r["chrom"], r["chrom_start"], r["chrom_end"]))
            meta.append((k, r["strand"]))
    rf = os.path.join(OUT, "groups", aid + ".regions.txt")
    io.open(rf, "w", encoding="utf8", newline=NL).write(NL.join(regs) + NL)
    r = subprocess.run(["samtools", "faidx", "-r", rf, FA], capture_output=True, text=True)
    seqs, cur = [], []
    for L in r.stdout.splitlines():
        if L.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(L.strip())
    if cur:
        seqs.append("".join(cur))
    os.remove(rf)
    bygrp = collections.defaultdict(list)
    for (k, st), s in zip(meta, seqs):
        bygrp[k].append(rc(s.upper()) if st == "-" else s.upper())
    for k in sorted(bygrp):
        us = bygrp[k]
        ml = collections.Counter(len(u) for u in us).most_common(1)[0][0]
        pool = [u for u in us if len(u) == ml]
        con = "".join(collections.Counter(u[i] for u in pool).most_common(1)[0][0] for i in range(ml)) if pool else us[0]
        p = os.path.join(OUT, "groups", "%s.g%d.unit.fa" % (aid, k + 1))
        io.open(p, "w", encoding="utf8", newline=NL).write(">%s_g%d n=%d len=%d%s%s%s" % (aid, k + 1, len(us), ml, NL, con, NL))
        paths.append(os.path.basename(p))
    SUM.write(TAB.join([SP, aid, str(n), str(len(cnt)), str(len(poly)), "여러 종류", str(len(big)),
                        ",".join(str(x) for x in sizes), "%.1f" % (100.0 * sizes[0] / n), " | ".join(defs[:3]),
                        "그룹 유닛 파일 " + ",".join(paths)]) + NL)
    print("%s : units %d, polymorphic %d -> groups %d %s" % (aid, n, len(poly), len(big), sizes))
SUM.close()
