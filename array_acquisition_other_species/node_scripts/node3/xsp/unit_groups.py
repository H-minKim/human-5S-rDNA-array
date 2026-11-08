# -*- coding: utf-8 -*-
# array 안 유닛이 한 종류인가 여러 종류인가 — 그룹핑하고 그룹마다 참조 유닛을 만든다
#   사용자 요구 : array 가 하나뿐인 종(유인원)도 "array 안이 같은 유닛의 반복인지" 확인하고, 아니면 그룹을 나눈다
#   방법 : G8 지도의 유닛별 변이목록을 쓴다
#     · 다형성 변이 = 유닛의 10~90% 에 있는 변이 (드문 변이는 개별 유닛 잡음, 90% 이상은 참조 유닛이 틀린 것)
#     · 다형성 변이가 3 개 미만이면 한 종류로 본다
#     · 아니면 유닛을 다형성 변이 벡터로 묶는다 (같은 패턴끼리 → 거리 15% 이내 패턴 병합)
#     · 유닛 5 개 이상 + 구분 변이 3 개 이상인 묶음만 '그룹'으로 인정, 나머지는 기타
#   그룹마다 : 유닛 서열을 모아 자리별 다수결로 그룹 참조 유닛 생성 + 그룹 안 일치도(blastn)
#   인자 : 종키 fasta 지도폴더 출력폴더
#   산출 : groups/<array_id>.groups.tsv · groups/<array_id>.g<k>.unit.fa · UNIT_GROUPS_<종>.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
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
        print("%s : 유닛 %d · 다형성 변이 %d -> 한 종류" % (aid, n, len(poly)))
        continue
    idx = dict((v, i) for i, v in enumerate(poly))
    vec = [tuple(1 if v in s else 0 for v in poly) for s in V]
    pat = collections.Counter(vec)
    reps = [p for p, c in pat.most_common()]
    thr = max(1, int(0.15 * len(poly)))
    groups = []          # [(대표패턴, [유닛 인덱스])]
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
        print("%s : 유닛 %d · 다형성 %d · 묶음 1 -> 한 종류" % (aid, n, len(poly)))
        continue
    # 그룹마다 구분 변이와 참조 유닛
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
    print("%s : 유닛 %d · 다형성 %d -> 그룹 %d 개 %s" % (aid, n, len(poly), len(big), sizes))
SUM.close()
