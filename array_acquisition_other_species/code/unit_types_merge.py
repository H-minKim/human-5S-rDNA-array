# -*- coding: utf-8 -*-
# 종별 UNIT_TYPES_*.tsv 를 한 표로 합치고, 그룹을 가르는 사건이 유전자부인지 스페이서부인지 적는다
#   G8 틀에서 유닛 = [스페이서][유전자] 이고 유전자가 끝에 온다 → 자리 > 유닛길이-120 이면 유전자부
#   그룹 대표 서열끼리의 길이차도 함께 적는다
#   인자 : 00_catalog 폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import glob
import itertools

TAB = chr(9)
NL = chr(10)
BASE = sys.argv[1]
UT = os.path.join(BASE, "unit_types")
GENE = 120
rows = []
for f in sorted(glob.glob(os.path.join(UT, "UNIT_TYPES_*.tsv"))):
    rows += list(csv.DictReader(io.open(f, encoding="utf8"), delimiter=TAB))


def seqof(p):
    L = io.open(p, encoding="utf8").read().split(NL)
    return "".join(x.strip() for x in L[1:])


extra = ["block1_region", "group_len_range", "group_pair_min_diff"]
HDR = list(rows[0].keys()) + extra
out = io.open(os.path.join(BASE, "UNIT_TYPE_CALL.tsv"), "w", encoding="utf8", newline=NL)
out.write(TAB.join(HDR) + NL)
for r in sorted(rows, key=lambda x: (x["species_key"], x["array_id"])):
    d = dict(r)
    d.update(dict((k, "-") for k in extra))
    try:
        ul = int(r["refunit_len"])
    except ValueError:
        ul = 0
    if r["block1_events"] not in ("-", ""):
        ing = out_g = 0
        for e in r["block1_events"].split(";"):
            p = int(e.split("-")[0])
            if ul and p > ul - GENE:
                ing += 1
            else:
                out_g += 1
        d["block1_region"] = "스페이서 %d · 유전자 %d" % (out_g, ing)
    fs = sorted(glob.glob(os.path.join(UT, "groups", r["array_id"] + ".grp*.unit.fa")))
    if len(fs) >= 2:
        ss = [seqof(p) for p in fs]
        ln = [len(s) for s in ss]
        d["group_len_range"] = "%d-%d" % (min(ln), max(ln))
        best = None
        for a, b in itertools.combinations(ss, 2):
            v = sum(1 for x, y in zip(a, b) if x != y) if len(a) == len(b) else abs(len(a) - len(b))
            best = v if best is None else min(best, v)
        d["group_pair_min_diff"] = best
    out.write(TAB.join(str(d.get(k, "-")) for k in HDR) + NL)
out.close()
import collections
c = collections.Counter(r["unit_type_verdict"].split("(")[0] for r in rows)
for k, v in c.most_common():
    print("%-24s %3d" % (k, v))
print("%-24s %3d" % ("합계", len(rows)))
