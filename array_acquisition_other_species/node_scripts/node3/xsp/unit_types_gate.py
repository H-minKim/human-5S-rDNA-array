# -*- coding: utf-8 -*-
# 유닛 그룹을 확실한 차이가 있을 때만 남긴다 (gate)
# "확실한 차이가 있을 때만 한 array 안에서 unit 별 그룹을 나눠라"
#   통과 조건 (모두 만족)
#     G1 홀/짝 교차검증 ARI >= 0.70        — 변이 절반만으로도 같은 분할이 재현된다
#     G2 귀무 대조 우연 블록 0             — 무작위로는 나오지 않는다
# G3 polymorphism density < 50/kb: the group call is withheld when the reference unit itself is doubtful
#     G4 그룹마다 유닛 >= 10 개 & >= 5%    — 못 미치는 그룹은 '미분류' 로 떨구고 남은 것만 센다
#     G5 남은 그룹 >= 2 개이고, 그룹 대표끼리 실제 차이(blastn 불일치+갭) >= 10 자리
#     G6 남은 그룹의 대표가 제 구성원을 대표한다 — 구성원 일치도 중앙값 >= 98%
#   떨어진 그룹의 참조 유닛 파일은 dropped/ 로 옮긴다(지우지 않는다)
#   인자 : 00_catalog 폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import re
import glob
import shutil
import subprocess
import itertools

TAB = chr(9)
NL = chr(10)
BASE = sys.argv[1]
UT = os.path.join(BASE, "unit_types")
GD = os.path.join(UT, "groups")
DROP = os.path.join(UT, "groups_dropped")
os.makedirs(DROP, exist_ok=True)
rows = list(csv.DictReader(io.open(os.path.join(BASE, "UNIT_TYPE_CALL.tsv"), encoding="utf8"), delimiter=TAB))


def seqof(p):
    L = io.open(p, encoding="utf8").read().split(NL)
    return "".join(x.strip() for x in L[1:])


def hdr(p):
    return io.open(p, encoding="utf8").read().split(NL)[0]


def member_ident(p):
    h = hdr(p)
    m = re.search(r"member_median_identity=([0-9.]+)", h)
    return float(m.group(1)) if m else 0.0


def real_diff(pa, pb):
    """두 그룹 대표가 실제로 몇 자리 다른가 — blastn 최고 HSP 의 불일치 + 갭 자리"""
    r = subprocess.run("blastn -task blastn -word_size 11 -dust no -query %s -subject %s "
                       "-outfmt '6 mismatch gaps bitscore' 2>/dev/null" % (pa, pb),
                       shell=True, capture_output=True, text=True)
    best = None
    for L in r.stdout.splitlines():
        mm, gp, bs = L.split(TAB)
        if best is None or float(bs) > best[0]:
            best = (float(bs), int(mm) + int(gp))
    return best[1] if best else 10 ** 6


HDR = list(rows[0].keys()) + ["gate_result", "gate_fail", "final_n_groups", "final_group_sizes",
                              "final_group_files", "dropped_groups", "group_pair_real_diff",
                              "group_ref_member_identity_min"]
out = io.open(os.path.join(BASE, "UNIT_TYPE_FINAL.tsv"), "w", encoding="utf8", newline=NL)
out.write(TAB.join(HDR) + NL)
npass = 0
for r in sorted(rows, key=lambda x: (x["species_key"], x["array_id"])):
    d = dict(r)
    d.update({"gate_result": "-", "gate_fail": "-", "final_n_groups": "1",
              "final_group_sizes": r.get("n_units_OK", "-"), "final_group_files": "-", "dropped_groups": "-"})
    if not r["unit_type_verdict"].startswith(chr(9733)):
        d["gate_result"] = "해당없음(여러 종류 아님)"
        out.write(TAB.join(str(d.get(k, "-")) for k in HDR) + NL)
        continue
    n = int(r["n_units_OK"])
    fails = []
    try:
        if float(r["holdout_ARI"]) < 0.70:
            fails.append("G1 ARI %s < 0.70" % r["holdout_ARI"])
    except ValueError:
        fails.append("G1 ARI 없음")
    if r["null_blocks_max"] not in ("0", "-"):
        fails.append("G2 귀무 블록 %s" % r["null_blocks_max"])
    try:
        if float(r["poly_per_kb"]) >= 50:
            fails.append("G3 다형성 밀도 %s/kb" % r["poly_per_kb"])
    except ValueError:
        pass
    sizes = [int(x) for x in r["group_sizes"].split(",")]
    fs = sorted(glob.glob(os.path.join(GD, r["array_id"] + ".grp*.unit.fa")),
                key=lambda p: int(p.split(".grp")[1].split(".")[0]))
    keep, drop = [], []
    for i, (sz, f) in enumerate(zip(sizes, fs)):
        (keep if (sz >= 10 and sz >= 0.05 * n) else drop).append((i, sz, f))
    if len(keep) < 2:
        fails.append("G4 유닛 >=10 & >=5%% 인 그룹이 %d 개" % len(keep))
    else:
        md = min(real_diff(a[2], b[2]) for a, b in itertools.combinations(keep, 2))
        if md < 10:
            fails.append("G5 그룹 대표 실제 차이 %d 자리 < 10" % md)
        lo = min(member_ident(f) for _i, _s, f in keep)
        if lo < 98.0:
            fails.append("G6 그룹 대표 구성원 일치도 %.2f%% < 98%%" % lo)
    if fails:
        d.update({"gate_result": "불통과 → 한 종류로 다룬다", "gate_fail": " · ".join(fails),
                  "dropped_groups": ",".join(os.path.basename(f) for _i, _s, f in keep + drop)})
        for _i, _s, f in keep + drop:
            if os.path.exists(f):
                shutil.move(f, os.path.join(DROP, os.path.basename(f)))
    else:
        npass += 1
        for _i, _s, f in drop:
            if os.path.exists(f):
                shutil.move(f, os.path.join(DROP, os.path.basename(f)))
        d.update({"gate_result": "통과 → 그룹별 참조 유닛 사용",
                  "group_pair_real_diff": md, "group_ref_member_identity_min": "%.2f" % lo,
                  "final_n_groups": len(keep),
                  "final_group_sizes": ",".join(str(s) for _i, s, _f in keep),
                  "final_group_files": ",".join(os.path.basename(f) for _i, _s, f in keep),
                  "dropped_groups": ",".join(os.path.basename(f) for _i, _s, f in drop) or "-"})
    out.write(TAB.join(str(d.get(k, "-")) for k in HDR) + NL)
out.close()
print("gate 통과 array %d" % npass)
for r in csv.DictReader(io.open(os.path.join(BASE, "UNIT_TYPE_FINAL.tsv"), encoding="utf8"), delimiter=TAB):
    if r["gate_result"].startswith("해당없음"):
        continue
    print("  %-46s %-28s %s" % (r["array_id"][:46], r["gate_result"],
                                r["gate_fail"] if r["gate_fail"] != "-" else
                                "그룹 %s (%s)" % (r["final_n_groups"], r["final_group_sizes"])))
