# -*- coding: utf-8 -*-
# 제브라피시 5S array 전수 표 — 염색체 · 크기 · 유닛 · 메틸 · flank 대비
#   제브라피시는 5S 가 두 갈래로 알려져 있다
#     chr4  = 모계형(난자 전용 · 체세포에서 침묵)
#     chr18 = 체세포형(배아 발생 중 켜져 체세포에서 발현)
#   우리 자료는 성체 뇌(체세포) 6 개체다. array 마다 그 둘 중 어디에 해당하는지,
#   그리고 array 가 제 주변(flank 100 kb)보다 높은지 낮은지를 한 표에 놓는다.
#   인자 : data_v2 폴더 출력.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import glob
import collections

TAB = chr(9)
NL = chr(10)
B, OUT = sys.argv[1], sys.argv[2]
CH = {}
for i, acc in enumerate(["CM1111%02d.1" % (79 + k) for k in range(25)]):
    CH[acc] = "chr%d" % (i + 1)
# 실제 접근번호는 CM111179.1 부터 연속 — 위 생성식과 맞춘다
CH = dict(("CM%d.1" % (111179 + k), "chr%d" % (k + 1)) for k in range(25))
# array 요약(메틸)
MET = collections.defaultdict(list)
for r in csv.DictReader(io.open(os.path.join(B, "03_reads", "zebrafish_GRCz12tu",
                                             "5S_array_summary.tsv"), encoding="utf8"), delimiter=TAB):
    if r["context"] == "CG":
        MET[r["array_id"]].append((r["sample"], float(r["frac_5mC"]), int(r["calls"]), int(r["N"])))
# flank
AVF = collections.defaultdict(lambda: collections.defaultdict(dict))
p = os.path.join(B, "09_array_vs_flank", "AVF.tsv")
if os.path.exists(p):
    for r in csv.DictReader(io.open(p, encoding="utf8"), delimiter=TAB):
        if r["species_key"] != "zebrafish_GRCz12tu" or r["frac_5mC"] == "-":
            continue
        AVF[r["array_id"]][r["sample"]][r["region"]] = (float(r["frac_5mC"]), int(r["n_c"]))
# 참조유닛 제원
PR = {}
for r in csv.DictReader(io.open(os.path.join(B, "00_catalog", "REFUNIT_PARAMS_MERGED.tsv"),
                                encoding="utf8"), delimiter=TAB):
    if r["species"] == "zebrafish_GRCz12tu":
        PR[r["array_id"]] = r
# 유닛 종류 판정
UT = {}
p = os.path.join(B, "00_catalog", "UNIT_REF_TABLE.tsv")
if os.path.exists(p):
    for r in csv.DictReader(io.open(p, encoding="utf8"), delimiter=TAB):
        if r["species_key"] == "zebrafish_GRCz12tu":
            UT[r["array_id"]] = r
HDR = ["array_id", "chrom", "chrom_name", "array_start", "array_end", "span_bp", "refunit_len",
       "g8_units", "unit_type_verdict", "n_samples", "meth_mean", "meth_min", "meth_max",
       "cpg_calls_total", "flank_mean", "array_minus_flank", "locus_type_lit"]
o = io.open(OUT, "w", encoding="utf8", newline=NL)
o.write(TAB.join(HDR) + NL)
rows = []
for aid in sorted(set(list(MET) + list(PR))):
    pr = PR.get(aid, {})
    acc = aid.split("__")[1].split("_")[0]
    ch = CH.get(acc, acc)
    m = MET.get(aid, [])
    vals = [x[1] for x in m]
    calls = sum(x[2] for x in m)
    n = m[0][3] if m else (pr.get("n_units_in_range") or "-")
    fl = []
    for sam, d in AVF.get(aid, {}).items():
        v = [d[k][0] for k in ("flank_left_100kb", "flank_right_100kb")
             if k in d and d[k][1] >= 300]
        if v:
            fl.append(sum(v) / len(v))
    lt = ("chr4 모계형(난자·체세포 침묵)" if ch == "chr4" else
          ("chr18 체세포형(발현)" if ch == "chr18" else "문헌 분류 없음"))
    row = [aid, acc, ch, pr.get("array_start", "-"), pr.get("array_end", "-"),
           (int(pr["array_end"]) - int(pr["array_start"]) + 1) if pr.get("array_end") else "-",
           pr.get("refunit_len", "-"), n,
           UT.get(aid, {}).get("unit_type_verdict", "-"), len(m),
           "%.3f" % (sum(vals) / len(vals)) if vals else "-",
           "%.3f" % min(vals) if vals else "-", "%.3f" % max(vals) if vals else "-",
           calls, "%.3f" % (sum(fl) / len(fl)) if fl else "-",
           "%+.3f" % (sum(vals) / len(vals) - sum(fl) / len(fl)) if (vals and fl) else "-", lt]
    o.write(TAB.join(str(x) for x in row) + NL)
    rows.append(row)
o.close()
rows.sort(key=lambda r: (r[2] != "chr18", -(int(r[7]) if str(r[7]).isdigit() else 0)))
print("%-30s %-7s %8s %7s %6s %8s %8s %9s  %s"
      % ("array", "염색체", "span", "유닛길이", "유닛", "메틸", "flank", "array-fl", "문헌"))
for r in rows:
    print("%-30s %-7s %8s %7s %6s %8s %8s %9s  %s"
          % (r[0].split("__")[1][:30], r[2], r[5], r[6], r[7], r[10], r[14], r[15], r[16]))
