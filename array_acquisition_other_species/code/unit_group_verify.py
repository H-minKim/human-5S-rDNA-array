# -*- coding: utf-8 -*-
# 2 단계 관문 — 유닛 하나하나에서 두 그룹이 실제로 갈리는가
#
#   1 단계 관문(unit_types_gate.py)은 array 통계로 판정했다. 그런데 array 통계가 통과해도
#   유닛 단위로 보면 두 그룹이 안 갈릴 수 있다(사용자 지적). 그래서 유닛마다
#   서로 독립인 두 증거를 매기고, 둘이 얼마나 맞는지로 다시 판정한다.
#     증거 A 변이 사건 패턴 (ev_group) · 증거 B 그룹 참조 정렬 점수 (aln_group)
#
#   G7 두 증거 일치율 >= 95%        — 유닛마다 같은 답이 나오는가
#   G8 불일치 유닛 <= 5%
#   G9 0 근처 빈칸 >= 0.10        — |margin_frac| 의 최소값. 진짜 두 종류면 점수차가 이봉이라
#                                    0 근처가 비어 있다. 연속체면 0 을 가로지르는 유닛이 있다.
#   G10 미정 유닛 <= 5%             — 어느 그룹으로도 못 간 유닛. 두 종류가 진짜면 거의 없어야 한다
#                                    (실측 고릴라 0.0% · 수마트라오랑 0.7% vs 침팬지 10.7%)
#
#   떨어지면 그 array 는 한 종류로 되돌리고(그룹 참조를 groups_dropped/ 로) array 참조만 쓴다.
#   인자 : 00_catalog 폴더 유닛판정표폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import glob
import shutil
import collections

TAB = chr(9)
NL = chr(10)
B, UGC = sys.argv[1], sys.argv[2]
UT = os.path.join(B, "unit_types")
GD = os.path.join(UT, "groups")
DROP = os.path.join(UT, "groups_dropped")
os.makedirs(DROP, exist_ok=True)
rows = list(csv.DictReader(io.open(os.path.join(B, "UNIT_TYPE_FINAL.tsv"), encoding="utf8"), delimiter=TAB))
stat = {}
for f in sorted(glob.glob(os.path.join(UGC, "UNIT_GROUP_CALL_*.tsv"))):
    rr = list(csv.DictReader(io.open(f, encoding="utf8"), delimiter=TAB))
    aid = rr[0]["array_id"]
    n = len(rr)
    fin = collections.Counter(r["final_group"] for r in rr)
    both = [r for r in rr if r["ev_group"] in ("grp1", "grp2") and r["aln_group"] in ("grp1", "grp2")]
    ag = sum(1 for r in both if r["ev_group"] == r["aln_group"])
    mf = sorted(abs(float(r["margin_frac"])) for r in rr if r["margin_frac"] not in ("-", ""))
    stat[aid] = dict(n=n, n_both=len(both), agree=ag,
                     agree_pct=(100.0 * ag / len(both)) if both else 0.0,
                     n_conf=sum(v for k, v in fin.items() if k.startswith("grp")),
                     n_undecided=fin["미정"], n_conflict=fin["불일치"],
                     gap=mf[0] if mf else 0.0)
HDR = list(rows[0].keys()) + ["unit_level_agree_pct", "unit_level_conflict", "unit_level_undecided",
                              "margin_gap_at_zero", "stage2_result", "stage2_fail", "final_unit_types"]
out = io.open(os.path.join(B, "UNIT_TYPE_FINAL2.tsv"), "w", encoding="utf8", newline=NL)
out.write(TAB.join(HDR) + NL)
keepn = 0
show = []
for r in sorted(rows, key=lambda x: (x["species_key"], x["array_id"])):
    d = dict(r)
    d.update({"unit_level_agree_pct": "-", "unit_level_conflict": "-", "unit_level_undecided": "-",
              "margin_gap_at_zero": "-", "stage2_result": "해당없음", "stage2_fail": "-",
              "final_unit_types": "1"})
    s = stat.get(r["array_id"])
    if r["gate_result"].startswith("통과") and s:
        fails = []
        if s["agree_pct"] < 95.0:
            fails.append("G7 두 증거 일치율 %.1f%% < 95%%" % s["agree_pct"])
        if s["n_conflict"] > 0.05 * s["n"]:
            fails.append("G8 불일치 %d/%d (%.1f%%) > 5%%" % (s["n_conflict"], s["n"],
                                                          100.0 * s["n_conflict"] / s["n"]))
        if s["gap"] < 0.10:
            fails.append("G9 0 근처 빈칸 %.3f < 0.10 (연속체)" % s["gap"])
        if s["n_undecided"] > 0.05 * s["n"]:
            fails.append("G10 미정 %d/%d (%.1f%%) > 5%%" % (s["n_undecided"], s["n"],
                                                          100.0 * s["n_undecided"] / s["n"]))
        d.update({"unit_level_agree_pct": "%.1f" % s["agree_pct"], "unit_level_conflict": s["n_conflict"],
                  "unit_level_undecided": s["n_undecided"], "margin_gap_at_zero": "%.3f" % s["gap"]})
        if fails:
            d.update({"stage2_result": "불통과 → 한 종류로 되돌림", "stage2_fail": " · ".join(fails),
                      "final_unit_types": "1"})
            for f in r["final_group_files"].split(","):
                p = os.path.join(GD, f)
                if os.path.exists(p):
                    shutil.move(p, os.path.join(DROP, f))
        else:
            keepn += 1
            d.update({"stage2_result": "통과 → 유닛마다 그룹 참조 사용",
                      "final_unit_types": r["final_n_groups"]})
        show.append((r, d, s))
    out.write(TAB.join(str(d.get(k, "-")) for k in HDR) + NL)
out.close()
print("2 단계 통과 array %d / %d" % (keepn, len(show)))
print("%-24s %5s %8s %8s %8s %8s  %s" % ("종", "유닛", "일치율", "불일치", "미정", "0빈칸", "결과"))
for r, d, s in show:
    print("%-24s %5d %7.1f%% %8d %8d %8.3f  %s %s"
          % (r["species_key"][:24], s["n"], s["agree_pct"], s["n_conflict"], s["n_undecided"], s["gap"],
             d["stage2_result"], d["stage2_fail"] if d["stage2_fail"] != "-" else ""))
