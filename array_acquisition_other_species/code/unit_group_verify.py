# -*- coding: utf-8 -*-
# The second-stage gate: do the two groups actually separate unit by unit?
#
#   The first-stage gate (unit_types_gate.py) judges from array-level statistics. Those can pass
#   while the two groups still fail to separate unit by unit. So each unit is given
#   two independent lines of evidence, and the call is remade from how well the two agree.
#     Evidence A is the variant event pattern (ev_group); evidence B is the alignment score against the group references (aln_group).
#
#   G7 the two lines of evidence agree on at least 95% of units
#   G8 at most 5% of units disagree
#   G9 the gap around 0 is at least 0.10, as the minimum of |margin_frac|. With two real types the score difference is bimodal
#                                    and empty near 0; in a continuum some units straddle 0.
#   G10 at most 5% of units are undecided, that is assigned to neither group; with two real types there should be almost none
#                                    (measured: gorilla 0.0%, Sumatran orangutan 0.7%, against chimpanzee 10.7%).
#
#   On failure the array reverts to a single type: the group references are moved to groups_dropped/ and only the array reference is used.
#   Arguments: the 00_catalog folder, the folder of per-unit call tables.
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
