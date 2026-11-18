# -*- coding: utf-8 -*-
# Per array: the number of unit types, and for each type the length of its reference unit, with the period and the unit count.
#   With one type there is a single array reference unit; with several there are group reference units.
#   Argument: the 00_catalog folder.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import io
import os
import re
import sys
import csv

TAB = chr(9)
NL = chr(10)
B = sys.argv[1]
GD = os.path.join(B, "unit_types", "groups")
FN = "UNIT_TYPE_FINAL2.tsv" if os.path.exists(os.path.join(B, "UNIT_TYPE_FINAL2.tsv")) else "UNIT_TYPE_FINAL.tsv"
rows = list(csv.DictReader(io.open(os.path.join(B, FN), encoding="utf8"), delimiter=TAB))


def head(p):
    return io.open(p, encoding="utf8").read().split(NL)[0]


def tag(h, k, d="-"):
    m = re.search(k + r"=([0-9.]+)", h)
    return m.group(1) if m else d


HDR = ["species_key", "array_id", "n_units_OK", "unit_type_verdict", "n_unit_types",
       "ref_name", "ref_len", "ref_period", "ref_n_units", "ref_make_method", "ref_member_identity"]
out = io.open(os.path.join(B, "UNIT_REF_TABLE.tsv"), "w", encoding="utf8", newline=NL)
out.write(TAB.join(HDR) + NL)
show = []
for r in sorted(rows, key=lambda x: (x["species_key"], x["array_id"])):
    v = r["unit_type_verdict"]
    base = dict(species_key=r["species_key"], array_id=r["array_id"], n_units_OK=r["n_units_OK"])
    if r["gate_result"].startswith("통과") and r.get("stage2_result", "").startswith("통과"):
        fs = r["final_group_files"].split(",")
        base["unit_type_verdict"] = "여러 종류"
        base["n_unit_types"] = len(fs)
        ent = []
        for f in fs:
            h = head(os.path.join(GD, f))
            ent.append(dict(base, ref_name=f.split(".")[-3], ref_len=tag(h, "len"),
                            ref_period=tag(h, "period"), ref_n_units=tag(h, "n"),
                            ref_make_method=(re.search(r"method=(\w+)", h).group(1)
                                             if "method=" in h else "-"),
                            ref_member_identity=tag(h, "member_median_identity")))
        for e in ent:
            out.write(TAB.join(str(e.get(k, "-")) for k in HDR) + NL)
        show.append((base, ent))
        continue
    if v.startswith("제외"):
        base.update(unit_type_verdict=v, n_unit_types="0")
    elif v.startswith("판정불가"):
        base.update(unit_type_verdict="판정불가(OK 유닛<20)", n_unit_types="?")
    else:
        base.update(unit_type_verdict="한 종류", n_unit_types="1")
    rl = r["refunit_len"]
    if rl in ("-", ""):
        # an array absent from the reference-unit parameter table, such as a zebrafish .k2 sub-map, is measured from the fasta directly
        fp = os.path.join(os.path.dirname(B), "01_canonical", r["species_key"],
                          r["array_id"] + ".unit.fa")
        if os.path.exists(fp):
            rl = str(len("".join(x.strip() for x in
                                 io.open(fp, encoding="utf8").read().split(NL)[1:])))
    e = dict(base, ref_name="arrayref", ref_len=rl, ref_period="-",
             ref_n_units=r["n_units_OK"], ref_make_method="modal_consensus", ref_member_identity="-")
    if r["gate_result"].startswith("불통과"):
        e["unit_type_verdict"] = "한 종류(1 단계 탈락)"
        e["n_unit_types"] = "1"
    elif r.get("stage2_result", "").startswith("사용자"):
        e["unit_type_verdict"] = "한 종류(그룹핑 철회)"
        e["n_unit_types"] = "1"
    elif r.get("stage2_result", "").startswith("불통과"):
        e["unit_type_verdict"] = "한 종류(2 단계 탈락·연속체)"
        e["n_unit_types"] = "1"
    out.write(TAB.join(str(e.get(k, "-")) for k in HDR) + NL)
    show.append((e, None))
out.close()
print("%-22s %-32s %5s %-16s %4s %-9s %6s %6s %5s %-16s %s"
      % ("종", "array", "유닛", "판정", "kind", "참조", "길이", "주기", "유닛", "method", "구성원일치도"))
for b, ent in show:
    if ent:
        for i, e in enumerate(ent):
            print("%-22s %-32s %5s %-16s %4s %-9s %6s %6s %5s %-16s %s"
                  % (b["species_key"][:22] if i == 0 else "", b["array_id"].split("__")[1][:32] if i == 0 else "",
                     b["n_units_OK"] if i == 0 else "", b["unit_type_verdict"] if i == 0 else "",
                     b["n_unit_types"] if i == 0 else "", e["ref_name"], e["ref_len"], e["ref_period"],
                     e["ref_n_units"], e["ref_make_method"], e["ref_member_identity"]))
    else:
        print("%-22s %-32s %5s %-16s %4s %-9s %6s %6s %5s %-16s %s"
              % (b["species_key"][:22], b["array_id"].split("__")[1][:32], b["n_units_OK"],
                 b["unit_type_verdict"], b["n_unit_types"], b["ref_name"], b["ref_len"],
                 b["ref_period"], b["ref_n_units"], b["ref_make_method"], b["ref_member_identity"]))
