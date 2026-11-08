# -*- coding: utf-8 -*-
# 실측 : ① array−flank(ONT, AVF.tsv) 종·array 별 ② 종별 array 수 정의별 ③ 기본값 정의 확인
# 결과 : data_v2/00_catalog/CMD01_*.tsv (새 파일만 씀)
# MIT License · Copyright (c) 2026 Hungmin Kim
import collections
import csv
import glob
import io
import os
import statistics

R = "$RDNA_OTHERSP_DIR/data_v2/"
OUT = R + "00_catalog/"
TAB = "\t"


def rows(p):
    with io.open(p, encoding="utf8") as f:
        return list(csv.DictReader(f, delimiter=TAB))


# ── ① array − flank ────────────────────────────────────────────────
avf = rows(R + "09_array_vs_flank/AVF.tsv")
d = collections.defaultdict(dict)   # (sp, sample, array) -> region -> (n_c, n_m)
for r in avf:
    k = (r["species_key"], r["sample"], r["array_id"])
    d[k][r["region"]] = (int(r["n_c"]), int(r["n_5mC"]), int(r["reads"]))
o = io.open(OUT + "CMD01_ARRAY_MINUS_FLANK.tsv", "w", encoding="utf8", newline="\n")
o.write(TAB.join(["species", "sample", "array_id", "array_frac", "flank_frac(L+R)", "array_minus_flank", "array_reads"]) + "\n")
per_sp = collections.defaultdict(list)
for k, reg in sorted(d.items()):
    if "array" not in reg:
        continue
    a = reg["array"]
    fl = [v for kk, v in reg.items() if kk.startswith("flank")]
    if not fl or a[0] == 0:
        continue
    fc = sum(x[0] for x in fl)
    fm = sum(x[1] for x in fl)
    if fc == 0:
        continue
    af = a[1] / a[0]
    ff = fm / fc
    o.write(TAB.join([k[0], k[1], k[2], "%.3f" % af, "%.3f" % ff, "%+.3f" % (af - ff), str(a[2])]) + "\n")
    per_sp[k[0]].append((k[2], k[1], af - ff))
o.close()
print("① 종별 array−flank (array×시료)")
for sp, v in sorted(per_sp.items()):
    arrays = sorted(set(x[0] for x in v))
    vals = [x[2] for x in v]
    print("  %-26s array %2d · 행 %3d · 범위 %+.3f ~ %+.3f · 중앙 %+.3f" % (sp, len(arrays), len(v), min(vals), max(vals), statistics.median(vals)))

# ── ② array 수 정의별 ───────────────────────────────────────────────
cnt = collections.defaultdict(dict)


def add(col, sp, n):
    cnt[sp][col] = n


for r in rows(OUT + "5S_array_catalog.tsv"):
    pass
cat = rows(OUT + "5S_array_catalog.tsv")
spk = [k for k in cat[0].keys() if k in ("species_key", "species", "assembly_key")][0]
c = collections.Counter(r[spk] for r in cat)
for sp, n in c.items():
    add("catalog_rows", sp, n)
st = rows(OUT + "5S_array_structure.tsv")
spk2 = [k for k in st[0].keys() if k in ("species_key", "species")][0]
for sp, n in collections.Counter(r[spk2] for r in st).items():
    add("structure_rows", sp, n)
for dd in glob.glob(R + "01_canonical/*/"):
    sp = os.path.basename(os.path.normpath(dd))
    add("canonical_fa", sp, len([f for f in glob.glob(dd + "*.fa") if not f.endswith(".unit.fa")]))
ur = rows(OUT + "UNIT_REF_TABLE.tsv")
for sp, n in collections.Counter(r["species_key"] for r in ur).items():
    add("unit_ref_rows", sp, n)
for sp, n in collections.Counter(r["species_key"] for r in ur if r["unit_type_verdict"].startswith("한 종류")).items():
    add("unit_ref_usable(한종류)", sp, n)
for dd in glob.glob(R + "02_consensus_map/*/"):
    sp = os.path.basename(os.path.normpath(dd))
    add("g8map_files", sp, len(glob.glob(dd + "*.g8map.tsv")))
for r in rows(OUT + "TIER_STATUS_FULL.tsv"):
    add("tier_n_arrays", r["species_key"], r["n_arrays"])
    add("tier_reads_placed", r["species_key"], r["reads_placed"])
for sp, v in per_sp.items():
    add("avf_arrays", sp, len(set(x[0] for x in v)))
for dd in glob.glob(R + "03_reads/*/"):
    sp = os.path.basename(os.path.normpath(dd))
    p = dd + "5S_array_summary.tsv"
    if os.path.exists(p):
        rr = rows(p)
        add("readunit_arrays", sp, len(set(r.get("array_id", "") for r in rr)))
cols = ["catalog_rows", "structure_rows", "canonical_fa", "unit_ref_rows", "unit_ref_usable(한종류)", "g8map_files",
        "tier_n_arrays", "avf_arrays", "readunit_arrays", "tier_reads_placed"]
o = io.open(OUT + "CMD01_ARRAY_COUNT_DEFS.tsv", "w", encoding="utf8", newline="\n")
o.write(TAB.join(["species"] + cols) + "\n")
print("② array 수 정의별")
for sp in sorted(cnt):
    o.write(TAB.join([sp] + [str(cnt[sp].get(k, "-")) for k in cols]) + "\n")
    print("  %-26s " % sp + " ".join("%s=%s" % (k.split("_")[0] + "_" + k.split("_")[-1], cnt[sp].get(k, "-")) for k in cols))
o.close()
