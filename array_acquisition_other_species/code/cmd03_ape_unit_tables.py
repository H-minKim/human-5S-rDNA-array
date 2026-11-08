# -*- coding: utf-8 -*-
# 유인원 5S array 유닛별 변이표 (바깥군 대조용)
#   변이 토큰은 전부 **그 종의 정본 유닛**(02_reference_units) 대비다. 사람 정본 대비가 아니다.
#   ① READLEVEL : read(=분자) × 유닛. 원본 06_units_from_reads/*.read_unit_seq.tsv 의 varlist 를 그대로 옮긴다
#      (mm2 로 종 정본 유닛에 정렬한 뒤 뽑은 토큰 · 형식 `<위치><염기>` 치환 · `<위치>I<n>`·`<위치>D<n>` 삽입/결실)
#   ② CONSENSUS : 유닛 자리(array_pos) × read 합의. 07_consensus_from_reads 의 varlist_cons · varlist_strict
#   ③ COORDS : array_pos ↔ 염색체 좌표 · 가닥 (방향 규약 근거)
#   ④ REFUNITS : 종별 정본 유닛 FASTA(+ 사람 CHM13)
#   유닛 수 규약 : 온전 유닛만 · array_pos 1 = 첫 온전 유닛(v11 과 같음). 양 끝 부분 유닛은 세지 않았다.
# MIT License · Copyright (c) 2026 Hungmin Kim
import csv
import glob
import io
import os

R = "$RDNA_OTHERSP_DIR/data_v2/"
OUT = R + "14_ape_unit_variants/"
os.makedirs(OUT, exist_ok=True)
T = "\t"
APES = ["bonobo_mPanPan1", "chimp_mPanTro3", "gorilla_mGorGor1", "orang_bornean_mPonPyg2",
        "orang_sumatran_mPonAbe1", "siamang_mSymSyn1"]
REFLEN = {}
for r in csv.DictReader(io.open(R + "02_reference_units/REFERENCE_UNITS_INDEX.tsv", encoding="utf8"), delimiter=T):
    REFLEN[r["species_key"]] = r["refunit_len"]

# ① read 수준
o = io.open(OUT + "APE_UNIT_VARIANTS_READLEVEL.tsv", "w", encoding="utf8", newline="\n")
o.write(T.join(["species", "sample", "molecule_id", "slot", "array_id", "array_pos", "unit_len", "canon_len",
                "spans_unit", "frac_covered", "mapq", "n_var", "variant"]) + "\n")
n1 = 0
for sp in APES:
    for p in sorted(glob.glob(R + "06_units_from_reads/%s__*.read_unit_seq.tsv" % sp)):
        for r in csv.DictReader(io.open(p, encoding="utf8"), delimiter=T):
            o.write(T.join([sp, r["sample"], r["read_id"], "read(분자 · 하플로타입 위상 없음)", r["array_id"],
                            r["array_pos"], r["unit_len"], REFLEN.get(sp, ""), r["spans_unit"], r["frac_covered"],
                            r["mapq"], r["n_var"], r["varlist"]]) + "\n")
            n1 += 1
o.close()

# ② 유닛 자리 합의 수준
o = io.open(OUT + "APE_UNIT_VARIANTS_CONSENSUS.tsv", "w", encoding="utf8", newline="\n")
o.write(T.join(["species", "sample", "slot", "array_id", "array_pos", "canon_len", "depth_median", "n_reads_span",
                "n_var_cons", "variant_cons", "n_var_strict", "variant_strict", "n_var_asm", "concordance_strict"]) + "\n")
n2 = 0
for sp in APES:
    for p in sorted(glob.glob(R + "07_consensus_from_reads/%s__*.unit_consensus.tsv" % sp)):
        for r in csv.DictReader(io.open(p, encoding="utf8"), delimiter=T):
            o.write(T.join([sp, r["sample"], "조립본 1 벌(하플로타입 분리 없음)", r["array_id"], r["array_pos"],
                            r["refunit_len"], r["depth_median"], r["n_reads_span"], r["n_var_cons"], r["varlist_cons"],
                            r["n_var_strict"], r["varlist_strict"], r["n_var_asm"], r["concordance_strict"]]) + "\n")
            n2 += 1
o.close()

# ③ 좌표 (방향 근거) — 사람 CHM13 포함
o = io.open(OUT + "APE_UNIT_COORDS.tsv", "w", encoding="utf8", newline="\n")
o.write(T.join(["species", "array_id", "array_pos", "unit_len", "chrom", "chrom_start", "chrom_end",
                "array_strand", "canon_start", "canon_end"]) + "\n")
n3 = 0
for sp in APES + ["human_CHM13v2"]:
    for p in sorted(glob.glob(R + "06_units_from_canonical/%s__*.units.tsv" % sp)):
        for r in csv.DictReader(io.open(p, encoding="utf8"), delimiter=T):
            o.write(T.join([sp, r["array_id"], r["array_pos"], r["unit_len"], r["chrom"], r["chrom_start"],
                            r["chrom_end"], r["strand"], r["canon_start"], r["canon_end"]]) + "\n")
            n3 += 1
o.close()

# ④ 종별 정본 유닛 FASTA
keep, name = [], None
for L in io.open(R + "02_reference_units/REFERENCE_UNITS_ALL.fa", encoding="utf8"):
    if L.startswith(">"):
        name = L[1:].split()[0].split("__")[0]
    if name in APES or name == "human_CHM13v2":
        keep.append(L.rstrip("\n"))
io.open(OUT + "APE_REFERENCE_UNITS.fa", "w", encoding="utf8", newline="\n").write("\n".join(keep) + "\n")
print("READLEVEL %d · CONSENSUS %d · COORDS %d · FASTA %d 줄" % (n1, n2, n3, len(keep)))
