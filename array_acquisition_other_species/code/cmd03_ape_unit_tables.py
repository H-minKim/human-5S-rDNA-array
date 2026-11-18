# -*- coding: utf-8 -*-
# Per-unit variant tables for the great ape 5S arrays, for use as an outgroup comparison.
#   Every variant token is relative to that species' own reference unit, from 02_reference_units, and never to the human canonical unit.
#   1) READLEVEL: read, that is molecule, by unit. The varlist of the source 06_units_from_reads/*.read_unit_seq.tsv is carried across unchanged
#      (the tokens taken after aligning to that species' reference unit with minimap2, in the form <position><base> for a substitution and <position>I<n> or <position>D<n> for an insertion or deletion).
#   2) CONSENSUS: unit position (array_pos) by read consensus, taking varlist_cons and varlist_strict from 07_consensus_from_reads.
#   3) COORDS: array_pos against chromosome coordinate and strand, which is the basis of the orientation convention.
#   4) REFUNITS: the reference unit FASTA per species, plus human CHM13.
#   Unit-count convention: complete units only, with array_pos 1 the first complete unit, as in v11. The partial units at either end are not counted.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import csv
import glob
import io
import os

R = os.environ["RDNA_OTHERSP_DIR"] + "/data_v2/"
OUT = R + "14_ape_unit_variants/"
os.makedirs(OUT, exist_ok=True)
T = "\t"
APES = ["bonobo_mPanPan1", "chimp_mPanTro3", "gorilla_mGorGor1", "orang_bornean_mPonPyg2",
        "orang_sumatran_mPonAbe1", "siamang_mSymSyn1"]
REFLEN = {}
for r in csv.DictReader(io.open(R + "02_reference_units/REFERENCE_UNITS_INDEX.tsv", encoding="utf8"), delimiter=T):
    REFLEN[r["species_key"]] = r["refunit_len"]

# 1) the read level
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

# 2) the per-unit-position consensus level
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

# 3) the coordinates, the basis of the orientation, including human CHM13
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

# 4) the reference unit FASTA per species
keep, name = [], None
for L in io.open(R + "02_reference_units/REFERENCE_UNITS_ALL.fa", encoding="utf8"):
    if L.startswith(">"):
        name = L[1:].split()[0].split("__")[0]
    if name in APES or name == "human_CHM13v2":
        keep.append(L.rstrip("\n"))
io.open(OUT + "APE_REFERENCE_UNITS.fa", "w", encoding="utf8", newline="\n").write("\n".join(keep) + "\n")
print("READLEVEL %d, CONSENSUS %d, COORDS %d, FASTA %d lines" % (n1, n2, n3, len(keep)))
