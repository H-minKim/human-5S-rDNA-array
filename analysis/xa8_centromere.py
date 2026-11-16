# -*- coding: utf-8 -*-
# Measuring the centromere position.
#
#   The available data differs between species, so the measurement is made by three methods, and which method was used is recorded in the table.
#
#   1) The five great apes: an annotated centromeric region.
#      The Eichler Lab T2T primate centromeric_regions RepeatMasker track, in chromosome coordinates.
#      https://eichlerlab.gs.washington.edu/help/glogsdon/t2t_primates/frozen_repaired_assemblies_repeatmasker/centromeric_regions/
#      The extent of the Satellite interval on chr1, with the two haplotypes merged: the NCBI title does not say which haplotype our accession is, so both were merged.
#
#   2) Arabidopsis Col-CC (GCA_028009825.2): measured directly from the sequence.
#      The NCBI genomic.fna.gz was downloaded and, in each 50 kb window, the autocorrelation at periods of 160-200 bp was measured over the first 5 kb;
#      windows with an identity of at least 0.80 and a period of 170-186 bp, that is CEN180 satellite, were collected and the longest run taken as the centromere.
#      The measured period is 178 bp on all five chromosomes, which matches CEN180.
#
#   3) Zebrafish GRCz12tu (GCA_049306965.1): a peak in satellite density.
#      Two separate steps use two records of this assembly:
#        a) the assembly the reads were aligned to is the GenBank GCA_049306965.1;
#           `node2/zfwgbs/zfB_setup.sh:10-12` fetches that accession through the NCBI datasets API
#           (`.../genome/accession/GCA_049306965.1/download?include_annotation_type=GENOME_FASTA`);
#        b) the satellite track is taken from the UCSC GenArk hub for GCF_049306965.2,
#           GenArk having no hub for GCA_049306965.1.
#      Downloaded from:
#        https://hgdownload.soe.ucsc.edu/hubs/GCF/049/306/965/GCF_049306965.2/GCF_049306965.2.repeatMasker.out.gz
#      The track was matched by chromosome name rather than by accession, which is why the accession column of the table holds entries of the form chr1.
#      The two records differ in release, .1 against .2, so a) and b) are stated separately in the Methods.
#      From that hub's repeatMasker.out.gz the Satellite annotations are collected into 100 kb bins,
#      and the longest run of bins at least 50% satellite is taken. Over that run the satellite density is 40 to 259 times the genome median.
#
#   4) Zebra finch bTaeGut7.mat (GCA_048771995.1): measured directly from the sequence.
#      The assembly downloaded and used is the GenBank GCA_048771995.1, not the RefSeq record:
#      of the 103,811 contig accessions across the 156 result tables, every one is CM, CP or OZ, and none is NC_, the RefSeq form.
#      A GCA accession and a GCF accession of the same assembly can carry different release numbers
#      (GRCz12tu is GCA_049306965.1 against GCF_049306965.2).
#      A search over periods up to 800 bp with a threshold of 0.70-0.80 misses the zebra finch centromeric satellite.
#      The literature on the complete zebra finch genome gives the centromeric satellites as Tgut716A, 716 bp, and Tgut191A, 191 bp.
#      The period range is therefore 150-1600 bp and the threshold 0.55; every 10 kb window is scanned,
#      and only windows whose period is a multiple of 191 or 716 are taken as centromeric satellite and joined where within 150 kb of each other.
#      The result: chr2 78.05-78.41 Mb (monomer 191), chr4 55.87-56.42 Mb (191, partly 716) and chr9 0.01-0.20 Mb (716).
#      chr9 lies at the very start of the chromosome, at the terminus, which matches the published case of Tgut716A abutting the telomere directly.
#
#   Output: tables/CENTROMERE.tsv
import io, os, sys
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'Environment variables')" % k)
    return v


NEW = _need("RDNA_OUT_DIR")
sys.path.insert(0, NEW + "/code")
import arch  # noqa
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
REASON = "중심절 — 유인원(주석)·애기장대(CEN180 실측)·제브라피시(위성밀도) · 금화조 없음(2026-10-06)"
M1 = "주석된 중심절 영역 (Eichler T2T centromeric_regions · Satellite 범위 · 두 하플로 합침)"
M2 = "서열에서 측정 (CEN180 주기 178 bp 탠덤 덩어리 · 일치율>=0.80)"
M3 = "위성 밀도 봉우리 (UCSC GenArk repeatMasker Satellite >=50%/100kb 연속 · 유전체 중앙의 40~259 배)"
M4 = ("서열에서 측정 (Tgut191A 191 bp · Tgut716A 716 bp 의 배수 주기 · 일치율>=0.55 · 10 kb 창) — 앞선 '못 찾음' 은 주기 범위를 400~800 bp 까지만 본 탐지 한계였다(2026-10-06 정정)")
ROWS = [
    ("수마트라오랑우탄", "CM054679.2", "chr1", 101370000, 110470000, M1),
    ("보르네오오랑우탄", "CM054630.2", "chr1", 101160000, 106560000, M1),
    ("고릴라", "CM055446.2", "chr1", 107540000, 117380000, M1),
    ("보노보", "CM055472.2", "chr1", 100210000, 101060000, M1),
    ("침팬지", "CM054434.2", "chr1", 99830000, 107540000, M1),
    ("애기장대", "CP116280.1", "chr1", 14850000, 16550000, M2),
    ("애기장대", "CP116281.2", "chr2", 9400000, 11600000, M2),
    ("애기장대", "CP116282.1", "chr3", 13750000, 15750000, M2),
    ("애기장대", "CP116283.2", "chr4", 8400000, 9900000, M2),
    ("애기장대", "CP116284.1", "chr5", 13400000, 14000000, M2),
    ("제브라피시", "chr1", "chr1", 30200000, 32000000, M3),
    ("제브라피시", "chr3", "chr3", 46100000, 48000000, M3),
    ("제브라피시", "chr4", "chr4", 27600000, 28700000, M3),
    ("제브라피시", "chr7", "chr7", 54600000, 56100000, M3),
    ("제브라피시", "chr10", "chr10", 19100000, 21700000, M3),
    ("제브라피시", "chr13", "chr13", 16000000, 17500000, M3),
    ("제브라피시", "chr14", "chr14", 22500000, 24400000, M3),
    ("제브라피시", "chr18", "chr18", 33200000, 34500000, M3),
    ("금화조", "CM109790.1", "chr2", 78050000, 78410000, "M4"),
    ("금화조", "CM109792.1", "chr4", 55870000, 56420000, "M4"),
    ("금화조", "CM109798.1", "chr9", 10000, 200000, "M4"),
]
T = pd.DataFrame(ROWS, columns=["종", "접근번호", "염색체", "시작", "끝", "방법"])
T["방법"] = T.방법.replace({"M4": M4})
T["시작Mb"] = (T.시작 / 1e6).round(2)
T["끝Mb"] = (T.끝 / 1e6).round(2)
T["길이Mb"] = ((T.끝 - T.시작) / 1e6).round(2)
OUT = T
TMP = os.environ.get("TEMP", ".") + "/xa8"
os.makedirs(TMP, exist_ok=True)
OUT.to_csv(TMP + "/CENTROMERE.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/CENTROMERE.tsv", "tables/CENTROMERE.tsv", REASON)
print(OUT[["종", "염색체", "시작Mb", "끝Mb", "길이Mb"]].to_string(index=False))
print("\n방법 : ①유인원 %s\n       ②애기장대 %s\n       ③제브라피시 %s" % (M1, M2, M3))
