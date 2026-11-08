# -*- coding: utf-8 -*-
# 중심절 위치 실측
#
#   종마다 자료가 달라 방법을 셋으로 나눠 실측 했고, 어느 방법인지 표에 함께 적는다.
#
#   ① 유인원 5 종 — 주석된 중심절 영역
#      Eichler Lab T2T 영장류 centromeric_regions RepeatMasker(염색체 좌표)
#      https://eichlerlab.gs.washington.edu/help/glogsdon/t2t_primates/frozen_repaired_assemblies_repeatmasker/centromeric_regions/
#      chr1 의 Satellite 구간 범위(두 하플로 합침). 우리 접근번호가 어느 하플로인지 NCBI 제목에 없어 둘을 합쳤다.
#
#   ② 애기장대 Col-CC (GCA_028009825.2) — 서열에서 직접 측정
#      NCBI genomic.fna.gz 를 받아 50 kb 창마다 앞 5 kb 로 주기 160~200 bp 자기상관을 재고,
#      일치율>=0.80 이고 주기 170~186 bp 인 창(= CEN180 위성)을 모아 가장 긴 덩어리 를 중심절로 잡았다.
#      실측 주기는 다섯 염색체 모두 178 bp (CEN180 과 같다).
#
#   ③ 제브라피시 GRCz12tu (GCA_049306965.1) — 위성 밀도 봉우리
#      정정 : 이 줄은 `GCF_049306965.2`(RefSeq·판 .2)라고 적고 있었는데 틀렸다.
#      실제로 받아 쓴 것은 GenBank `GCA_049306965.1` 이다 — 상류 노드 스크립트
# `node2/zfwgbs/zfB_setup.sh:10-12` fetches that accession through the NCBI datasets API
#      (`.../genome/accession/GCA_049306965.1/download?include_annotation_type=GENOME_FASTA`).
#      원고와 카탈로그도 `GCA_049306965.1` 로 적는다.
#      덧붙임 : 그때 *"이 주석만 어긋나 있었다"* 라고 적었는데 그것이 틀렸다 —
#      아래 ④ 의 금화조 줄이 똑같은 흠이었고 네 줄 아래였는데 못 보았다(2026-10-28 정정).
#      정정 : 이 ③ 은 두 걸음 을 한 줄로 적고 있었다. 가려 적는다.
#        (가) 읽기를 붙인 조립본 은 위의 GenBank `GCA_049306965.1` 이 맞다.
#        (나) 위성 트랙을 받은 자리 는 UCSC GenArk 의 `GCF_049306965.2` 허브 다.
#      10-24 에 (가)와 (나)를 한 걸음으로 보고 `GCF` 를 `GCA` 로 함께 바꾼 것이 과한 고침이었다.
#      근거(2026-11-06 눌러 확인) : GenArk 에 `GCA_049306965.1` 허브가 없다.
#        `/hubs/GCA/049/306/` 아래에 `905/` 하나뿐이고 `.../965/` 는 404 다.
#        `/hubs/GCF/049/306/965/GCF_049306965.2/` 는 200 이고 그 안에 repeatMasker.out.gz 가 있다.
#      받은 자리 :
#        https://hgdownload.soe.ucsc.edu/hubs/GCF/049/306/965/GCF_049306965.2/GCF_049306965.2.repeatMasker.out.gz
#      트랙은 접근번호가 아니라 염색체 이름 으로 짝지었다. 그래서 표의 접근번호 칸에 `chr1` 꼴이 들어간다.
#      두 레코드는 판이 `.1` 과 `.2` 로 다르니 Methods 에 (가)와 (나)를 따로 적는다.
#      표의 「방법」 글자(아래 M3)는 고치지 않는다 . 발행된 CENTROMERE.tsv 를 바꾸지 않기 위함이다.
#      그 허브의 repeatMasker.out.gz 에서 Satellite 를 100 kb 칸으로 모아
#      위성이 칸의 50% 이상인 연속 구간 가운데 가장 긴 것. 그 구간의 위성 밀도는 유전체 중앙값의 40~259 배.
#
#   ④ 금화조 bTaeGut7.mat (GCA_048771995.1) — 서열에서 직접 측정
#      정정 : 이 줄은 `GCF_048771995.1`(RefSeq)이라고 적고 있었는데 틀렸다.
#      받아 쓴 것은 GenBank `GCA_048771995.1` 이다 — 출력물이 증거다 :
#      결과표 156 개의 contig 접근번호 103,811 개가 `CM`/`CP`/`OZ` 뿐이고 `NC_`(RefSeq)가 0 개다.
#      `GCA`→`GCF` 로 접두만 바꾸면 안 된다 — 같은 조립본이라도 판이 다르다
#      (`GRCz12tu` 는 `GCA_049306965.1` ↔ `GCF_049306965.2`).
# The earlier negative result came from the limits of the search, not from absence :
#        주기를 400(뒤에는 800) bp 까지만 보고 문턱을 0.70~0.80 으로 둬서, 금화조 중심절 위성을 놓쳤다.
#      문헌(완전 금화조 유전체) : 중심절 위성은 Tgut716A(716 bp) 와 Tgut191A(191 bp) 다.
#      그래서 주기 범위를 150~1600 bp 로 넓히고 문턱을 0.55 로 낮춰 10 kb 창을 빠짐없이 훑은 뒤,
#      주기가 191 또는 716 의 배수인 창 만 중심절 위성으로 보고 150 kb 안에 있는 것끼리 이었다.
#      결과 : chr2 78.05~78.41 Mb(단량체 191) · chr4 55.87~56.42 Mb(191, 일부 716) · chr9 0.01~0.20 Mb(716)
#      — chr9 은 염색체 맨 앞(말단) 에 붙어 있어 문헌이 적은 "Tgut716A 가 텔로미어에 바로 붙는 경우" 와 맞는다.
#
#   산출 : tables/CENTROMERE.tsv
import io, os, sys
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
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
