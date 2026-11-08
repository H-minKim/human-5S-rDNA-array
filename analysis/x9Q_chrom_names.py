# -*- coding: utf-8 -*-
# 접근번호 → 염색체 번호 표
#   자료 폴더에는 유인원·애기장대의 염색체 이름이 없고 접근번호만 있었다(제브라·금화조는 이름이 있다).
#   NCBI eutils esummary(db=nuccore) 로 확인해 적는다 — 추정하지 않는다. 조회일 2026-10-06.
#     https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=<accession>&retmode=json
#   검산 : 같은 조회에서 받은 서열 길이(slen)가 우리 표의 염색체길이와 글자까지 같아야 한다.
#   산출 : tables/CHROM_NAMES.tsv
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
REASON = "접근번호→염색체 번호(NCBI esummary 2026-10-06)"
# (접근번호, 종, NCBI 가 적은 organism, 염색체, NCBI slen)
ROWS = [("CM054679.2", "수마트라오랑우탄", "Pongo abelii", "chr1", 233795519),
        ("CM054630.2", "보르네오오랑우탄", "Pongo pygmaeus", "chr1", 230807400),
        ("CM055446.2", "고릴라", "Gorilla gorilla gorilla", "chr1", 243847345),
        ("CM055472.2", "보노보", "Pan paniscus", "chr1", 227705052),
        ("CM054434.2", "침팬지", "Pan troglodytes", "chr1", 231299448),
        ("CP116280.1", "애기장대", "Arabidopsis thaliana", "chr1", 32637894),
        ("CP116281.2", "애기장대", "Arabidopsis thaliana", "chr2", 27782540),
        ("CP116282.1", "애기장대", "Arabidopsis thaliana", "chr3", 26149336),
        ("CP116283.2", "애기장대", "Arabidopsis thaliana", "chr4", 25768865),
        ("CP116284.1", "애기장대", "Arabidopsis thaliana", "chr5", 30142610)]
T = pd.DataFrame(ROWS, columns=["접근번호", "종", "NCBI_organism", "염색체번호", "NCBI_slen"])
T["출처"] = "NCBI eutils esummary db=nuccore (조회 2026-10-06)"
M = pd.read_csv(NEW + "/tables/FIGDATA_A1_chrom_map.tsv", sep="\t")
LEN = dict(zip(M.염색체이름.astype(str), M.염색체길이))
T["우리표_길이"] = T.접근번호.map(LEN)
T["길이일치"] = [("확인" if (pd.notna(a) and pd.notna(b) and int(a) == int(b)) else
                ("우리표에 없음" if pd.isna(b) else "조회값 없음")) for a, b in zip(T.NCBI_slen, T.우리표_길이)]
TMP = os.environ.get("TEMP", ".") + "/x9Q"
os.makedirs(TMP, exist_ok=True)
T.to_csv(TMP + "/CHROM_NAMES.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/CHROM_NAMES.tsv", "tables/CHROM_NAMES.tsv", REASON)
print(T[["접근번호", "종", "염색체번호", "NCBI_slen", "우리표_길이", "길이일치"]].to_string(index=False))
print()
print("유인원 5 종은 모두 chromosome 1 — 사람 5S array 가 있는 염색체와 같다")
