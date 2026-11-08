# -*- coding: utf-8 -*-
# 참조유닛 유전자부 · 스페이서부 쌍별 일치도 (금화조 8 array · 제브라 chr18 vs chr4)
#   입력 : 02_reference_units/REFERENCE_UNITS_ALL.fa (헤더 gene=a-b, 유닛 틀 [스페이서][유전자])
#   일치도 = 전역 정렬(Biopython PairwiseAligner, 일치 +1 · 불일치 −1 · 틈 열기 −2 · 늘이기 −0.5)에서
#            일치 열 수 / 정렬 열 수(틈 포함)
#   CpG 수 : 각 부분의 CG 이핵산 수(메틸 측정 가능 자리 수 참고)
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import itertools
from Bio import Align

R = "$RDNA_OTHERSP_DIR/data_v2/"
OUT = R + "13_array_level_diff/REFUNIT_IDENTITY.tsv"
T = "\t"
seqs = {}
name = None
for L in io.open(R + "02_reference_units/REFERENCE_UNITS_ALL.fa", encoding="utf8"):
    L = L.strip()
    if L.startswith(">"):
        f = L[1:].split()
        name = f[0]
        meta = dict(x.split("=", 1) for x in f[1:])
        seqs[name] = {"meta": meta, "s": ""}
    elif name:
        seqs[name]["s"] += L.upper()

al = Align.PairwiseAligner()
al.mode = "global"
al.match_score, al.mismatch_score = 1, -1
al.open_gap_score, al.extend_gap_score = -2, -0.5


def ident(a, b):
    aln = al.align(a, b)[0]
    s1, s2 = aln[0], aln[1]
    m = sum(1 for x, y in zip(s1, s2) if x == y and x != "-")
    return m / len(s1)


def parts(k):
    s = seqs[k]["s"]
    g0, g1 = (int(x) for x in seqs[k]["meta"]["gene"].split("-"))
    return s[:g0 - 1], s[g0 - 1:g1]


FL = {"CM109790.1_96427419_96511119": "chr2A_mat", "CM109751.1_96489183_96572883": "chr2A_pat",
      "CM109790.1_97181170_97425313": "chr2B_mat", "CM109751.1_97242934_97487077": "chr2B_pat",
      "CM109792.1_6514976_6563661": "chr4_mat", "CM109753.1_6524260_6552707": "chr4_pat",
      "CM109798.1_6937339_6952337": "chr9_mat", "CM109759.1_7075841_7101562": "chr9_pat"}
o = io.open(OUT, "w", encoding="utf8", newline="\n")
o.write(T.join(["species", "a", "b", "class", "spacer_ident", "gene_ident", "spacer_len_a", "spacer_len_b",
                "spacer_CpG_a", "spacer_CpG_b", "gene_CpG_a", "gene_CpG_b"]) + "\n")
fk = sorted((k for k in seqs if k.startswith("zfinch_v03_dip__")), key=lambda k: FL[k.split("__")[1]])
for a, b in itertools.combinations(fk, 2):
    la, lb = FL[a.split("__")[1]], FL[b.split("__")[1]]
    ca, cb = la.split("_")[0], lb.split("_")[0]
    cls = "같은 좌위 mat/pat" if ca == cb else ("chr9 vs 나머지" if "chr9" in (ca, cb) else "chr9 밖 좌위끼리")
    (sa, ga), (sb, gb) = parts(a), parts(b)
    o.write(T.join(["zfinch", la, lb, cls, "%.3f" % ident(sa, sb), "%.3f" % ident(ga, gb), str(len(sa)), str(len(sb)),
                    str(sa.count("CG")), str(sb.count("CG")), str(ga.count("CG")), str(gb.count("CG"))]) + "\n")
zk = [k for k in seqs if k.startswith("zebrafish_GRCz12tu__")]
z18 = [k for k in zk if "CM111196.1" in k]
OKZ = set()
for L in io.open(R + "02_reference_units/REFERENCE_UNITS_INDEX.tsv", encoding="utf8"):
    f = L.rstrip("\n").split(T)
    if f[0] == "zebrafish_GRCz12tu" and f[17] == "CONFIRMED" and f[18] == "한 종류":
        OKZ.add(f[1])
z4 = [k for k in zk if "CM111182.1" in k and k in OKZ]   # 한 종류 판정 통과 chr4 array 만(유닛 틀 확실)
for a in z18:
    for b in z4:
        (sa, ga), (sb, gb) = parts(a), parts(b)
        # 스페이서 길이가 크게 다르면(짧은 쪽/긴 쪽 < 0.8) 일치도 대신 '길이다름'
        sid = "%.3f" % ident(sa, sb) if min(len(sa), len(sb)) / max(len(sa), len(sb)) >= 0.8 else "길이다름"
        o.write(T.join(["zebrafish", "chr18:" + a.split("__")[1], "chr4:" + b.split("__")[1], "chr18 vs chr4",
                        sid, "%.3f" % ident(ga, gb), str(len(sa)), str(len(sb)),
                        str(sa.count("CG")), str(sb.count("CG")), str(ga.count("CG")), str(gb.count("CG"))]) + "\n")
o.close()
print("완료", len(fk), len(z18), len(z4))
