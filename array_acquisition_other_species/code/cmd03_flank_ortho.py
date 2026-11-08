# -*- coding: utf-8 -*-
# "array_pos 1 쪽 flank" 가 사람의 어느 쪽 flank 와 상동인가 — 공유 k-mer 로 판정
#   정본 fasta 는 [flank5 20kb][array][flank3 20kb] 이고 유전자(전사) 방향으로 저장돼 있다.
#   그러므로 flank5 = array_pos 1 쪽이다(06_units_from_canonical 의 canon_start 로 확인: pos1 ≈ 22 kb).
#   판정 : 종 flank5 의 k-mer 집합을 사람 flank5 / flank3 집합과 견준다(양 가닥 · 자카드와 적중수 둘 다).
# A call needs one side clearly larger than the other; otherwise it is recorded as unknown.
# MIT License · Copyright (c) 2026 Hungmin Kim
import glob
import io
import os
import sys

K = 21
FL = 20000
S = os.path.dirname(os.path.abspath(__file__))
SRC = sys.argv[1]          # 정본 fasta 사본 폴더 (C:)
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fa(p):
    s = []
    for L in io.open(p, encoding="utf8"):
        if not L.startswith(">"):
            s.append(L.strip().upper())
    return "".join(s)


def kset(s):
    out = set()
    for i in range(len(s) - K + 1):
        w = s[i:i + K]
        if "N" in w:
            continue
        out.add(min(w, rc(w)))        # 가닥 중립
    return out


seqs = {}
for p in sorted(glob.glob(SRC + "/*.fa")):
    sp = os.path.basename(p).split("__")[0]
    s = read_fa(p)
    seqs[sp] = (kset(s[:FL]), kset(s[-FL:]))
    print("%-26s len %d" % (sp, len(s)))

h5, h3 = seqs["human_CHM13v2"]
o = io.open(S + "/FLANK_ORTHOLOGY.tsv", "w", encoding="utf8", newline="\n")
o.write("\t".join(["species", "side_tested", "n_kmer", "hit_human_flank5(array_pos1쪽)", "hit_human_flank3(반대쪽)",
                   "jaccard_f5", "jaccard_f3", "ratio_f5_over_f3", "call"]) + "\n")
for sp in sorted(seqs):
    if sp == "human_CHM13v2":
        continue
    for side, ks in (("flank5(array_pos 1 쪽)", seqs[sp][0]), ("flank3(반대쪽)", seqs[sp][1])):
        a, b = len(ks & h5), len(ks & h3)
        ja = a / max(1, len(ks | h5))
        jb = b / max(1, len(ks | h3))
        if max(a, b) < 300:
            call = "판정불가(공유 k-mer 너무 적음 · 300 미만)"
        elif a >= 3 * max(1, b):
            call = "사람 flank5 쪽(= 사람 array_pos 1 쪽)"
        elif b >= 3 * max(1, a):
            call = "사람 flank3 쪽(= 사람 반대쪽)"
        else:
            call = "모름(양쪽이 비슷)"
        o.write("\t".join([sp, side, str(len(ks)), str(a), str(b), "%.4f" % ja, "%.4f" % jb,
                           "%.2f" % (a / max(1, b)), call]) + "\n")
o.close()
print(io.open(S + "/FLANK_ORTHOLOGY.tsv", encoding="utf8").read())
