# -*- coding: utf-8 -*-
# 종마다 5S 유전자 22-mer 4 개 + 역상보 = k-mer 8 줄 만들기 (선별용)
#   재료 = array_struct.py 가 뽑은 <종>.units.fa 의 유닛들. 유닛 앞쪽 119 bp 가 유전자(hit 구간)다
#   여러 유닛에서 다수결을 취해 개별 유닛의 변이를 피한다
#   인자 : units.fa · 유전자 길이(기본 119) · 출력파일
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import sys
import collections

NL = chr(10)
FA, GL, OUT = sys.argv[1], int(sys.argv[2] if len(sys.argv) > 2 else 119), sys.argv[3]
seqs, cur = [], []
for L in io.open(FA, encoding="utf8"):
    if L.startswith(">"):
        if cur:
            seqs.append("".join(cur))
        cur = []
    else:
        cur.append(L.strip())
if cur:
    seqs.append("".join(cur))
genes = [s[:GL].upper() for s in seqs if len(s) >= GL]
cons = []
for i in range(GL):
    c = collections.Counter(g[i] for g in genes if g[i] in "ACGT")
    cons.append(c.most_common(1)[0][0] if c else "N")
gene = "".join(cons)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
ks = [gene[i:i + 22] for i in (12, 39, 62, 84)]
o = io.open(OUT, "w", encoding="utf8", newline=NL)
for k in ks:
    if "N" in k:
        continue
    o.write(k + NL)
    o.write("".join(COMP[x] for x in reversed(k)) + NL)
o.close()
print("유닛 %d 개에서 합의 유전자 %d bp · k-mer %d 줄" % (len(genes), len(gene), 2 * len([k for k in ks if "N" not in k])))
print(gene)
