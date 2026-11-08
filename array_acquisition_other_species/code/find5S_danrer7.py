# -*- coding: utf-8 -*-
# danRer7 염색체에서 모계형(chr4)·체세포형(chr18) 5S 유전자를 서열로 찾는다(좌표 옮기기 대신 서열 검색).
#   두 유형 119 bp 를 양가닥으로 미끄러뜨려 불일치 ≤12 인 자리를 모으고, 더 가까운 유형으로 배정한다.
#   빠르게 : 두 유형이 공유하는 뒤쪽 60 bp(54–119, 차이 없음) 의 16-mer 씨앗으로 후보만 본다.
#   인자 : chr.fa.gz 출력.tsv
import gzip
import io
import sys

TYPES = {}
name = None
for L in io.open(sys.argv[3] if len(sys.argv) > 3 else "../gsc/zf5S_types.fa"):
    L = L.strip()
    if L.startswith(">"):
        name = L[1:]
    elif L:
        TYPES[name] = L
RC = str.maketrans("ACGTN", "TGCAN")
seq = []
chrom = None
with gzip.open(sys.argv[1], "rt") as f:
    for L in f:
        if L.startswith(">"):
            chrom = L[1:].strip().split()[0]
        else:
            seq.append(L.strip().upper())
g = "".join(seq)
N = len(g)
mat = TYPES["maternal_chr4"]
som = TYPES["somatic_chr18"]
seeds = set()
for t in (mat, som):
    for i in range(60, 119 - 16 + 1):
        seeds.add(t[i:i + 16])
SEEDPOS = mat.find("GGGAGACCGCCTGGGA")
assert SEEDPOS == som.find("GGGAGACCGCCTGGGA") and SEEDPOS >= 0
out = io.open(sys.argv[2], "a", encoding="utf8", newline="\n")
seen = set()
hits = 0


def mm(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


for strand in "+-":
    s = g if strand == "+" else g.translate(RC)[::-1]
    i = 0
    while True:
        # 씨앗 하나(가장 보존된 곳)만으로 빠르게 찾는다
        j = s.find("GGGAGACCGCCTGGGA", i)
        if j < 0:
            break
        i = j + 1
        # 씨앗 위치는 서열에서 직접 계산한다(가정 금지 — 처음 87 로 넣어 전부 탈락했다)
        st = j - SEEDPOS
        if st < 0 or st + 119 > len(s):
            continue
        w = s[st:st + 119]
        dm, ds = mm(w, mat), mm(w, som)
        best = min(dm, ds)
        if best > 12:
            continue
        typ = "maternal" if dm < ds else ("somatic" if ds < dm else "tie")
        if strand == "+":
            a, b = st + 1, st + 119
        else:
            a, b = N - (st + 119) + 1, N - st
        if (a, b) in seen:
            continue
        seen.add((a, b))
        hits += 1
        out.write("\t".join([chrom, str(a), str(b), strand, typ, str(dm), str(ds)]) + "\n")
out.close()
print(chrom, "히트", hits)
