# -*- coding: utf-8 -*-
# 번호 옮기기(relabel) 되돌림 — 옛 build_ao_final.py 와 같은 규칙 (2026-09-13)
#   왜 : ao_g8cut 은 저장 배치와 어긋나 보이는 read 를 저장 번호로 '옮겨' 준다(relabeled_d).
# The number is the real position; forcing the stored value is wrong (measured over 66 individuals, 48 against 19).
#        옛 파이프라인은 build_ao_final.py 에서 이 옮김을 되돌렸다 :
#            r[7] = r[7] - d   (array_pos)
#            r[8] = r[8] - d   (G6 번호)
#        내가 이 단계를 빠뜨려 새 판에 57 read 의 번호가 옮겨진 채 들어갔다.
#   같이 하는 것 : 수리 슬롯 2 곳의 좌표틀 보정(fix_shift2 와 같은 규칙) — 둘은 서로 다른 원인이라 함께 적용한다.
#   원본은 덮지 않는다.
#   인자 : 원 유닛표(ao_g8cut 산출) · AO_G8_READS.tsv · 출력
import io
import re
import sys

NL = chr(10)
TAB = chr(9)
INP, READF, OUTP = sys.argv[1], sys.argv[2], sys.argv[3]
SHIFT = {}

RE = re.compile(r"relabeled_d=([+-]\d+)")
D = {}
for L in io.open(READF, encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) < 11:
        continue
    m = RE.search(r[10])
    if m:
        D[(r[0], r[1])] = int(m.group(1))
print("relabel 된 read %d" % len(D))

n = rel = sh = 0
O = io.open(OUTP, "w", encoding="utf8", newline=NL)
for L in io.open(INP, encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    n += 1
    if len(r) >= 15:
        d = D.get((r[0], r[2]))
        if d and r[7].lstrip("-").isdigit() and r[8].lstrip("-").isdigit():
            r[7] = str(int(r[7]) - d)
            r[8] = str(int(r[8]) - d)
            rel += 1
        k = (r[0], r[3])
        if k in SHIFT and r[7].lstrip("-").isdigit() and int(r[7]) > SHIFT[k]:
            r[7] = str(int(r[7]) - 1)
            sh += 1
    O.write(TAB.join(r) + NL)
O.close()
print("%s -> %s · 전체 %d · relabel 되돌린 유닛 %d · 수리틀 보정 유닛 %d" % (INP, OUTP, n, rel, sh))
