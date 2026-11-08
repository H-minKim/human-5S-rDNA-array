# -*- coding: utf-8 -*-
# 타 생물종 5S : 정본(flank-array-flank 대표서열) + consensus map(유닛별 변이표) — 사람 작업과 같은 꼴로
#   정본   : array 양옆 FL bp 를 붙인 서열. array 의 유전자 방향으로 맞춰 둔다(- 가닥이면 역상보)
#   유닛틀 : 사람 정본과 같게 "스페이서 다음 유전자" 순서 — 유닛 i = 앞 유전자 끝+1 .. 이 유전자 끝
#            (사람 canon_v2 는 AluY..NTS2 다음 gene 2122-2240 로 끝난다)
#   유닛번호: array 의 유전자 방향으로 1 부터. 사람의 '원위=1' 규칙은 다른 종에 근거가 없어 쓰지 않는다.
#            대신 유닛마다 염색체 좌표를 같이 적어 어느 끝인지 언제든 되짚게 한다
#   변이   : 사람과 같은 방법 — minimap2 -cx asm5 --cs 로 종별 대표유닛에 맞춰 뽑고 pos+염기 / posD n / posI n 으로 적는다
#   인자 : 종키 fasta tblout struct표 출력폴더 [flank bp=20000]
#   산출 : canon/<species>__<chrom>_<start>_<end>.fa   (정본)
#          map/<species>__<chrom>_<start>_<end>.map.tsv (consensus map)
#          canon_unit/<species>__....unit.fa            (그 array 의 대표 유닛)
#          MAP_SUMMARY.tsv                          (array 한 줄)
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
SP, FA, TBL, STRUCT, OUT = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
FL = int(sys.argv[6]) if len(sys.argv) > 6 else 20000
for d in ("canon", "map", "canon_unit"):
    os.makedirs(os.path.join(OUT, d), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(COMP.get(x, "N") for x in reversed(s))


LEN = {}
for L in io.open(FA + ".fai", encoding="utf8"):
    f = L.split(TAB)
    LEN[f[0]] = int(f[1])
HIT = collections.defaultdict(list)
for L in io.open(TBL, encoding="utf8"):
    if L.startswith("#"):
        continue
    f = L.split()
    a, b = int(f[7]), int(f[8])
    HIT[f[0]].append((min(a, b), max(a, b), f[9]))
for c in HIT:
    HIT[c].sort()
ARR = []
f = io.open(STRUCT, encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ix = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["species"]] != SP:
        continue
    ARR.append((r[ix["chrom"]], int(r[ix["array_start"]]), int(r[ix["array_end"]]), r[ix["strand_major"]][0], int(r[ix["copies"]])))


def faidx(reg):
    r = subprocess.run(["samtools", "faidx", FA, reg], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


def cs_to_vars(cs, rs):
    # minimap2 cs(short) -> 사람 지도와 같은 토큰 : 123A(치환) · 123D2(결실) · 123I3(삽입) · 1-기반 대표유닛 좌표
    out, p, i = [], rs, 0
    while i < len(cs):
        c = cs[i]
        if c == ":":
            j = i + 1
            while j < len(cs) and cs[j].isdigit():
                j += 1
            p += int(cs[i + 1:j])
            i = j
        elif c == "*":
            out.append("%d%s" % (p + 1, cs[i + 2].upper()))
            p += 1
            i += 3
        elif c == "-":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            n = j - i - 1
            out.append("%dD%d" % (p + 1, n))
            p += n
            i = j
        elif c == "+":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            out.append("%dI%d" % (p + 1, j - i - 1))
            i = j
        elif c == "~":
            j = i + 1
            while j < len(cs) and not (cs[j] in ":*-+~"):
                j += 1
            i = j
        else:
            i += 1
    return out


SUM = io.open(os.path.join(OUT, "MAP_SUMMARY_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
SUM.write(TAB.join(["species", "array_id", "chrom", "array_start", "array_end", "strand", "copies", "units_mapped",
                    "canon_unit_len", "units_identical_pct", "var_sites", "var_per_unit_med", "canon_fa", "map_tsv"]) + NL)
for (c, s0, e1, strand, cop) in ARR:
    aid = "%s__%s_%d_%d" % (SP, c, s0, e1)
    L0, R1 = max(1, s0 - FL), min(LEN[c], e1 + FL)
    seq = faidx("%s:%d-%d" % (c, L0, R1))
    if not seq:
        continue
    hits = [h for h in HIT[c] if s0 <= h[0] and h[1] <= e1]
    # 작업 방향 = array 유전자 방향. - 가닥이면 전체를 뒤집는다
    if strand == "-":
        wseq = rc(seq)
        gen = sorted([(len(seq) - (h[1] - L0) - 1, len(seq) - (h[0] - L0) - 1) for h in hits])
    else:
        wseq = seq
        gen = sorted([(h[0] - L0, h[1] - L0) for h in hits])
    io.open(os.path.join(OUT, "canon", aid + ".fa"), "w", encoding="utf8", newline=NL).write(
        ">%s flank%d_array_flank%d strand=%s chrom=%s array=%d-%d%s%s%s" % (aid, FL, FL, strand, c, s0, e1, NL, wseq, NL))
    units, qf = [], io.open(os.path.join(OUT, "map", aid + ".units.fa"), "w", encoding="utf8", newline=NL)
    for i in range(1, len(gen)):
        us, ue = gen[i - 1][1] + 1, gen[i][1]
        if ue <= us:
            continue
        u = wseq[us:ue + 1]
        units.append((i, us, ue, u, gen[i][0] - us))
        qf.write(">u%d%s%s%s" % (i, NL, u, NL))
    qf.close()
    if not units:
        continue
    ln = collections.Counter(len(u[3]) for u in units).most_common(1)[0][0]
    same = [u[3] for u in units if len(u[3]) == ln]
    canon = "".join(collections.Counter(x[i] for x in same).most_common(1)[0][0] for i in range(ln)) if same else units[0][3]
    cp = os.path.join(OUT, "canon_unit", aid + ".unit.fa")
    io.open(cp, "w", encoding="utf8", newline=NL).write(">%s_canon_unit len=%d%s%s%s" % (aid, ln, NL, canon, NL))
    r = subprocess.run("minimap2 -cx asm5 --cs -t 4 %s %s 2>/dev/null" % (cp, os.path.join(OUT, "map", aid + ".units.fa")),
                       shell=True, capture_output=True, text=True)
    best = {}
    for L in r.stdout.splitlines():
        f2 = L.split(TAB)
        if len(f2) < 12:
            continue
        q, rs, bs = f2[0], int(f2[7]), int(f2[13].split(":")[-1]) if f2[13].startswith("AS") else 0
        cs = ""
        for x in f2[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        if q not in best or bs > best[q][0]:
            best[q] = (bs, rs, cs, int(f2[9]), int(f2[10]))
    M = io.open(os.path.join(OUT, "map", aid + ".map.tsv"), "w", encoding="utf8", newline=NL)
    M.write(TAB.join(["species", "array_id", "chrom", "strand", "array_pos", "chrom_start", "chrom_end", "unit_len",
                      "gene_start_in_unit", "gene_len", "aln_len", "identity_pct", "n_var", "varlist"]) + NL)
    nv, ident, allv = [], 0, collections.Counter()
    for (i, us, ue, u, gs) in units:
        b = best.get("u%d" % i)
        v = cs_to_vars(b[2], b[1]) if b else []
        idp = 100.0 * b[3] / b[4] if b and b[4] else 0.0
        if strand == "-":
            cs_, ce_ = R1 - ue, R1 - us
        else:
            cs_, ce_ = L0 + us, L0 + ue
        M.write(TAB.join([SP, aid, c, strand, str(i), str(cs_), str(ce_), str(len(u)), str(gs + 1),
                          str(gen[i][1] - gen[i][0] + 1), str(b[4]) if b else "0", "%.2f" % idp, str(len(v)),
                          ",".join(v) if v else "-"]) + NL)
        nv.append(len(v))
        ident += (len(v) == 0)
        for x in v:
            allv[x] += 1
    M.close()
    nv.sort()
    SUM.write(TAB.join([SP, aid, c, str(s0), str(e1), strand, str(cop), str(len(units)), str(ln),
                        "%.1f" % (100.0 * ident / len(units)), str(len(allv)), str(nv[len(nv) // 2]),
                        os.path.join("canon", aid + ".fa"), os.path.join("map", aid + ".map.tsv")]) + NL)
    print("%s %s : 유닛 %d · 대표유닛 %d bp · 변이없는 유닛 %.1f%% · 변이자리 %d" % (SP, aid, len(units), ln, 100.0 * ident / len(units), len(allv)))
SUM.close()
