# -*- coding: utf-8 -*-
# 타 생물종 5S : array 마다 참조(대표) 유닛을 만든다 — 사람 G8v4 방법을 종별로 환산
#   틀 : 유닛 = 유전자 끝 +1 .. 다음 유전자 끝  (유닛이 [스페이서][유전자] 로 끝난다 · 사람 정본유닛과 같은 위상)
#   유전자 경계 : 그 array 의 유전자 합의로 blastn 한 뒤 질의 좌표로 전체 길이까지 확장
#                 (+ 가닥 : 시작 = sstart-(qstart-1) · 끝 = send+(L-qend) / - 가닥은 반대)
#                 사람 실측 : HSP 는 유전자 끝 염기가 변이면 한두 염기 앞에서 멈춘다
#   종별 환산 : PAD = 5% (최소 20) · CTX = 9% (최소 40) · LO/HI = 주기의 -6% / +7% (최소 ±30bp)
#   대표 유닛 : 최빈 길이 유닛들의 자리별 다수결 → 유전자가 유닛 끝에 오는지 서열로 확인
#   산출 : refunit/<array_id>.unit.fa · refunit/<array_id>.gene.fa · REFUNIT_PARAMS.tsv
#   인자 : 종키 fasta tblout struct표 출력폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
SP, FA, TBL, STRUCT, OUT = sys.argv[1:6]
os.makedirs(os.path.join(OUT, "refunit"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(COMP.get(x, "N") for x in reversed(s))


def faidx(reg):
    r = subprocess.run(["samtools", "faidx", FA, reg], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


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
    if r[ix["species"]] == SP:
        ARR.append((r[ix["chrom"]], int(r[ix["array_start"]]), int(r[ix["array_end"]]), r[ix["strand_major"]][0]))
P = io.open(os.path.join(OUT, "REFUNIT_PARAMS_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
P.write(TAB.join(["species", "array_id", "chrom", "array_start", "array_end", "strand", "n_gene_hsp", "period_med",
                  "LO", "HI", "PAD", "CTX", "n_units_in_range", "n_break", "refunit_len", "gene_len",
                  "gene_start_in_refunit", "gene_at_end", "gene_identity_pct", "units_ge99_pct_blastn",
                  "units_median_identity", "refunit_usable", "note"]) + NL)
for (c, s0, e1, strand) in ARR:
    aid = "%s__%s_%d_%d" % (SP, c, s0, e1)
    hits = [h for h in HIT[c] if s0 <= h[0] and h[1] <= e1]
    if len(hits) < 3:
        continue
    mar = 3000
    L0, R1 = max(1, s0 - mar), min(LEN[c], e1 + mar)
    seq = faidx("%s:%d-%d" % (c, L0, R1))
    if not seq:
        continue
    # 1) 그 array 의 유전자 합의 (cmsearch hit 서열의 자리별 다수결 · 최빈 길이만)
    gs = []
    for (a, b, st) in hits:
        g = seq[a - L0:b - L0 + 1]
        gs.append(rc(g) if st == "-" else g)
    gl = collections.Counter(len(x) for x in gs).most_common(1)[0][0]
    same = [x for x in gs if len(x) == gl]
    gene = "".join(collections.Counter(x[i] for x in same).most_common(1)[0][0] for i in range(gl))
    gp = os.path.join(OUT, "refunit", aid + ".gene.fa")
    io.open(gp, "w", encoding="utf8", newline=NL).write(">%s_gene len=%d%s%s%s" % (aid, gl, NL, gene, NL))
    # 2) blastn 으로 유전자 HSP -> 질의 좌표로 전체 길이까지 확장
    rp = os.path.join(OUT, "refunit", aid + ".region.fa")
    io.open(rp, "w", encoding="utf8", newline=NL).write(">region%s%s%s" % (NL, seq, NL))
    rb = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1e-5 -query %s -subject %s "
                        "-outfmt '6 qstart qend sstart send pident bitscore' 2>/dev/null" % (gp, rp),
                        shell=True, capture_output=True, text=True)
    G = []
    for L in rb.stdout.splitlines():
        q1, q2, s1, s2, pid, bs = L.split(TAB)
        q1, q2, s1, s2 = int(q1), int(q2), int(s1), int(s2)
        if s1 <= s2:
            st, a, b = "+", s1 - (q1 - 1), s2 + (gl - q2)
        else:
            st, a, b = "-", s2 - (gl - q2), s1 + (q1 - 1)
        G.append((a, b, st, float(bs)))
    G.sort()
    keep, last = [], -1
    for (a, b, st, bs) in G:
        if a > last:
            keep.append((a, b, st))
            last = b
    major = collections.Counter(x[2] for x in keep).most_common(1)[0][0] if keep else strand
    if major == "-":
        n = len(seq)
        wseq = rc(seq)
        gg = sorted([(n - b + 1, n - a + 1) for (a, b, st) in keep if st == major])
    else:
        wseq = seq
        gg = sorted([(a, b) for (a, b, st) in keep if st == major])
    if len(gg) < 3:
        P.write(TAB.join([SP, aid, c, str(s0), str(e1), strand, str(len(keep)), "NA", "NA", "NA", "NA", "NA",
                          "0", "NA", "NA", str(gl), "NA", "NA", "NA", "NA", "NA", "N", "유전자 HSP 부족"]) + NL)
        continue
    per = [gg[i + 1][1] - gg[i][1] for i in range(len(gg) - 1)]
    per.sort()
    pm = per[len(per) // 2]
    # 사람 기준(주기 2240 에 LO 2100 · HI 2400 = -6.25% / +7.1%)을 종별로 환산하되 더 너그러운 쪽을 쓴다
    LO = min(pm - 30, int(round(pm * 0.94)))
    HI = max(pm + 30, int(round(pm * 1.07)))
    PAD = max(20, int(round(pm * 0.05)))
    CTX = max(40, int(round(pm * 0.09)))
    units, brk = [], 0
    for i in range(len(gg) - 1):
        a, b = gg[i][1] + 1, gg[i + 1][1]
        n = b - a + 1
        if LO <= n <= HI:
            units.append(wseq[a - 1:b])
        else:
            brk += 1
    if not units:
        P.write(TAB.join([SP, aid, c, str(s0), str(e1), major, str(len(gg)), str(pm), str(LO), str(HI), str(PAD),
                          str(CTX), "0", str(brk), "NA", str(gl), "NA", "NA", "NA", "NA", "NA", "N", "범위 안 유닛 없음"]) + NL)
        continue
    ml = collections.Counter(len(u) for u in units).most_common(1)[0][0]
    pool = [u for u in units if len(u) == ml]
    ref = "".join(collections.Counter(u[i] for u in pool).most_common(1)[0][0] for i in range(ml))
    up = os.path.join(OUT, "refunit", aid + ".unit.fa")
    io.open(up, "w", encoding="utf8", newline=NL).write(">%s_refunit len=%d gene_len=%d%s%s%s" % (aid, ml, gl, NL, ref, NL))
    # 3) 검증 : 유전자가 유닛 끝에 오는가 (서열로) · 유닛들이 대표와 얼마나 같은가
    tail = ref[-gl:]
    idp = 100.0 * sum(1 for i in range(gl) if tail[i] == gene[i]) / gl
    rb2 = subprocess.run("blastn -task blastn -word_size 11 -dust no -query %s -subject %s "
                         "-outfmt '6 sstart send pident bitscore' 2>/dev/null" % (gp, up),
                         shell=True, capture_output=True, text=True)
    gstart = "NA"
    bestbs = -1
    for L in rb2.stdout.splitlines():
        s1, s2, pid, bs = L.split(TAB)
        if float(bs) > bestbs:
            bestbs = float(bs)
            gstart = min(int(s1), int(s2))
    atend = "Y" if gstart != "NA" and abs((ml - gl + 1) - gstart) <= 2 else "N"
    # 동질성은 정렬(blastn)로 잰다 — 자리별 직접 비교는 indel 하나에 무너진다
    ufp = os.path.join(OUT, "refunit", aid + ".units.fa")
    with io.open(ufp, "w", encoding="utf8", newline=NL) as g2:
        for k, u in enumerate(units):
            g2.write(">u%d%s%s%s" % (k + 1, NL, u, NL))
    rb3 = subprocess.run("blastn -task blastn -word_size 11 -dust no -query %s -subject %s "
                         "-outfmt '6 qseqid pident length qlen bitscore' 2>/dev/null" % (ufp, up),
                         shell=True, capture_output=True, text=True)
    bb = {}
    for L in rb3.stdout.splitlines():
        q, pid, aln, qlen, bs = L.split(TAB)
        if int(aln) < 0.5 * int(qlen):
            continue
        if q not in bb or float(bs) > bb[q][0]:
            bb[q] = (float(bs), float(pid))
    pids = sorted(v[1] for v in bb.values())
    medid = "%.2f" % pids[len(pids) // 2] if pids else "NA"
    ge = sum(1 for v in bb.values() if v[1] >= 99.0)
    gepct = "%.1f" % (100.0 * ge / len(units)) if units else "NA"
    try:
        os.remove(ufp)
    except OSError:
        pass
    # 쓸 수 있는 대표인가 : 유닛 5 개 이상 + 유전자가 끝 + 유닛 절반 이상이 대표에 정렬
    why = []
    if len(units) < 5:
        why.append("유닛<5")
    if atend != "Y":
        why.append("유전자가 끝 아님")
    if len(bb) < 0.5 * len(units):
        why.append("정렬된 유닛<50%")
    usable = "Y" if not why else "N"
    P.write(TAB.join([SP, aid, c, str(s0), str(e1), major, str(len(gg)), str(pm), str(LO), str(HI), str(PAD), str(CTX),
                      str(len(units)), str(brk), str(ml), str(gl), str(gstart), atend, "%.1f" % idp,
                      gepct, medid, usable, ";".join(why)]) + NL)
    print("%s : 유전자HSP %d · 주기 %d (LO %d HI %d) · 유닛 %d (끊김 %d) · 대표 %d bp · 유전자 끝 %s (자리 %s · 일치 %.1f%%)"
          % (aid, len(gg), pm, LO, HI, len(units), brk, ml, atend, gstart, idp))
    for x in (rp,):
        try:
            os.remove(x)
        except OSError:
            pass
P.close()
