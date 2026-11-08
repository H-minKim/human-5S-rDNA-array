# -*- coding: utf-8 -*-
# 타 생물종 5S : G8v4 방식 유닛 지도 (array 마다 제 참조 유닛으로)
#   사람 consmapG8v4.py 의 규칙을 그대로 옮기고, 종별 값(참조유닛·PAD·CTX·LO/HI)만 REFUNIT_PARAMS.tsv 에서 읽는다
#     · 유닛 경계 = 유전자 끝 +1 .. 다음 유전자 끝 (유전자 HSP 를 질의 좌표로 전체 길이까지 확장)
#     · 자를 때 앞뒤 PAD 여유 · 정렬 대상은 [참조유닛 끝 CTX][참조유닛][참조유닛 앞 CTX]
#     · 채움 금지 — 변이는 정렬(cs/BTOP)로만. 가운데 유닛 구간을 다 못 덮으면 NO_FULL
#     · 간격이 LO~HI 밖이면 유닛을 만들지 않고 끊김(n_break)
#     · 양끝 조각은 PARTIAL 로 적고 변이를 부르지 않는다
#   산출 : g8map/<array_id>.g8map.tsv · G8MAP_SUM.tsv
#   인자 : 종키 fasta tblout REFUNIT_PARAMS.tsv 참조유닛폴더 출력폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
SP, FA, TBL, PARAM, RDIR, OUT = sys.argv[1:7]
os.makedirs(os.path.join(OUT, "g8map"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))


def rdfa(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()


def faidx(reg):
    r = subprocess.run(["samtools", "faidx", FA, reg], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


LEN = dict((L.split(TAB)[0], int(L.split(TAB)[1])) for L in io.open(FA + ".fai", encoding="utf8"))
PR = [r for r in csv.DictReader(io.open(PARAM, encoding="utf8"), delimiter=TAB)
      if r["species"] == SP and r["refunit_usable"] == "Y"]
SUM = io.open(os.path.join(OUT, "G8MAP_SUM_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
SUM.write(TAB.join(["species", "array_id", "chrom", "strand", "refunit_len", "gene_len", "N_units", "n_break",
                    "verify_OK", "verify_NO_FULL", "verify_NOALN", "n_var_total", "n_var_gene_total",
                    "var_sites", "median_var_per_unit"]) + NL)


def cs_vars(cs, rs):
    out, p, i = [], rs, 0
    while i < len(cs):
        ch = cs[i]
        if ch == ":":
            j = i + 1
            while j < len(cs) and cs[j].isdigit():
                j += 1
            p += int(cs[i + 1:j])
            i = j
        elif ch == "*":
            out.append((p + 1, "S", cs[i + 2].upper()))
            p += 1
            i += 3
        elif ch == "-":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            out.append((p + 1, "D", j - i - 1))
            p += j - i - 1
            i = j
        elif ch == "+":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            out.append((p + 1, "I", j - i - 1))
            i = j
        else:
            i += 1
    return out


def btop_vars(btop, ss):
    out, p, i = [], ss - 1, 0
    while i < len(btop):
        if btop[i].isdigit():
            j = i
            while j < len(btop) and btop[j].isdigit():
                j += 1
            p += int(btop[i:j])
            i = j
        else:
            q, s = btop[i], btop[i + 1]
            i += 2
            if q != "-" and s != "-":
                out.append((p + 1, "S", q.upper()))
                p += 1
            elif q == "-":
                out.append((p + 1, "D", 1))
                p += 1
            else:
                out.append((p + 1, "I", 1))
    return out


def fmt(tok):
    # 사람 지도와 같은 표기 : 1096G · 1476D10 · 1275I2 (연속 결실/삽입은 합친다)
    out, i = [], 0
    while i < len(tok):
        p, k, v = tok[i]
        if k == "S":
            out.append("%d%s" % (p, v))
            i += 1
        elif k == "D":
            n = v
            while i + 1 < len(tok) and tok[i + 1][1] == "D" and tok[i + 1][0] == p + n:
                n += tok[i + 1][2]
                i += 1
            out.append("%dD%d" % (p, n))
            i += 1
        else:
            n = v
            while i + 1 < len(tok) and tok[i + 1][1] == "I" and tok[i + 1][0] == p:
                n += tok[i + 1][2]
                i += 1
            out.append("%dI%d" % (p, n))
            i += 1
    return out


for r in PR:
    aid = r["array_id"]
    c, s0, e1 = r["chrom"], int(r["array_start"]), int(r["array_end"])
    strand = r["strand"]
    U, GL = int(r["refunit_len"]), int(r["gene_len"])
    PAD, CTX = int(r["PAD"]), int(r["CTX"])
    LO, HI = int(r["LO"]), int(r["HI"])
    ref = rdfa(os.path.join(RDIR, aid + ".unit.fa"))
    gene = rdfa(os.path.join(RDIR, aid + ".gene.fa"))
    if len(ref) != U:
        continue
    ctxp = os.path.join(OUT, "g8map", aid + ".ctx.fa")
    io.open(ctxp, "w", encoding="utf8", newline=NL).write(">ctx%s%s%s" % (NL, ref[-CTX:] + ref + ref[:CTX], NL))
    gp = os.path.join(RDIR, aid + ".gene.fa")
    mar = 3 * U + 3000
    L0, R1 = max(1, s0 - mar), min(LEN[c], e1 + mar)
    seq = faidx("%s:%d-%d" % (c, L0, R1))
    rp = os.path.join(OUT, "g8map", aid + ".region.fa")
    io.open(rp, "w", encoding="utf8", newline=NL).write(">region%s%s%s" % (NL, seq, NL))
    rb = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1e-5 -query %s -subject %s "
                        "-outfmt '6 qstart qend sstart send bitscore' 2>/dev/null" % (gp, rp),
                        shell=True, capture_output=True, text=True)
    G = []
    for L in rb.stdout.splitlines():
        q1, q2, s1, s2, bs = L.split(TAB)
        q1, q2, s1, s2 = int(q1), int(q2), int(s1), int(s2)
        if s1 <= s2:
            G.append((s1 - (q1 - 1), s2 + (GL - q2), "+"))
        else:
            G.append((s2 - (GL - q2), s1 + (q1 - 1), "-"))
    G.sort()
    keep, last = [], -1
    for (a, b, st) in G:
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
        continue
    # 유닛 자르기 (PAD 여유) + 끊김
    qs, meta, brk = [], [], 0
    pos = 0
    for i in range(len(gg) - 1):
        a, b = gg[i][1] + 1, gg[i + 1][1]
        if not (LO <= b - a + 1 <= HI):
            brk += 1
            continue
        pos += 1
        la, lb = max(1, a - PAD), min(len(wseq), b + PAD)
        meta.append((pos, a, b, la, lb, b - a + 1))
        qs.append(">u%d%s%s%s" % (pos, NL, wseq[la - 1:lb], NL))
    if not meta:
        continue
    qp = os.path.join(OUT, "g8map", aid + ".units.fa")
    io.open(qp, "w", encoding="utf8", newline=NL).write("".join(qs))
    best = {}
    rm = subprocess.run("minimap2 -cx asm20 --cs -t 4 %s %s 2>/dev/null" % (ctxp, qp), shell=True, capture_output=True, text=True)
    for L in rm.stdout.splitlines():
        f = L.split(TAB)
        if len(f) < 12:
            continue
        q, rs, re_ = f[0], int(f[7]), int(f[8])
        bs, cs = 0, ""
        for x in f[12:]:
            if x.startswith("AS:i:"):
                bs = int(x[5:])
            elif x.startswith("cs:Z:"):
                cs = x[5:]
        if q not in best or bs > best[q][0]:
            best[q] = (bs, "mm2", cs_vars(cs, rs), rs + 1, re_)
    miss = [m for m in meta if ("u%d" % m[0]) not in best]
    if miss:
        mp = os.path.join(OUT, "g8map", aid + ".miss.fa")
        with io.open(mp, "w", encoding="utf8", newline=NL) as g2:
            for m in miss:
                g2.write(">u%d%s%s%s" % (m[0], NL, wseq[m[3] - 1:m[4]], NL))
        rbb = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1 -query %s -subject %s "
                             "-outfmt '6 qseqid sstart send btop bitscore' 2>/dev/null" % (mp, ctxp),
                             shell=True, capture_output=True, text=True)
        for L in rbb.stdout.splitlines():
            f = L.split(TAB)
            if len(f) < 5:
                continue
            q, ss, se, btop, bs = f[0], int(f[1]), int(f[2]), f[3], float(f[4])
            if ss > se:
                continue
            if q not in best or bs > best[q][0]:
                best[q] = (bs, "blastn", btop_vars(btop, ss), ss, se)
        try:
            os.remove(mp)
        except OSError:
            pass
    M = io.open(os.path.join(OUT, "g8map", aid + ".g8map.tsv"), "w", encoding="utf8", newline=NL)
    M.write(TAB.join(["species", "array_id", "chrom", "strand", "array_pos", "N", "chrom_start", "chrom_end",
                      "unit_len", "gap_to_next", "aln_method", "cov_start", "cov_end", "n_var", "n_var_gene",
                      "varlist", "verify"]) + NL)
    N = len(meta)
    cnt = collections.Counter()
    nv, nvg, allv = [], 0, collections.Counter()
    for k, (p, a, b, la, lb, ulen) in enumerate(meta):
        bb = best.get("u%d" % p)
        gap = (meta[k + 1][1] - b - 1) if k + 1 < len(meta) else -1
        cs_, ce_ = (R1 - (lb - 0), R1 - (la - 0)) if major == "-" else (L0 + a - 1, L0 + b - 1)
        if major == "-":
            cs_, ce_ = R1 - b + 1, R1 - a + 1
        if not bb:
            M.write(TAB.join([SP, aid, c, major, str(p), str(N), str(cs_), str(ce_), str(ulen), str(gap),
                              "none", "NA", "NA", "NA", "NA", "NOALN", "NOALN"]) + NL)
            cnt["NOALN"] += 1
            continue
        # ctx 좌표 -> 참조유닛 좌표 (가운데 유닛만)
        toks = [(pp - CTX, kk, vv) for (pp, kk, vv) in bb[2] if CTX < pp <= CTX + U]
        cov0, cov1 = bb[3] - CTX, bb[4] - CTX
        full = (cov0 <= 1 and cov1 >= U)
        ver = "OK" if full else "NO_FULL"
        tl = fmt(toks)
        ng = sum(1 for (pp, kk, vv) in toks if pp > U - GL)
        M.write(TAB.join([SP, aid, c, major, str(p), str(N), str(cs_), str(ce_), str(ulen), str(gap), bb[1],
                          str(max(1, cov0)), str(min(U, cov1)), str(len(tl)), str(ng),
                          ",".join(tl) if tl else "-", ver]) + NL)
        cnt[ver] += 1
        nv.append(len(tl))
        nvg += ng
        for x in tl:
            allv[x] += 1
    M.close()
    nv.sort()
    SUM.write(TAB.join([SP, aid, c, major, str(U), str(GL), str(N), str(brk), str(cnt["OK"]), str(cnt["NO_FULL"]),
                        str(cnt["NOALN"]), str(sum(nv)), str(nvg), str(len(allv)),
                        str(nv[len(nv) // 2]) if nv else "NA"]) + NL)
    print("%s : 유닛 %d (끊김 %d) · OK %d · NO_FULL %d · NOALN %d · 변이자리 %d · 유전자부 변이 %d"
          % (aid, N, brk, cnt["OK"], cnt["NO_FULL"], cnt["NOALN"], len(allv), nvg))
    for x in (rp, qp, ctxp):
        try:
            os.remove(x)
        except OSError:
            pass
SUM.close()
