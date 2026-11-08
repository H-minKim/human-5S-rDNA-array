# -*- coding: utf-8 -*-
# 타 생물종 5S : 정본(flank-array-flank) + consensus map(유닛별 변이표) — 사람 작업과 같은 꼴 · 2판
#   1판의 결함 정정 : minimap2 가 못 맞춘 유닛을 "변이 0" 으로 세었다 -> 이제 blastn 으로 되받고, 그래도 안 되면 NOALN 으로 남긴다
#   유닛틀 : 사람 정본과 같게 "스페이서 다음 유전자" (유닛이 유전자로 끝난다) — 사람에 돌리면 2,240 bp · 유전자 2122 재현
#   유닛번호: array 의 유전자 방향으로 1 부터 + 유닛마다 염색체 좌표 병기 (사람의 '원위=1' 규칙은 다른 종에 근거 없음)
#   변이   : 1 차 minimap2 -cx asm20 --cs (사람 방법) · 2 차 blastn -task blastn (짧거나 많이 갈린 유닛)
#            토큰은 사람 지도와 같게 1484T(치환) · 437D1(결실) · 1056I3(삽입) · 좌표는 대표유닛 1-기반
#   인자 : 종키 fasta tblout struct표 출력폴더 [flank bp=20000]
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


def merge(tok):
    # 1 bp 씩 나온 결실/삽입을 사람 지도처럼 붙인다 : 100D1,101D1 -> 100D2
    out, i = [], 0
    while i < len(tok):
        t = tok[i]
        if "D" in t and t.split("D")[1].isdigit():
            p = int(t.split("D")[0])
            n = int(t.split("D")[1])
            while i + 1 < len(tok) and "D" in tok[i + 1] and tok[i + 1].split("D")[0].isdigit() and int(tok[i + 1].split("D")[0]) == p + n:
                n += int(tok[i + 1].split("D")[1])
                i += 1
            out.append("%dD%d" % (p, n))
        elif "I" in t and t.split("I")[1].isdigit():
            p = int(t.split("I")[0])
            n = int(t.split("I")[1])
            while i + 1 < len(tok) and "I" in tok[i + 1] and tok[i + 1].split("I")[0].isdigit() and int(tok[i + 1].split("I")[0]) == p:
                n += int(tok[i + 1].split("I")[1])
                i += 1
            out.append("%dI%d" % (p, n))
        else:
            out.append(t)
        i += 1
    return out


def cs_to_vars(cs, rs):
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
        else:
            i += 1
    return merge(out)


def btop_to_vars(btop, ss):
    # blast BTOP -> 같은 토큰. ss = 대표유닛(subject) 시작 (1-기반)
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
                out.append("%d%s" % (p + 1, q.upper()))
                p += 1
            elif q == "-":
                out.append("%dD1" % (p + 1))
                p += 1
            else:
                out.append("%dI1" % (p + 1))
    return merge(out)


SUM = io.open(os.path.join(OUT, "MAP_SUMMARY_%s.tsv" % SP), "w", encoding="utf8", newline=NL)
SUM.write(TAB.join(["species", "array_id", "chrom", "array_start", "array_end", "strand", "copies", "units_total",
                    "units_mm2", "units_blastn", "units_noaln", "canon_unit_len", "units_identical_pct_of_aligned",
                    "var_sites", "var_per_unit_med", "canon_fa", "map_tsv"]) + NL)
for (c, s0, e1, strand, cop) in ARR:
    aid = "%s__%s_%d_%d" % (SP, c, s0, e1)
    L0, R1 = max(1, s0 - FL), min(LEN[c], e1 + FL)
    seq = faidx("%s:%d-%d" % (c, L0, R1))
    if not seq:
        continue
    hits = [h for h in HIT[c] if s0 <= h[0] and h[1] <= e1]
    if strand == "-":
        wseq = rc(seq)
        gen = sorted([(len(seq) - (h[1] - L0) - 1, len(seq) - (h[0] - L0) - 1) for h in hits])
    else:
        wseq = seq
        gen = sorted([(h[0] - L0, h[1] - L0) for h in hits])
    io.open(os.path.join(OUT, "canon", aid + ".fa"), "w", encoding="utf8", newline=NL).write(
        ">%s flank%d_array_flank%d strand=%s chrom=%s array=%d-%d%s%s%s" % (aid, FL, FL, strand, c, s0, e1, NL, wseq, NL))
    units = []
    qp = os.path.join(OUT, "map", aid + ".units.fa")
    qf = io.open(qp, "w", encoding="utf8", newline=NL)
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
    best = {}
    r = subprocess.run("minimap2 -cx asm20 --cs -t 4 %s %s 2>/dev/null" % (cp, qp), shell=True, capture_output=True, text=True)
    for L in r.stdout.splitlines():
        f2 = L.split(TAB)
        if len(f2) < 12:
            continue
        q, rs = f2[0], int(f2[7])
        bs = 0
        cs = ""
        for x in f2[12:]:
            if x.startswith("AS:i:"):
                bs = int(x[5:])
            elif x.startswith("cs:Z:"):
                cs = x[5:]
        if q not in best or bs > best[q][0]:
            best[q] = (bs, "mm2", cs_to_vars(cs, rs), int(f2[9]), int(f2[10]))
    # 2 차 : minimap2 가 못 맞춘 유닛만 blastn
    miss = [u for u in units if ("u%d" % u[0]) not in best]
    if miss:
        mp = os.path.join(OUT, "map", aid + ".miss.fa")
        with io.open(mp, "w", encoding="utf8", newline=NL) as g:
            for (i, us, ue, u, gs) in miss:
                g.write(">u%d%s%s%s" % (i, NL, u, NL))
        rb = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1 -query %s -subject %s "
                            "-outfmt '6 qseqid sstart send pident length btop bitscore' 2>/dev/null" % (mp, cp),
                            shell=True, capture_output=True, text=True)
        for L in rb.stdout.splitlines():
            f3 = L.split(TAB)
            if len(f3) < 7:
                continue
            q, ss, se, pid, aln, btop, bs = f3[0], int(f3[1]), int(f3[2]), float(f3[3]), int(f3[4]), f3[5], float(f3[6])
            if q in best and bs <= best[q][0]:
                continue
            if ss > se:
                best[q] = (bs, "blastn_rev", [], aln, aln)
            else:
                best[q] = (bs, "blastn", btop_to_vars(btop, ss), aln, aln)
        try:
            os.remove(mp)
        except OSError:
            pass
    M = io.open(os.path.join(OUT, "map", aid + ".map.tsv"), "w", encoding="utf8", newline=NL)
    M.write(TAB.join(["species", "array_id", "chrom", "strand", "array_pos", "chrom_start", "chrom_end", "unit_len",
                      "gene_start_in_unit", "gene_len", "aln_method", "aln_len", "identity_pct", "n_var", "varlist"]) + NL)
    nv, ident, allv = [], 0, collections.Counter()
    cnt = collections.Counter()
    for (i, us, ue, u, gs) in units:
        b = best.get("u%d" % i)
        if b:
            meth, v, alen = b[1], b[2], b[4]
            cnt[meth.split("_")[0]] += 1
            idps = "%.2f" % (100.0 * b[3] / b[4]) if b[4] else "NA"
            if meth == "blastn_rev":
                vstr, nvstr = "REV", "NA"
            else:
                vstr = ",".join(v) if v else "-"
                nvstr = str(len(v))
                nv.append(len(v))
                ident += (len(v) == 0)
                for x in v:
                    allv[x] += 1
        else:
            meth, vstr, nvstr, idps, alen = "none", "NOALN", "NA", "NA", 0
            cnt["none"] += 1
        cs_, ce_ = (R1 - ue, R1 - us) if strand == "-" else (L0 + us, L0 + ue)
        M.write(TAB.join([SP, aid, c, strand, str(i), str(cs_), str(ce_), str(len(u)), str(gs + 1),
                          str(gen[i][1] - gen[i][0] + 1), meth, str(alen), idps, nvstr, vstr]) + NL)
    M.close()
    nv.sort()
    na = len(nv)
    SUM.write(TAB.join([SP, aid, c, str(s0), str(e1), strand, str(cop), str(len(units)), str(cnt["mm2"]),
                        str(cnt["blastn"]), str(cnt["none"]), str(ln),
                        ("%.1f" % (100.0 * ident / na)) if na else "NA", str(len(allv)),
                        str(nv[na // 2]) if na else "NA", os.path.join("canon", aid + ".fa"),
                        os.path.join("map", aid + ".map.tsv")]) + NL)
    print("%s : 유닛 %d (mm2 %d · blastn %d · 실패 %d) · 대표유닛 %d bp · 변이없음 %s · 변이자리 %d"
          % (aid, len(units), cnt["mm2"], cnt["blastn"], cnt["none"], ln,
             ("%.1f%%" % (100.0 * ident / na)) if na else "NA", len(allv)))
SUM.close()
