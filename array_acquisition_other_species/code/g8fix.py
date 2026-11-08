# -*- coding: utf-8 -*-
# G8v4 규칙 유닛 지도를 정본 fasta 만으로 만든다 (게놈 파일 불필요) + 유전자 k 개짜리 유닛 지원
#   species_g8map.py 의 규칙을 그대로 옮긴다 :
#     · 유닛 경계 = 유전자 끝 +1 .. k 개 뒤 유전자 끝 (유전자 HSP 를 질의 좌표로 전체 길이까지 확장)
#     · 자를 때 앞뒤 PAD 여유 · 정렬 대상은 [참조유닛 끝 CTX][참조유닛][참조유닛 앞 CTX]
#     · 채움 금지 — 변이는 정렬(cs/BTOP)로만. 가운데 유닛을 다 못 덮으면 NO_FULL
#     · 간격이 LO~HI 밖이면 유닛을 만들지 않고 끊김(n_break)
#   달라진 것 : 참조 유닛을 최빈 길이가 아니라 중앙값 길이로 고른다 (퍼진 array 서 임의 선택 방지)
#   인자 : 정본fa cms_tbl k 출력폴더
#   산출 : <out>/<array_id>.g8map.tsv · .unit.fa · .gene.fa · PARAMS.tsv (덧붙임)
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import statistics
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
CANFA, TBL, K, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
# 온전한 유전자만 경계로 쓴다 (5S 3' 조각이 섞인 array 에서 주기가 반토막 난다 — 애기장대 chr4 15.92 Mb 실측)
MINCOV = float(sys.argv[5]) if len(sys.argv) > 5 else 0.80
# 참조 유닛 고르는 법 : consensus(같은 길이 유닛 자리별 다수결 · 기본) / medoid(가장 중심에 있는 실제 유닛)
#   균질하지 않은 array 에서는 consensus 가 한두 유닛만 모아 참조만 이질적인 사태가 난다
#   (붉은털원숭이 20.4 Mb 실측 : 26 유닛 전부가 같은 353 자리에서 참조와 달랐다)
REFMODE = sys.argv[6] if len(sys.argv) > 6 else 'consensus'
os.makedirs(OUT, exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def rc(s):
    return "".join(COMP.get(x, "N") for x in reversed(s))


def merge(tok):
    out, i = [], 0
    while i < len(tok):
        p, kk, vv = tok[i]
        if kk == "D":
            n = vv
            while i + 1 < len(tok) and tok[i + 1][1] == "D" and tok[i + 1][0] == p + n:
                n += tok[i + 1][2]
                i += 1
            out.append((p, "D", n))
        elif kk == "I":
            n = vv
            while i + 1 < len(tok) and tok[i + 1][1] == "I" and tok[i + 1][0] == p:
                n += tok[i + 1][2]
                i += 1
            out.append((p, "I", n))
        else:
            out.append((p, kk, vv))
        i += 1
    return out


def fmt(tok):
    o = []
    for (p, kk, vv) in tok:
        o.append("%d%s" % (p, vv) if kk == "S" else "%d%s%d" % (p, kk, vv))
    return o


def cs_vars(cs, rs):
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
            out.append((p + 1, "S", cs[i + 2].upper()))
            p += 1
            i += 3
        elif c == "-":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            n = j - i - 1
            out.append((p + 1, "D", n))
            p += n
            i = j
        elif c == "+":
            j = i + 1
            while j < len(cs) and cs[j].isalpha():
                j += 1
            out.append((p + 1, "I", j - i - 1))
            i = j
        else:
            i += 1
    return merge(out)


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
    return merge(out)


aid = os.path.basename(CANFA).replace(".fa", "")
SP = aid.split("__")[0]
h = io.open(CANFA, encoding="utf8")
hd = h.readline().strip()
wseq = "".join(x.strip() for x in h).upper()
kv = dict(x.split("=", 1) for x in hd.split() if "=" in x)
a0, a1 = [int(x) for x in kv["array"].split("-")]
chrom, strand = kv["chrom"], kv["strand"]
span = a1 - a0 + 1
FL = (len(wseq) - span) // 2


def to_chrom(p):
    # 정본창 1-기반 p -> 염색체 좌표
    return (a1 + FL) - p + 1 if strand == "-" else (a0 - FL) + p - 1


def to_win(cpos):
    return (a1 + FL) - cpos + 1 if strand == "-" else cpos - (a0 - FL) + 1


# 1) 그 array 의 유전자 합의 (cmsearch hit · 정본창 방향으로)
gs = []
for L in io.open(TBL, encoding="utf8"):
    if L.startswith("#"):
        continue
    ff = L.split()
    if len(ff) < 10 or ff[0] != chrom:
        continue
    x, y = int(ff[7]), int(ff[8])
    lo, hi = min(x, y), max(x, y)
    if not (a0 <= lo and hi <= a1):
        continue
    w1, w2 = sorted([to_win(lo), to_win(hi)])
    g = wseq[w1 - 1:w2]
    gs.append(g if ff[9] == strand else rc(g))
if len(gs) < 3:
    sys.exit("유전자 3 개 미만 : " + aid)
gl = collections.Counter(len(x) for x in gs).most_common(1)[0][0]
same = [x for x in gs if len(x) == gl]
gene = "".join(collections.Counter(x[i] for x in same).most_common(1)[0][0] for i in range(gl))
gp = os.path.join(OUT, aid + ".gene.fa")
io.open(gp, "w", encoding="utf8", newline=NL).write(">%s_gene len=%d%s%s%s" % (aid, gl, NL, gene, NL))
# 2) 유전자 HSP -> 질의 좌표로 전체 길이 확장 (정본창 안, array ± 유닛 1 개까지만)
rp = os.path.join(OUT, aid + ".region.fa")
io.open(rp, "w", encoding="utf8", newline=NL).write(">region%s%s%s" % (NL, wseq, NL))
rb = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1e-5 -query %s -subject %s "
                    "-outfmt '6 qstart qend sstart send bitscore length' 2>/dev/null" % (gp, rp),
                    shell=True, capture_output=True, text=True)
G = []
ndrop = 0
for L in rb.stdout.splitlines():
    q1, q2, s1, s2, bs, aln = L.split(TAB)
    q1, q2, s1, s2, aln = int(q1), int(q2), int(s1), int(s2), int(aln)
    if aln < MINCOV * gl:
        ndrop += 1
        continue
    if s1 <= s2:
        G.append((s1 - (q1 - 1), s2 + (gl - q2), "+"))
    else:
        G.append((s2 - (gl - q2), s1 + (q1 - 1), "-"))
G.sort()
keep, last = [], -1
for (a, b, st) in G:
    if a > last:
        keep.append((a, b, st))
        last = b
major = collections.Counter(x[2] for x in keep).most_common(1)[0][0] if keep else "+"
print("   유전자 HSP %d 개 · 조각(덮힘<%.2f) %d 개 버림" % (len(G), MINCOV, ndrop))
gg = sorted([(a, b) for (a, b, st) in keep if st == major])
if len(gg) < 3:
    sys.exit("유전자 HSP 부족 : " + aid)
# 3) k 개 단위 주기와 파라미터
per = sorted([gg[i + K][1] - gg[i][1] for i in range(len(gg) - K)])
pm = per[len(per) // 2]
LO = min(pm - 30, int(round(pm * 0.94)))
HI = max(pm + 30, int(round(pm * 1.07)))
PAD = max(20, int(round(pm * 0.05)))
CTX = max(40, int(round(pm * 0.09)))
# 위상 맞춤 걸음 : k 개 창이 범위 밖이면 k 칸이 아니라 1 칸만 나아간다
#   (이중 유닛 array 는 곳곳에서 위상이 밀린다 — 붉은털원숭이 20.4 Mb 는 AB 교대 중 2 곳이 AA·BB.
#    k 칸씩 건너뛰면 그 뒤가 전부 반대 위상으로 잡혀 절반이 NO_FULL 이 된다)
def walk():
    cuts, nbrk, i = [], 0, 0
    while i + K < len(gg):
        a, b = gg[i][1] + 1, gg[i + K][1]
        if LO <= b - a + 1 <= HI:
            cuts.append((a, b))
            i += K
        else:
            nbrk += 1
            i += 1
    return cuts, nbrk


cuts0, brk0 = walk()
raw = [wseq[a - 1:b] for (a, b) in cuts0]
if not raw:
    sys.exit("범위 안 유닛 없음 : " + aid)
lens = sorted(len(u) for u in raw)
med = statistics.median(lens)
tgt = min(lens, key=lambda x: (abs(x - med), x))
pool = [u for u in raw if len(u) == tgt]
if len(pool) < 3:
    band = [u for u in raw if abs(len(u) - med) <= max(20, 0.02 * med)]
    if band:
        tgt = collections.Counter(len(u) for u in band).most_common(1)[0][0]
        pool = [u for u in band if len(u) == tgt]
if REFMODE == "medoid":
    # 가장 중심에 있는 유닛 : 다른 유닛들과의 blastn 비트점수 합이 가장 큰 것
    cand = [u for u in raw if abs(len(u) - med) <= max(30, 0.01 * med)]
    if len(cand) < 2:
        cand = raw
    cp0 = os.path.join(OUT, aid + ".cand.fa")
    with io.open(cp0, "w", encoding="utf8", newline=NL) as g0:
        for i2, u in enumerate(cand):
            g0.write(">c%d%s%s%s" % (i2, NL, u, NL))
    rr = subprocess.run("blastn -task blastn -word_size 16 -dust no -evalue 1e-10 -query %s -subject %s "
                        "-outfmt '6 qseqid sseqid bitscore pident length' 2>/dev/null" % (cp0, cp0),
                        shell=True, capture_output=True, text=True)
    sc = collections.Counter()
    pid = collections.defaultdict(list)
    for L5 in rr.stdout.splitlines():
        q, sj, bs, pp, aln = L5.split(TAB)
        if q == sj:
            continue
        sc[q] += float(bs)
        pid[q].append(float(pp))
    if sc:
        bestq = max(sc, key=lambda k2: sc[k2])
        ref = cand[int(bestq[1:])]
        medp = sorted(pid[bestq])
        print("   medoid 참조 : %s · 다른 유닛과 일치도 중앙 %.2f%%" % (bestq, medp[len(medp) // 2]))
    else:
        ref = cand[0]
    try:
        os.remove(cp0)
    except OSError:
        pass
else:
    ref = "".join(collections.Counter(u[i] for u in pool).most_common(1)[0][0] for i in range(tgt))
U = len(ref)
up = os.path.join(OUT, aid + ".unit.fa")
io.open(up, "w", encoding="utf8", newline=NL).write(
    ">%s_refunit len=%d gene_len=%d genes_per_unit=%d source=%s%s%s%s" % (aid, U, gl, K, REFMODE, NL, ref, NL))
ctxp = os.path.join(OUT, aid + ".ctx.fa")
io.open(ctxp, "w", encoding="utf8", newline=NL).write(">ctx%s%s%s" % (NL, ref[-CTX:] + ref + ref[:CTX], NL))
# 4) 유닛 자르기(PAD) + 끊김
meta, qs = [], []
cuts, brk = walk()
for pos, (a, b) in enumerate(cuts, 1):
    la, lb = max(1, a - PAD), min(len(wseq), b + PAD)
    meta.append((pos, a, b, la, lb, b - a + 1))
    qs.append(">u%d%s%s%s" % (pos, NL, wseq[la - 1:lb], NL))
qp = os.path.join(OUT, aid + ".units.fa")
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
    mp = os.path.join(OUT, aid + ".miss.fa")
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
M = io.open(os.path.join(OUT, aid + ".g8map.tsv"), "w", encoding="utf8", newline=NL)
M.write(TAB.join(["species", "array_id", "chrom", "strand", "array_pos", "N", "chrom_start", "chrom_end",
                  "unit_len", "gap_to_next", "aln_method", "cov_start", "cov_end", "n_var", "n_var_gene",
                  "varlist", "verify"]) + NL)
N = len(meta)
cnt = collections.Counter()
nv, nvg, allv = [], 0, collections.Counter()
for kk, (p, a, b, la, lb, ulen) in enumerate(meta):
    bb = best.get("u%d" % p)
    gap = (meta[kk + 1][1] - b - 1) if kk + 1 < len(meta) else -1
    cs_, ce_ = sorted([to_chrom(a), to_chrom(b)])
    if not bb:
        M.write(TAB.join([SP, aid, chrom, strand, str(p), str(N), str(cs_), str(ce_), str(ulen), str(gap),
                          "none", "NA", "NA", "NA", "NA", "NOALN", "NOALN"]) + NL)
        cnt["NOALN"] += 1
        continue
    toks = [(pp - CTX, k2, v2) for (pp, k2, v2) in bb[2] if CTX < pp <= CTX + U]
    cov0, cov1 = bb[3] - CTX, bb[4] - CTX
    ver = "OK" if (cov0 <= 1 and cov1 >= U) else "NO_FULL"
    tl = fmt(toks)
    ng = sum(1 for (pp, k2, v2) in toks if pp > U - gl)
    M.write(TAB.join([SP, aid, chrom, strand, str(p), str(N), str(cs_), str(ce_), str(ulen), str(gap), bb[1],
                      str(max(1, cov0)), str(min(U, cov1)), str(len(tl)), str(ng),
                      ",".join(tl) if tl else "-", ver]) + NL)
    cnt[ver] += 1
    nv.append(len(tl))
    nvg += ng
    for x in tl:
        allv[x] += 1
M.close()
nv.sort()
pf = os.path.join(OUT, "PARAMS.tsv")
newf = not os.path.exists(pf)
with io.open(pf, "a", encoding="utf8", newline=NL) as g3:
    if newf:
        g3.write(TAB.join(["species", "array_id", "chrom", "array_start", "array_end", "strand", "genes_per_unit",
                           "gene_hsp_kept", "period_med", "LO", "HI", "PAD", "CTX", "refunit_len", "gene_len",
                           "N_units", "n_break", "verify_OK", "verify_NO_FULL", "verify_NOALN",
                           "n_var_total", "n_var_gene_total", "var_sites", "median_var_per_unit"]) + NL)
    g3.write(TAB.join([SP, aid, chrom, str(a0), str(a1), strand, str(K), str(len(gg)), str(pm), str(LO), str(HI),
                       str(PAD), str(CTX), str(U), str(gl), str(N), str(brk), str(cnt["OK"]),
                       str(cnt["NO_FULL"]), str(cnt["NOALN"]), str(sum(nv)), str(nvg), str(len(allv)),
                       str(nv[len(nv) // 2]) if nv else "NA"]) + NL)
print("%s : k %d · 유전자HSP %d · 주기 %d · 참조 %d bp · 유닛 %d (끊김 %d) · OK %d · NO_FULL %d · NOALN %d · 변이자리 %d"
      % (aid, K, len(gg), pm, U, N, brk, cnt["OK"], cnt["NO_FULL"], cnt["NOALN"], len(allv)))
for x in (rp, qp, ctxp):
    try:
        os.remove(x)
    except OSError:
        pass
