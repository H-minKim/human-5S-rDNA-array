# -*- coding: utf-8 -*-
# 타 생물종 5S : 참조 유닛이 여럿인 array 의 G8 유닛 지도 (유닛마다 제일 맞는 참조를 쓴다)
#
# species_g8map.py (참조 하나) 를 그대로 옮기고, 참조가 여럿일 때 달라지는 곳만 고쳤다.
#
# 경계가 복잡해지는 지점은 딱 하나다 — "이 간격을 유닛으로 받을 것인가" 의 길이 창
#   유닛 경계 자체는 유전자 끝에서 나오므로 참조와 무관하다(유전자끝+1 .. 다음 유전자끝).
#   그런데 그 간격을 유닛으로 받을지 말지는 길이 창 LO~HI 로 정하고, 이 창은 참조의 주기에서 온다.
#   참조가 둘이면 창도 둘이다. 그래서 이 코드는
#     1) 그룹마다 제 주기(구성원 유닛 길이 중앙값)로 창 [LO_k, HI_k] 를 따로 만든다
#     2) 간격이 어느 한 창에라도 들면 유닛으로 받는다 (창의 합집합) · 어느 창이 받았는지 적는다
#     3) 배수 함정 검사 : 어떤 창의 HI 가 다른 창 LO 의 2 배 이상이면, 유전자를 하나 놓쳐
#        두 유닛이 합쳐진 조각이 '다른 종류의 정상 유닛' 으로 통과할 수 있다 → 경고를 적는다
#     4) 창이 겹치는 구간도 적는다 — 이 구간의 간격은 어느 종류인지 길이로는 못 가린다
#        (길이로 못 가려도 되는 이유 : 종류는 서열 정렬로 정하지 길이로 정하지 않는다)
#     5) 참조 하나(array 참조)만 썼을 때와 비교해 새로 살아난 유닛을 센다
#   길이 창은 합집합이므로 유닛이 사라지는 일은 없다(단일 창 ⊆ 합집합).
#
# 유닛마다 참조 고르기 : 모든 참조의 문맥판에 정렬해 점수 1 등. 확실성은 참조끼리 거리 Dmax 로 정규화
#   margin_frac >= 0.50 OK · 0.20~0.50 WEAK · < 0.20 AMBIG(→ array 참조 공통 프레임)
# 채움 금지 — 변이는 정렬(cs/BTOP)로만. 고른 참조의 유닛 구간을 다 못 덮으면 NO_FULL
#
# 인자 : array_id fasta gene.fa 참조폴더 그룹파일(쉼표) PAD CTX GL 출력폴더 array시작 array끝 chrom strand
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
(AID, FA, GENEF, RDIR, GFILES, PADS, CTXS, GLS, OUT, S0, E1, CHROM, STRAND) = sys.argv[1:14]
PAD, CTX, GL = int(PADS), int(CTXS), int(GLS)
s0, e1 = int(S0), int(E1)
os.makedirs(os.path.join(OUT, "g8map"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))


def rdfa(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()


def faidx(reg):
    r = subprocess.run(["samtools", "faidx", FA, reg], capture_output=True, text=True)
    return "".join(r.stdout.splitlines()[1:]).upper()


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


def fmt(tok):
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


def window(period):
    # 사람 G8v4 규칙의 종별 환산 : 주기 -6% / +7% 와 ±30 bp 중 더 너그러운 쪽
    return (int(min(period * 0.94, period - 30)), int(max(period * 1.07, period + 30)))


LEN = dict((L.split(TAB)[0], int(L.split(TAB)[1])) for L in io.open(FA + ".fai", encoding="utf8"))
refs = collections.OrderedDict()
refs["arrayref"] = rdfa(os.path.join(RDIR, AID + ".unit.fa"))
PERIOD = {"arrayref": len(refs["arrayref"])}
for f in [x for x in GFILES.split(",") if x]:
    k = os.path.basename(f).split(".")[-3]
    refs[k] = rdfa(f)
    h = io.open(f, encoding="utf8").read().split(NL)[0]
    PERIOD[k] = len(refs[k])
    for tok in h.split():
        if tok.startswith("period="):
            PERIOD[k] = int(tok[7:])
WIN = dict((k, window(PERIOD[k])) for k in refs)
GK = [k for k in refs if k != "arrayref"]
ULO = min(w[0] for w in WIN.values())
UHI = max(w[1] for w in WIN.values())
# 배수 함정 검사
warn = []
for a in WIN:
    for b in WIN:
        if a != b and WIN[a][1] >= 2 * WIN[b][0]:
            warn.append("%s 창 상한 %d 이 %s 창 하한 %d 의 2 배 이상" % (a, WIN[a][1], b, WIN[b][0]))
ov = (max(w[0] for w in WIN.values()), min(w[1] for w in WIN.values()))
ovtxt = "%d-%d" % ov if ov[0] <= ov[1] else "없음"

ctxp = os.path.join(OUT, "g8map", AID + ".ctx.fa")
with io.open(ctxp, "w", encoding="utf8", newline=NL) as fh:
    for k, s in refs.items():
        fh.write(">%s%s%s%s" % (k, NL, s[-CTX:] + s + s[:CTX], NL))


def aln_score(qseq, tseq):
    qp = os.path.join(OUT, "g8map", AID + ".q.fa")
    tp = os.path.join(OUT, "g8map", AID + ".t.fa")
    io.open(qp, "w", encoding="utf8", newline=NL).write(">q%s%s%s" % (NL, qseq, NL))
    io.open(tp, "w", encoding="utf8", newline=NL).write(
        ">t%s%s%s" % (NL, tseq[-CTX:] + tseq + tseq[:CTX], NL))
    r = subprocess.run("minimap2 -cx asm20 -t 2 %s %s 2>/dev/null" % (tp, qp),
                       shell=True, capture_output=True, text=True)
    b = 0
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        for x in (f[12:] if len(f) > 12 else []):
            if x.startswith("AS:i:"):
                b = max(b, int(x[5:]))
    os.remove(qp)
    os.remove(tp)
    return b


DM = 0
if len(GK) >= 2:
    DM = min((aln_score(refs[k], refs[k]) - max(aln_score(refs[k], refs[j]) for j in GK if j != k))
             for k in GK)

# ---- 유전자 위치로 경계 만들기 (참조와 무관) ----
mar = 3 * max(PERIOD.values()) + 3000
L0, R1 = max(1, s0 - mar), min(LEN[CHROM], e1 + mar)
seq = faidx("%s:%d-%d" % (CHROM, L0, R1))
rp = os.path.join(OUT, "g8map", AID + ".region.fa")
io.open(rp, "w", encoding="utf8", newline=NL).write(">region%s%s%s" % (NL, seq, NL))
rb = subprocess.run("blastn -task blastn -word_size 11 -dust no -evalue 1e-5 -query %s -subject %s "
                    "-outfmt '6 qstart qend sstart send bitscore' 2>/dev/null" % (GENEF, rp),
                    shell=True, capture_output=True, text=True)
G = []
for L in rb.stdout.splitlines():
    q1, q2, x1, x2, bs = L.split(TAB)
    q1, q2, x1, x2 = int(q1), int(q2), int(x1), int(x2)
    if x1 <= x2:
        G.append((x1 - (q1 - 1), x2 + (GL - q2), "+"))
    else:
        G.append((x2 - (GL - q2), x1 + (q1 - 1), "-"))
G.sort()
keep, last = [], -1
for (a, b, st) in G:
    if a > last:
        keep.append((a, b, st))
        last = b
major = collections.Counter(x[2] for x in keep).most_common(1)[0][0] if keep else STRAND
if major == "-":
    n = len(seq)
    wseq = rc(seq)
    gg = sorted([(n - b + 1, n - a + 1) for (a, b, st) in keep if st == major])
else:
    wseq = seq
    gg = sorted([(a, b) for (a, b, st) in keep if st == major])

qs, meta, brk = [], [], 0
pos = 0
n_rescued = 0
for i in range(len(gg) - 1):
    a, b = gg[i][1] + 1, gg[i + 1][1]
    L = b - a + 1
    acc = [k for k in refs if WIN[k][0] <= L <= WIN[k][1]]
    if not acc:
        brk += 1
        continue
    if "arrayref" not in acc:
        n_rescued += 1
    pos += 1
    la, lb = max(1, a - PAD), min(len(wseq), b + PAD)
    meta.append((pos, a, b, la, lb, L, ",".join(acc)))
    qs.append(">u%d%s%s%s" % (pos, NL, wseq[la - 1:lb], NL))
qp = os.path.join(OUT, "g8map", AID + ".units.fa")
io.open(qp, "w", encoding="utf8", newline=NL).write("".join(qs))

# ---- 유닛을 모든 참조에 정렬 ----
best = collections.defaultdict(dict)
rm = subprocess.run("minimap2 -cx asm20 --cs -N 10 -p 0.5 -t 4 %s %s 2>/dev/null" % (ctxp, qp),
                    shell=True, capture_output=True, text=True)
for L in rm.stdout.splitlines():
    f = L.split(TAB)
    if len(f) < 12:
        continue
    q, tgt, rs, re_ = f[0], f[5], int(f[7]), int(f[8])
    bs, cs = 0, ""
    for x in f[12:]:
        if x.startswith("AS:i:"):
            bs = int(x[5:])
        elif x.startswith("cs:Z:"):
            cs = x[5:]
    if tgt not in best[q] or bs > best[q][tgt][0]:
        best[q][tgt] = (bs, cs_vars(cs, rs), rs + 1, re_)

HDR = ["array_id", "chrom", "strand", "array_pos", "N", "chrom_start", "chrom_end", "unit_len", "gap_to_next",
       "len_window_accepted_by", "ref_used", "ref_best_raw", "score_best", "score_second", "margin_frac",
       "assign", "cov_start", "cov_end", "n_var", "n_var_gene", "varlist", "varlist_arrayref",
       "verify"] + ["score_%s" % k for k in refs]
M = io.open(os.path.join(OUT, "g8map", AID + ".g8mr.tsv"), "w", encoding="utf8", newline=NL)
M.write(TAB.join(HDR) + NL)
N = len(meta)
cnt, acnt, rcnt = collections.Counter(), collections.Counter(), collections.Counter()
for k, (p, a, b, la, lb, ulen, accby) in enumerate(meta):
    h = best.get("u%d" % p, {})
    gap = (meta[k + 1][1] - b - 1) if k + 1 < len(meta) else -1
    cs_, ce_ = (R1 - b + 1, R1 - a + 1) if major == "-" else (L0 + a - 1, L0 + b - 1)
    d = dict(array_id=AID, chrom=CHROM, strand=major, array_pos=p, N=N, chrom_start=cs_, chrom_end=ce_,
             unit_len=ulen, gap_to_next=gap, len_window_accepted_by=accby)
    for x in refs:
        d["score_%s" % x] = h[x][0] if x in h else 0

    def call(key):
        if key not in h:
            return None
        bs, toks, c0, c1 = h[key]
        U = len(refs[key])
        t = [(pp - CTX, kk, vv) for (pp, kk, vv) in toks if CTX < pp <= CTX + U]
        return (bs, fmt(t), max(1, c0 - CTX), min(U, c1 - CTX), (c0 - CTX <= 1 and c1 - CTX >= U),
                sum(1 for (pp, kk, vv) in t if pp > U - GL))

    ar = call("arrayref")
    d["varlist_arrayref"] = ",".join(ar[1]) if ar and ar[1] else ("-" if ar else "NOALN")
    cand = sorted(((h[x][0], x) for x in GK if x in h), reverse=True)
    if not cand:
        d.update(dict(ref_used="NOALN", ref_best_raw="NOALN", assign="NOALN", verify="NOALN",
                      score_best=0, score_second=0, margin_frac="-", cov_start="NA", cov_end="NA",
                      n_var="NA", n_var_gene="NA", varlist="NOALN"))
        cnt["NOALN"] += 1
    else:
        s1, k1 = cand[0]
        s2 = cand[1][0] if len(cand) > 1 else 0
        frac = (s1 - s2) / float(DM) if DM > 0 else 1.0
        lv = "OK" if frac >= 0.50 else ("WEAK" if frac >= 0.20 else "AMBIG")
        used = "arrayref" if lv == "AMBIG" else k1
        cc = call(used) or ar
        d.update(dict(ref_used=used, ref_best_raw=k1, score_best=s1, score_second=s2,
                      margin_frac="%.3f" % frac, assign=lv,
                      cov_start=cc[2], cov_end=cc[3], n_var=len(cc[1]), n_var_gene=cc[5],
                      varlist=",".join(cc[1]) if cc[1] else "-",
                      verify="OK" if cc[4] else "NO_FULL"))
        cnt[d["verify"]] += 1
        acnt[lv] += 1
        rcnt[used] += 1
    M.write(TAB.join(str(d.get(x, "-")) for x in HDR) + NL)
M.close()
for x in (rp, qp, ctxp):
    try:
        os.remove(x)
    except OSError:
        pass
S = io.open(os.path.join(OUT, "g8map", AID + ".g8mr.info.tsv"), "w", encoding="utf8", newline=NL)
S.write(TAB.join(["array_id", "n_units", "n_break", "n_rescued_by_group_window", "verify_OK", "verify_NO_FULL",
                  "verify_NOALN", "assign_OK", "assign_WEAK", "assign_AMBIG", "ref_pair_Dmax",
                  "windows", "window_overlap", "union_window", "multiple_period_warning"]) + NL)
S.write(TAB.join([AID, str(N), str(brk), str(n_rescued), str(cnt["OK"]), str(cnt["NO_FULL"]), str(cnt["NOALN"]),
                  str(acnt["OK"]), str(acnt["WEAK"]), str(acnt["AMBIG"]), str(DM),
                  " ".join("%s:%d-%d(주기%d)" % (k, WIN[k][0], WIN[k][1], PERIOD[k]) for k in refs),
                  ovtxt, "%d-%d" % (ULO, UHI), " · ".join(warn) if warn else "없음"]) + NL)
S.close()
print("%s : 유닛 %d (끊김 %d · 그룹 창이 살린 유닛 %d) · OK %d NO_FULL %d NOALN %d · 배정 %s · 참조 %s%s"
      % (AID, N, brk, n_rescued, cnt["OK"], cnt["NO_FULL"], cnt["NOALN"], dict(acnt), dict(rcnt),
         (" · 배수함정 " + warn[0]) if warn else ""))
