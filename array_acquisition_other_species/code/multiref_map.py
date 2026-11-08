# -*- coding: utf-8 -*-
# 그룹 참조가 여럿인 array : 유닛마다 가장 맞는 참조를 골라 변이를 다시 부른다
#
# 경계는 다시 자르지 않는다 (사용자 경고 대응)
#   G8 유닛 경계는 유전자 HSP 에서만 나온다(유전자끝+1 .. 다음 유전자끝). 참조 유닛은
#   (a) 길이 문턱 LO/HI 와 (b) 변이 호출용 정렬 에만 쓰인다. 그룹 참조 길이가 전부 LO/HI
#   안에 드는 것을 먼저 확인했으므로, 여기서는 기존 지도의 유닛 좌표를 그대로 읽어 쓴다.
#   따라서 참조가 늘어도 경계는 한 벌 그대로이고, 바뀌는 것은 '어느 참조로 변이를 불렀나' 뿐이다.
#
# 절차
#   1) 기존 g8map 에서 유닛 좌표(chrom_start/chrom_end/strand)를 읽는다
#   2) 그룹 참조마다 [끝 CTX][참조][앞 CTX] 문맥 참조를 만든다 (G8v4 규칙 · 채움 금지)
#   3) 유닛 서열을 한 번에 minimap2 로 모든 문맥 참조에 정렬 → 점수 최고를 선택 참조로
#   4) 선택 참조 프레임에서 cs 로 변이를 부른다. array 참조 프레임 변이는 그대로 보존
#   5) 2 등과 점수 차가 작으면 AMBIG 로 표시한다(어느 종류인지 못 가린 유닛)
#   6) 그룹 참조 <-> array 참조 자리 대응표(xwalk)를 따로 낸다 — 프레임이 다르므로 필수
#
# 인자 : array_id 지도파일 조립fasta 참조유닛폴더 그룹유닛파일(쉼표) CTX 출력폴더
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import re
import sys
import csv
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
AID, MAPF, FA, RDIR, GFILES, CTXS, OUT = sys.argv[1:8]
CTX = int(CTXS)
os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(OUT, "xwalk"), exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))


def rdfa(p):
    L = io.open(p, encoding="utf8").read().split(NL)
    return "".join(x.strip() for x in L[1:]).upper()


def cs_vars(cs, rs):
    # minimap2 cs 태그를 참조 자리 변이 토큰으로 (G8v4 와 같은 표기 · 채움 없음)
    p = rs
    out = []
    for m in re.finditer(r"([:*+-])([0-9]+|[a-z]+)", cs):
        op, v = m.group(1), m.group(2)
        if op == ":":
            p += int(v)
        elif op == "*":
            p += 1
            out.append("%d%s" % (p, v[1].upper()))
        elif op == "-":
            out.append("%dD%d" % (p + 1, len(v)))
            p += len(v)
        elif op == "+":
            out.append("%dI%d" % (p, len(v)))
    return out


rows = list(csv.DictReader(io.open(MAPF, encoding="utf8"), delimiter=TAB))
use = [r for r in rows if r["verify"] in ("OK", "NO_FULL") and r["chrom_start"].isdigit()]
GF = [x for x in GFILES.split(",") if x]
refs = collections.OrderedDict()
refs["arrayref"] = rdfa(os.path.join(RDIR, AID + ".unit.fa"))
for f in GF:
    refs[os.path.basename(f).split(".")[-3]] = rdfa(f)   # <array>.grpK.unit.fa -> grpK
cp = os.path.join(OUT, AID + ".ctxrefs.fa")
with io.open(cp, "w", encoding="utf8", newline=NL) as fh:
    for k, s in refs.items():
        fh.write(">%s%s%s%s" % (k, NL, s[-CTX:] + s + s[:CTX], NL))
reg = os.path.join(OUT, AID + ".regions.txt")
io.open(reg, "w", encoding="utf8", newline=NL).write(
    NL.join("%s:%s-%s" % (r["chrom"], r["chrom_start"], r["chrom_end"]) for r in use) + NL)
rr = subprocess.run(["samtools", "faidx", "-r", reg, FA], capture_output=True, text=True)
seqs, cur = [], []
for L in rr.stdout.splitlines():
    if L.startswith(">"):
        if cur:
            seqs.append("".join(cur))
        cur = []
    else:
        cur.append(L.strip())
if cur:
    seqs.append("".join(cur))
os.remove(reg)
up = os.path.join(OUT, AID + ".units.fa")
with io.open(up, "w", encoding="utf8", newline=NL) as fh:
    for i, (r, s) in enumerate(zip(use, seqs)):
        s = s.upper()
        fh.write(">u%d%s%s%s" % (i, NL, rc(s) if r["strand"] == "-" else s, NL))
# 참조끼리의 거리 Dmax : 참조를 짝마다 따로 정렬해 (자기점수 - 상대점수) 를 잰다.
#   유닛의 1 등-2 등 점수차를 이 Dmax 로 나누면 "어느 쪽인지 얼마나 확실한가" 가 된다.
#   고정 문턱(예 10 점)은 array 마다 참조 거리가 달라 자의적이다.
#   짝마다 따로 도는 이유 : arrayref 가 그룹 대표 중 하나와 서열이 같아서, 표적을 한꺼번에
#   주면 minimap2 가 같은 표적을 중복으로 보고 상대 점수를 아예 내놓지 않는다(실측 orang_bornean).
GK = [k for k in refs if k != "arrayref"]


def aln_score(qseq, tseq):
    qp = os.path.join(OUT, AID + ".q.fa")
    tp = os.path.join(OUT, AID + ".t.fa")
    io.open(qp, "w", encoding="utf8", newline=NL).write(">q%s%s%s" % (NL, qseq, NL))
    io.open(tp, "w", encoding="utf8", newline=NL).write(
        ">t%s%s%s" % (NL, tseq[-CTX:] + tseq + tseq[:CTX], NL))
    r = subprocess.run("minimap2 -cx asm20 -t 2 %s %s 2>/dev/null" % (tp, qp),
                       shell=True, capture_output=True, text=True)
    b = 0
    for L in r.stdout.splitlines():
        f = L.split(TAB)
        for x in f[12:] if len(f) > 12 else []:
            if x.startswith("AS:i:"):
                b = max(b, int(x[5:]))
    os.remove(qp)
    os.remove(tp)
    return b


DMAX = {}
for k in GK:
    own = aln_score(refs[k], refs[k])
    oth = [aln_score(refs[k], refs[j]) for j in GK if j != k]
    DMAX[k] = own - max(oth) if oth else 0
DM = min(v for v in DMAX.values()) if DMAX else 0
rm = subprocess.run("minimap2 -cx asm20 --cs -N 10 -p 0.5 -t 4 %s %s 2>/dev/null" % (cp, up),
                    shell=True, capture_output=True, text=True)
hit = collections.defaultdict(dict)
for L in rm.stdout.splitlines():
    f = L.split(TAB)
    if len(f) < 12:
        continue
    q, tgt, rs = f[0], f[5], int(f[7])
    bs, cs = 0, ""
    for x in f[12:]:
        if x.startswith("AS:i:"):
            bs = int(x[5:])
        elif x.startswith("cs:Z:"):
            cs = x[5:]
    if tgt not in hit[q] or bs > hit[q][tgt][0]:
        hit[q][tgt] = (bs, rs, cs, int(f[2]), int(f[3]))
HDR = ["species_key", "array_id", "chrom", "strand", "array_pos", "chrom_start", "chrom_end", "unit_len",
       "verify_orig", "ref_used", "ref_best_raw", "score_best", "score_second", "score_margin",
       "margin_frac", "ref_pair_Dmax", "assign",
       "n_var_bestref", "varlist_bestref", "varlist_arrayref"] + ["score_%s" % k for k in refs]
o = io.open(os.path.join(OUT, AID + ".multiref.tsv"), "w", encoding="utf8", newline=NL)
o.write(TAB.join(HDR) + NL)
cnt = collections.Counter()
for i, r in enumerate(use):
    h = hit.get("u%d" % i, {})
    d = dict((k, r.get(k, "-")) for k in ("chrom", "strand", "array_pos",
                                          "chrom_start", "chrom_end", "unit_len"))
    d["species_key"] = r.get("species", r.get("species_key", "-"))
    d.update({"array_id": AID, "verify_orig": r["verify"], "varlist_arrayref": r["varlist"],
              "ref_pair_Dmax": DM})
    for k in refs:
        d["score_%s" % k] = h[k][0] if k in h else 0
    cand = sorted(((h[k][0], k) for k in refs if k != "arrayref" and k in h), reverse=True)
    if not cand:
        d.update({"ref_used": "NOALN", "assign": "NOALN"})
        cnt["NOALN"] += 1
    else:
        s1, k1 = cand[0]
        s2 = cand[1][0] if len(cand) > 1 else 0
        bs, rs, cs, qs, ql = h[k1]
        vs = cs_vars(cs, rs - CTX) if cs else []
        vs = [v for v in vs if 1 <= int(re.match(r"([0-9]+)", v).group(1)) <= len(refs[k1])]
        frac = (s1 - s2) / float(DM) if DM > 0 else 0.0
        lv = "OK" if frac >= 0.50 else ("WEAK" if frac >= 0.20 else "AMBIG")
        if lv == "AMBIG":
            # 어느 종류인지 못 가린 유닛은 공통 프레임(array 참조)으로 둔다
            d.update({"ref_used": "arrayref", "n_var_bestref": "-", "varlist_bestref": "-"})
        else:
            d.update({"ref_used": k1, "n_var_bestref": len(vs),
                      "varlist_bestref": ",".join(vs) if vs else "-"})
        d.update({"ref_best_raw": k1, "score_best": s1, "score_second": s2, "score_margin": s1 - s2,
                  "margin_frac": "%.3f" % frac, "assign": lv})
        cnt[lv if lv != "OK" else k1] += 1
    o.write(TAB.join(str(d.get(k, "-")) for k in HDR) + NL)
o.close()
ar = os.path.join(OUT, AID + ".arrayref.fa")
io.open(ar, "w", encoding="utf8", newline=NL).write(">arrayref%s%s%s" % (NL, refs["arrayref"], NL))
for k in refs:
    if k == "arrayref":
        continue
    gp = os.path.join(OUT, AID + ".%s.fa" % k)
    io.open(gp, "w", encoding="utf8", newline=NL).write(">%s%s%s%s" % (k, NL, refs[k], NL))
    rr = subprocess.run("minimap2 -cx asm20 --cs -t 2 %s %s 2>/dev/null" % (ar, gp),
                        shell=True, capture_output=True, text=True)
    best = None
    for L in rr.stdout.splitlines():
        f = L.split(TAB)
        if len(f) < 12:
            continue
        sc, cs = 0, ""
        for x in f[12:]:
            if x.startswith("AS:i:"):
                sc = int(x[5:])
            elif x.startswith("cs:Z:"):
                cs = x[5:]
        if best is None or sc > best[0]:
            best = (sc, int(f[2]), int(f[7]), cs)
    w = io.open(os.path.join(OUT, "xwalk", "%s.%s_to_arrayref.tsv" % (AID, k)),
                "w", encoding="utf8", newline=NL)
    w.write(TAB.join(["grp_pos", "arrayref_pos", "op"]) + NL)
    if best:
        _sc, qp, tp, cs = best
        for m in re.finditer(r"([:*+-])([0-9]+|[a-z]+)", cs):
            op, v = m.group(1), m.group(2)
            if op == ":":
                for _ in range(int(v)):
                    qp += 1
                    tp += 1
                    w.write("%d%s%d%smatch%s" % (qp, TAB, tp, TAB, NL))
            elif op == "*":
                qp += 1
                tp += 1
                w.write("%d%s%d%smismatch%s" % (qp, TAB, tp, TAB, NL))
            elif op == "+":
                for _ in range(len(v)):
                    qp += 1
                    w.write("%d%s-%sins_in_grp%s" % (qp, TAB, TAB, NL))
            elif op == "-":
                for _ in range(len(v)):
                    tp += 1
                    w.write("-%s%d%sdel_in_grp%s" % (TAB, tp, TAB, NL))
    w.close()
    os.remove(gp)
os.remove(ar)
os.remove(up)
os.remove(cp)
print("%s : 유닛 %d · 참조거리 Dmax=%d · %s"
      % (AID, len(use), DM, " ".join("%s=%d" % (k, v) for k, v in cnt.most_common())))
