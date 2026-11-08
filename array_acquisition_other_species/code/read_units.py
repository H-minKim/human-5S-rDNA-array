# -*- coding: utf-8 -*-
# array 에 정렬된 read 를 유닛 단위로 나눈다 (정본 유닛 분할과 같은 경계)
#
#   경계는 G8 지도의 유닛 좌표(조립본 좌표)를 그대로 쓴다 — 정본을 자른 것과 같은 경계다.
#   read 마다 CIGAR 를 한 번 훑어 정렬 블록을 만들고, 유닛 구간과 겹치는 부분을 잘라낸다.
#
#   변이는 잘라낸 조각을 참조 유닛의 문맥판 `[끝 CTX][참조][앞 CTX]` 에 다시 정렬해 부른다
#     (G8 과 같은 프레임 · 채움 금지). 조립본 좌표로 부르면 프레임이 달라져 지도와 못 합친다.
#   보조/이차 정렬 제외(-F 0x900). mapq 는 적기만 하고 거르지 않는다(탠덤은 mapq0).
#   서열은 저장하지 않는다(용량) — 좌표·덮음·변이만 표로 낸다.
#
#   산출 : <out>/<species>__<sample>.read_unit_seq.tsv
#   인자 : 종키 샘플명 bam 지도폴더 참조유닛폴더 파라미터.tsv 출력폴더 [PAD기본=100]
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import re
import sys
import csv
import bisect
import subprocess
import collections
import pysam

NL = chr(10)
TAB = chr(9)
SP, SAMPLE, BAM, MAPD, RDIR, PARAM, OUT = sys.argv[1:8]
PADDEF = int(sys.argv[8]) if len(sys.argv) > 8 else 100
os.makedirs(OUT, exist_ok=True)
TMP = os.path.join(OUT, "_tmp_%s_%s" % (SP, SAMPLE))
os.makedirs(TMP, exist_ok=True)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
rc = lambda s: "".join(COMP.get(x, "N") for x in reversed(s))
PR = dict((r["array_id"], r) for r in csv.DictReader(io.open(PARAM, encoding="utf8"), delimiter=TAB))


def rdfa(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()


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


HDR = ["species_key", "sample", "read_id", "read_strand", "mapq", "array_id", "chrom", "array_strand",
       "array_pos", "unit_ref_start", "unit_ref_end", "unit_len", "cov_ref_start", "cov_ref_end",
       "covered_bp", "frac_covered", "q_start", "q_end", "seg_len", "spans_unit",
       "aln_to_refunit", "cov_start_in_ref", "cov_end_in_ref", "n_var", "varlist"]
o = io.open(os.path.join(OUT, "%s__%s.read_unit_seq.tsv" % (SP, SAMPLE)), "w", encoding="utf8", newline=NL)
o.write(TAB.join(HDR) + NL)
bam = pysam.AlignmentFile(BAM, "rb")
tot = collections.Counter()
for fn in sorted(os.listdir(MAPD)):
    if not fn.endswith(".g8map.tsv"):
        continue
    aid = fn[:-len(".g8map.tsv")]
    pr = PR.get(aid)
    rp = os.path.join(RDIR, aid + ".unit.fa")
    if not pr or not os.path.exists(rp):
        continue
    ref = rdfa(rp)
    U = len(ref)
    CTX = int(pr["CTX"])
    PAD = PADDEF
    ctxp = os.path.join(TMP, aid + ".ctx.fa")
    io.open(ctxp, "w", encoding="utf8", newline=NL).write(
        ">ctx%s%s%s" % (NL, ref[-CTX:] + ref + ref[:CTX], NL))
    units = []
    for r in csv.DictReader(io.open(os.path.join(MAPD, fn), encoding="utf8"), delimiter=TAB):
        if r["chrom_start"].isdigit():
            units.append((int(r["chrom_start"]), int(r["chrom_end"]), r["array_pos"],
                          r["chrom"], r["strand"]))
    if not units:
        continue
    chrom = units[0][3]
    astrand = units[0][4]
    # array 가 '-' 가닥이면 array_pos 순서와 조립본 좌표 순서가 반대다.
    #   구간 찾기는 좌표 오름차순으로 정렬한 목록에서 한다(bisect 전제).
    units.sort(key=lambda u: u[0])
    ustarts = [u[0] for u in units]
    maxlen = max(u[1] - u[0] + 1 for u in units)
    lo = min(u[0] for u in units)
    hi = max(u[1] for u in units)
    if chrom not in bam.references:
        continue
    segs, meta = [], []
    nseg = 0
    for rd in bam.fetch(chrom, max(0, lo - 1 - PAD), hi + PAD):
        if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
            continue
        q = rd.query_sequence
        if not q:
            continue
        # 정렬 블록 (ref 구간 -> query 구간)
        blocks = []
        rpos, qpos = rd.reference_start + 1, 0
        for op, ln in rd.cigartuples:
            if op in (0, 7, 8):
                blocks.append((rpos, rpos + ln - 1, qpos, qpos + ln - 1))
                rpos += ln
                qpos += ln
            elif op == 1:
                qpos += ln
            elif op in (2, 3):
                rpos += ln
            elif op == 4:
                qpos += ln
        if not blocks:
            continue
        bs = [b[0] for b in blocks]
        rs0, re0 = blocks[0][0], blocks[-1][1]
        i0 = bisect.bisect_right(ustarts, re0)
        for ui in range(max(0, bisect.bisect_left(ustarts, rs0 - maxlen - PAD)), i0):
            us, ue, apos, _c, _st = units[ui]
            if ue < rs0 or us > re0:
                continue
            w0, w1 = us - PAD, ue + PAD
            cov = 0
            qa = qb = None
            ca = cb = None
            for (b0, b1, p0, p1) in blocks:
                if b1 < w0 or b0 > w1:
                    continue
                x0, x1 = max(b0, w0), min(b1, w1)
                qq0, qq1 = p0 + (x0 - b0), p0 + (x1 - b0)
                if qa is None:
                    qa, ca = qq0, x0
                qb, cb = qq1, x1
                y0, y1 = max(b0, us), min(b1, ue)
                if y1 >= y0:
                    cov += y1 - y0 + 1
            if qa is None or cov < 20:
                continue
            seg = q[qa:qb + 1]
            if astrand == "-":
                seg = rc(seg)
            nseg += 1
            segs.append(">s%d%s%s%s" % (nseg, NL, seg, NL))
            meta.append((rd.query_name, "-" if rd.is_reverse else "+", rd.mapping_quality, aid, chrom,
                         astrand, apos, us, ue, ue - us + 1, ca, cb, cov,
                         "%.3f" % (cov / float(ue - us + 1)), qa, qb, len(seg),
                         "Y" if (ca <= us and cb >= ue) else "N"))
    if not segs:
        os.remove(ctxp)
        continue
    sp_ = os.path.join(TMP, aid + ".seg.fa")
    io.open(sp_, "w", encoding="utf8", newline=NL).write("".join(segs))
    best = {}
    rm = subprocess.run("minimap2 -cx map-ont --cs -t 8 %s %s 2>/dev/null" % (ctxp, sp_),
                        shell=True, capture_output=True, text=True)
    for L in rm.stdout.splitlines():
        f = L.split(TAB)
        if len(f) < 12:
            continue
        nm, rs, re_ = f[0], int(f[7]), int(f[8])
        sc, cs = 0, ""
        for x in f[12:]:
            if x.startswith("AS:i:"):
                sc = int(x[5:])
            elif x.startswith("cs:Z:"):
                cs = x[5:]
        if nm not in best or sc > best[nm][0]:
            best[nm] = (sc, rs, re_, cs)
    for i, m in enumerate(meta):
        b = best.get("s%d" % (i + 1))
        if not b:
            row = list(m) + ["NOALN", "NA", "NA", "NA", "NOALN"]
        else:
            _sc, rs, re_, cs = b
            toks = [(pp - CTX, kk, vv) for (pp, kk, vv) in cs_vars(cs, rs) if CTX < pp <= CTX + U]
            tl = fmt(toks)
            row = list(m) + ["mm2", str(max(1, rs + 1 - CTX)), str(min(U, re_ - CTX)), str(len(tl)),
                             ",".join(tl) if tl else "-"]
        o.write(TAB.join([SP, SAMPLE] + [str(x) for x in row]) + NL)
        tot["row"] += 1
        tot["noaln" if not b else "aln"] += 1
    os.remove(sp_)
    os.remove(ctxp)
    tot["array"] += 1
o.close()
bam.close()
try:
    os.rmdir(TMP)
except OSError:
    pass
print("%s %s : array %d · read×유닛 %d (정렬 %d · NOALN %d)"
      % (SP, SAMPLE, tot["array"], tot["row"], tot["aln"], tot["noaln"]))
