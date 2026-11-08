# -*- coding: utf-8 -*-
# 5S array : read x 유닛 메틸 행렬 — 종 공통판 (zf_units.py 일반화, 2026-09-11)
#   유닛 경계 = cmsearch 5S hit 시작점 (array 안 이웃 hit 시작 사이 = 유닛 하나) — bp/주기 환산 금지
#   맥락 = 참조 서열로 판정 (CG / CHG / CHH). + 가닥 C : 뒤 1 = G -> CG, 뒤 2 = G -> CHG, 아니면 CHH
#                                           - 가닥 (참조 G) : 앞 1 = C -> CG, 앞 2 = C -> CHG, 아니면 CHH
#   호출 = mod_qual >= 0.5 -> 그 변형 (m = 5mC · h = 5hmC). 한 자리에 m 행과 h 행이 따로 온다 (modkit extract full)
#   MAPQ 로 거르지 않는다 (탠덤 array 는 mapq0 이 정상) — mapq 열로 남긴다
#   입력 : --fa 참조 · --bed array BED(0-based) · --tbl cmsearch tblout(쉼표로 여럿) · --samples 쉼표
#          --bam 경로틀({s}=표본) · --meth modkit extract full tsv 경로틀({s}) · --out 산출 폴더
#   산출 : <표본>.read_unit.tsv · 5S_unit_summary.tsv · 5S_array_summary.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import bisect
import argparse
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
ap = argparse.ArgumentParser()
for k in ("--fa", "--bed", "--tbl", "--samples", "--bam", "--meth", "--out"):
    ap.add_argument(k, required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
ARR = []
for L in io.open(a.bed, encoding="utf8"):
    f = L.split()
    if len(f) >= 3:
        ARR.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else "%s:%s-%s" % tuple(f[:3])))
HIT = collections.defaultdict(list)
for t in a.tbl.split(","):
    for L in io.open(t, encoding="utf8"):
        if L.startswith("#"):
            continue
        f = L.split()
        x, y = int(f[7]), int(f[8])
        HIT[f[0]].append((min(x, y), f[9]))
UNITS = {}
SEQ = {}
PAD = 5
for i, (c, s, e, nm) in enumerate(ARR):
    st = sorted(set(x[0] for x in HIT[c] if s < x[0] <= e))
    cnt = collections.Counter(x[1] for x in HIT[c] if s < x[0] <= e)
    strand = cnt.most_common(1)[0][0] if cnt else "."
    UNITS[i] = (c, s, e, nm, st, strand)
    r = subprocess.run(["samtools", "faidx", a.fa, "%s:%d-%d" % (c, max(1, s + 1 - PAD), e + PAD)], capture_output=True, text=True)
    SEQ[i] = (max(0, s - PAD), "".join(r.stdout.splitlines()[1:]).upper())
BYC = collections.defaultdict(list)
for i, (c, s, e, nm, st, strand) in UNITS.items():
    BYC[c].append((s, e, i))


def locate(c, pos1):
    for s, e, i in BYC.get(c, ()):
        if s < pos1 <= e:
            st = UNITS[i][4]
            k = bisect.bisect_right(st, pos1) - 1
            if 0 <= k < len(st) - 1:
                return i, k + 1
            return i, None
    return None, None


def context(i, pos0, strand):
    off, sq = SEQ[i]
    j = pos0 - off
    if j < 2 or j + 2 >= len(sq):
        return "NA"
    if strand == "+":
        if sq[j] != "C":
            return "NA"
        return "CG" if sq[j + 1] == "G" else ("CHG" if sq[j + 2] == "G" else "CHH")
    if sq[j] != "G":
        return "NA"
    return "CG" if sq[j - 1] == "C" else ("CHG" if sq[j - 2] == "C" else "CHH")


US = io.open(os.path.join(a.out, "5S_unit_summary.tsv"), "w", encoding="utf8", newline=NL)
US.write(TAB.join(["sample", "array_idx", "array_name", "chrom", "array_start", "array_end", "strand", "unit_idx", "unit_start", "unit_end",
                   "context", "n_reads", "n_calls", "frac_5mC", "mean_p5mC", "frac_5hmC"]) + NL)
AS = io.open(os.path.join(a.out, "5S_array_summary.tsv"), "w", encoding="utf8", newline=NL)
AS.write(TAB.join(["sample", "array_idx", "array_name", "chrom", "array_start", "array_end", "n_units", "context", "reads_any", "reads_mapq_ge20",
                   "reads_ge10_units", "calls", "frac_5mC", "frac_5mC_mapq_ge20", "read_level_frac_5mC_median", "frac_5hmC"]) + NL)
for smp in a.samples.split(","):
    MQ = {}
    r = subprocess.run("samtools view -F 0x900 %s | cut -f1,5" % a.bam.format(s=smp), shell=True, capture_output=True, text=True)
    for L in r.stdout.splitlines():
        x, q = L.split(TAB)
        MQ[x] = int(q)
    # 흘려 읽기 (애기장대 모든 맥락은 수천만 행이라 자리를 다 메모리에 둘 수 없다)
    #   입력은 read_id · ref_position 순으로 정렬돼 있어야 한다 (호출 쪽에서 sort). 같은 자리의 m · h 행을 모아 한 자리로 센다
    #   정렬 검사 : 끝난 read 가 다시 나오거나 read 안에서 위치가 뒤로 가면 멈춘다 (조용한 이중 계산 방지)
    RU = collections.defaultdict(lambda: [0, 0, 0.0, 0])  # (read, array, unit, ctx) -> n, n_5mC, sum_p5mC, n_5hmC

    def add_site(rid, chrom, pos0, sd, d):
        i, u = locate(chrom, pos0 + 1)
        if i is None or u is None:
            return
        e = RU[(rid, i, u, context(i, pos0, sd))]
        pm = d.get("m", 0.0)
        e[0] += 1
        e[1] += (pm >= 0.5)
        e[2] += pm
        e[3] += (d.get("h", 0.0) >= 0.5)

    f = io.open(a.meth.format(s=smp), encoding="utf8")
    ix = dict((n, k) for k, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    cur, site, codes, last, seen = None, None, {}, -2, set()
    for L in f:
        x = L.rstrip(NL).split(TAB)
        code = x[ix["mod_code"]]
        if code not in ("m", "h"):
            continue
        pos0 = int(x[ix["ref_position"]])
        if pos0 < 0:
            continue
        rid = x[ix["read_id"]]
        if rid != cur:
            if site:
                add_site(cur, site[0], site[1], site[2], codes)
            if rid in seen:
                sys.exit("입력이 read_id · ref_position 순으로 정렬돼 있지 않다 (read 재등장): %s" % rid)
            if cur is not None:
                seen.add(cur)
            cur, site, codes, last = rid, None, {}, -2
        if pos0 < last:
            sys.exit("입력이 정렬돼 있지 않다 (read 안 위치 역행): %s" % rid)
        key = (x[ix["chrom"]], pos0, x[ix["ref_strand"]])
        if key != site:
            if site:
                add_site(cur, site[0], site[1], site[2], codes)
            site, codes = key, {}
        codes[code] = float(x[ix["mod_qual"]])
        last = pos0
    if site:
        add_site(cur, site[0], site[1], site[2], codes)
    O = io.open(os.path.join(a.out, smp + ".read_unit.tsv"), "w", encoding="utf8", newline=NL)
    O.write(TAB.join(["read_id", "mapq", "array_idx", "array_name", "chrom", "unit_idx", "context", "n_calls", "n_5mC", "frac_5mC", "mean_p5mC", "n_5hmC"]) + NL)
    UA = collections.defaultdict(lambda: [set(), 0, 0, 0.0, 0])
    AR = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0]))
    for (rid, i, u, cx), (n, m, sp, h) in sorted(RU.items(), key=lambda t: (t[0][1], t[0][2], t[0][3], t[0][0])):
        O.write(TAB.join([rid, str(MQ.get(rid, -1)), str(i), UNITS[i][3], UNITS[i][0], str(u), cx, str(n), str(m),
                          "%.3f" % (m / float(n)), "%.3f" % (sp / n), str(h)]) + NL)
        g = UA[(i, u, cx)]
        g[0].add(rid)
        g[1] += n
        g[2] += m
        g[3] += sp
        g[4] += h
        v = AR[(i, cx)][rid]
        v[0] += n
        v[1] += m
        v[2] += h
    O.close()
    for (i, u, cx), (rs, n, m, sp, h) in sorted(UA.items()):
        c, s, e, nm, st, strand = UNITS[i]
        US.write(TAB.join([smp, str(i), nm, c, str(s), str(e), strand, str(u), str(st[u - 1]), str(st[u] - 1), cx, str(len(rs)), str(n),
                           "%.3f" % (m / float(n)), "%.3f" % (sp / n), "%.3f" % (h / float(n))]) + NL)
    UC = collections.Counter((rid, i) for (rid, i, u, cx) in RU)
    for (i, cx) in sorted(AR):
        c, s, e, nm, st, strand = UNITS[i]
        rd = AR[(i, cx)]
        tot = sum(v[0] for v in rd.values())
        me = sum(v[1] for v in rd.values())
        hm = sum(v[2] for v in rd.values())
        hq = [v for rid, v in rd.items() if MQ.get(rid, 0) >= 20]
        th, mh = sum(v[0] for v in hq), sum(v[1] for v in hq)
        fr = sorted(v[1] / float(v[0]) for v in rd.values() if v[0] >= 5)
        AS.write(TAB.join([smp, str(i), nm, c, str(s), str(e), str(max(0, len(st) - 1)), cx, str(len(rd)), str(len(hq)),
                           str(sum(1 for (rid, j), k in UC.items() if j == i and k >= 10)), str(tot),
                           "%.3f" % (me / float(tot)) if tot else "NA", "%.3f" % (mh / float(th)) if th else "NA",
                           "%.3f" % fr[len(fr) // 2] if fr else "NA", "%.3f" % (hm / float(tot)) if tot else "NA"]) + NL)
    print("%s : read-유닛-맥락 %d 행 · read %d" % (smp, len(RU), len(set(k[0] for k in RU))))
US.close()
AS.close()
