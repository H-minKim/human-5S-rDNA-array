# -*- coding: utf-8 -*-
# read x 유닛 메틸 행렬 — G8 틀 (유닛 경계를 G8 지도에서 그대로 가져온다)
#   옛 판(units_matrix.py)은 유닛 경계를 cmsearch 유전자 시작점으로 잡았다 = 사람이 버린 G6 틀.
#   이제 경계는 G8 지도의 chrom_start..chrom_end (유전자 끝 → 다음 유전자 끝) 를 쓴다.
#   그래서 read 행렬과 consensus map 이 (array_id, array_pos) 하나로 조인된다.
#   맥락(CG/CHG/CHH)은 참조 서열로 판정 · 호출은 mod_qual >= 0.5 · MAPQ 로 거르지 않는다
#   입력 : --fa 참조 · --mapdir G8 지도 폴더(<array_id>.g8map.tsv) · --samples · --bam({s}) · --meth({s}) · --out
#   산출 : <표본>.read_unit.tsv · 5S_unit_summary.tsv · 5S_array_summary.tsv  (열쇠 = array_id · array_pos)
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import csv
import bisect
import argparse
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
ap = argparse.ArgumentParser()
for k in ("--fa", "--mapdir", "--samples", "--bam", "--meth", "--out"):
    ap.add_argument(k, required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
# 1) G8 지도에서 유닛 구간 읽기
UNITS = collections.defaultdict(list)   # chrom -> [(start, end, array_id, array_pos, N)]
AINFO = {}
for fn in sorted(os.listdir(a.mapdir)):
    if not fn.endswith(".g8map.tsv"):
        continue
    for r in csv.DictReader(io.open(os.path.join(a.mapdir, fn), encoding="utf8"), delimiter=TAB):
        s, e = int(r["chrom_start"]), int(r["chrom_end"])
        UNITS[r["chrom"]].append((min(s, e), max(s, e), r["array_id"], int(r["array_pos"]), int(r["N"])))
        AINFO[r["array_id"]] = (r["chrom"], r["strand"], int(r["N"]))
for c in UNITS:
    UNITS[c].sort()
STARTS = dict((c, [x[0] for x in v]) for c, v in UNITS.items())
PAD = 5
SEQ = {}
for c, v in UNITS.items():
    lo, hi = v[0][0] - PAD, v[-1][1] + PAD
    r = subprocess.run(["samtools", "faidx", a.fa, "%s:%d-%d" % (c, max(1, lo), hi)], capture_output=True, text=True)
    SEQ[c] = (max(1, lo), "".join(r.stdout.splitlines()[1:]).upper())


def locate(c, pos1):
    v = UNITS.get(c)
    if not v:
        return None
    i = bisect.bisect_right(STARTS[c], pos1) - 1
    if i < 0:
        return None
    s, e, aid, ap_, N = v[i]
    return (aid, ap_, N) if s <= pos1 <= e else None


def context(c, pos1, strand):
    off, sq = SEQ.get(c, (1, ""))
    j = pos1 - off
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
US.write(TAB.join(["sample", "array_id", "chrom", "strand", "N", "array_pos", "context", "n_reads", "n_calls",
                   "frac_5mC", "mean_p5mC", "frac_5hmC"]) + NL)
AS = io.open(os.path.join(a.out, "5S_array_summary.tsv"), "w", encoding="utf8", newline=NL)
AS.write(TAB.join(["sample", "array_id", "chrom", "strand", "N", "context", "reads_any", "reads_mapq_ge20",
                   "calls", "frac_5mC", "frac_5mC_mapq_ge20", "read_level_frac_5mC_median", "frac_5hmC"]) + NL)
for smp in a.samples.split(","):
    MQ = {}
    r = subprocess.run("samtools view -F 0x900 %s | cut -f1,5" % a.bam.format(s=smp), shell=True, capture_output=True, text=True)
    for L in r.stdout.splitlines():
        x, q = L.split(TAB)
        MQ[x] = int(q)
    RU = collections.defaultdict(lambda: [0, 0, 0.0, 0])

    def add(rid, key, d):
        aid, ap_, cx = key
        e = RU[(rid, aid, ap_, cx)]
        pm = d.get("m", 0.0)
        e[0] += 1
        e[1] += (pm >= 0.5)
        e[2] += pm
        e[3] += (d.get("h", 0.0) >= 0.5)

    f = io.open(a.meth.format(s=smp), encoding="utf8")
    ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
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
            if site and site[3]:
                add(cur, site[3], codes)
            if rid in seen:
                sys.exit("입력이 read_id · ref_position 순으로 정렬돼 있지 않다: %s" % rid)
            if cur is not None:
                seen.add(cur)
            cur, site, codes, last = rid, None, {}, -2
        if pos0 < last:
            sys.exit("입력이 정렬돼 있지 않다 (위치 역행): %s" % rid)
        key = (x[ix["chrom"]], pos0, x[ix["ref_strand"]])
        if site is None or key != site[:3]:
            if site and site[3]:
                add(cur, site[3], codes)
            u = locate(key[0], pos0 + 1)
            site = (key[0], key[1], key[2], (u[0], u[1], context(key[0], pos0 + 1, key[2])) if u else None)
            codes = {}
        if site[3] is None:
            last = pos0
            continue
        codes[code] = float(x[ix["mod_qual"]])
        last = pos0
    if site and site[3]:
        add(cur, site[3], codes)
    O = io.open(os.path.join(a.out, smp + ".read_unit.tsv"), "w", encoding="utf8", newline=NL)
    O.write(TAB.join(["read_id", "mapq", "array_id", "array_pos", "context", "n_calls", "n_5mC", "frac_5mC",
                      "mean_p5mC", "n_5hmC"]) + NL)
    UA = collections.defaultdict(lambda: [set(), 0, 0, 0.0, 0])
    AR = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0]))
    for (rid, aid, ap_, cx), (n, m, sp_, h) in sorted(RU.items(), key=lambda t: (t[0][1], t[0][2], t[0][3], t[0][0])):
        O.write(TAB.join([rid, str(MQ.get(rid, -1)), aid, str(ap_), cx, str(n), str(m), "%.3f" % (m / float(n)),
                          "%.3f" % (sp_ / n), str(h)]) + NL)
        g = UA[(aid, ap_, cx)]
        g[0].add(rid)
        g[1] += n
        g[2] += m
        g[3] += sp_
        g[4] += h
        v = AR[(aid, cx)][rid]
        v[0] += n
        v[1] += m
        v[2] += h
    O.close()
    for (aid, ap_, cx), (rs, n, m, sp_, h) in sorted(UA.items()):
        c, st, N = AINFO[aid]
        US.write(TAB.join([smp, aid, c, st, str(N), str(ap_), cx, str(len(rs)), str(n), "%.3f" % (m / float(n)),
                           "%.3f" % (sp_ / n), "%.3f" % (h / float(n))]) + NL)
    for (aid, cx) in sorted(AR):
        c, st, N = AINFO[aid]
        rd = AR[(aid, cx)]
        tot = sum(v[0] for v in rd.values())
        me = sum(v[1] for v in rd.values())
        hm = sum(v[2] for v in rd.values())
        hq = [v for rid, v in rd.items() if MQ.get(rid, 0) >= 20]
        th, mh = sum(v[0] for v in hq), sum(v[1] for v in hq)
        fr = sorted(v[1] / float(v[0]) for v in rd.values() if v[0] >= 5)
        AS.write(TAB.join([smp, aid, c, st, str(N), cx, str(len(rd)), str(len(hq)), str(tot),
                           "%.3f" % (me / float(tot)) if tot else "NA", "%.3f" % (mh / float(th)) if th else "NA",
                           "%.3f" % fr[len(fr) // 2] if fr else "NA", "%.3f" % (hm / float(tot)) if tot else "NA"]) + NL)
    print("%s : read-유닛-맥락 %d 행 · read %d" % (smp, len(RU), len(set(k[0] for k in RU))))
US.close()
AS.close()
