# -*- coding: utf-8 -*-
# 제브라피시 5S array : read x 유닛 CpG 메틸 행렬 (사람 HPRC 분석과 같은 형식)
#   유닛 경계 = GRCz12tu cmsearch 5S hit 시작점 (array 안에서 이웃 hit 시작 사이 = 유닛 하나) — bp/주기 환산 금지
#   입력 : modkit extract full (read · ref_position · mod_qual) · 개체 BAM(MAPQ) · array BED · hit 표
#   호출 : mod_qual >= 0.5 -> 메틸, < 0.5 -> 비메틸 (modkit 0.6.4 · C+m? 이라 m 확률만 있다)
#   MAPQ 로 거르지 않는다 (탠덤 array 는 mapq0 이 정상) — mapq 열로 남겨 나중에 고른다
#   산출 : <개체>.read_unit.tsv · 5S_unit_summary.tsv (개체 x array x 유닛) · 5S_array_summary.tsv
#   인자 : zf 폴더 · 개체 목록(쉼표)
import io
import os
import sys
import bisect
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
Z = sys.argv[1]
SAMPLES = sys.argv[2].split(",")
ARR = []
for L in io.open(os.path.join(Z, "grcz12_5S_arrays.bed"), encoding="utf8"):
    c, s, e = L.split()[:3]
    ARR.append((c, int(s), int(e)))
HIT = collections.defaultdict(list)
for L in io.open(os.path.join(Z, "zebrafish_GRCz12tu.tbl"), encoding="utf8"):
    if L.startswith("#"):
        continue
    f = L.split()
    a, b = int(f[7]), int(f[8])
    HIT[f[0]].append((min(a, b), f[9]))
UNITS = {}
for i, (c, s, e) in enumerate(ARR):
    st = sorted(x[0] for x in HIT[c] if s < x[0] <= e)
    strand = collections.Counter(x[1] for x in HIT[c] if s < x[0] <= e).most_common(1)[0][0]
    UNITS[i] = (c, s, e, st, strand)
BYC = collections.defaultdict(list)
for i, (c, s, e, st, strand) in UNITS.items():
    BYC[c].append((s, e, i))
for c in BYC:
    BYC[c].sort()


def locate(c, pos1):
    for s, e, i in BYC.get(c, ()):
        if s < pos1 <= e:
            st = UNITS[i][3]
            k = bisect.bisect_right(st, pos1) - 1
            if 0 <= k < len(st) - 1:
                return i, k + 1
            return i, None
    return None, None


US = io.open(os.path.join(Z, "5S_unit_summary.tsv"), "w", encoding="utf8", newline=NL)
US.write(TAB.join(["sample", "array_idx", "chrom", "array_start", "array_end", "strand", "unit_idx", "unit_start", "unit_end",
                   "n_reads", "n_cpg_calls", "frac_meth", "mean_prob"]) + NL)
AS = io.open(os.path.join(Z, "5S_array_summary.tsv"), "w", encoding="utf8", newline=NL)
AS.write(TAB.join(["sample", "array_idx", "chrom", "array_start", "array_end", "n_units", "reads_any", "reads_mapq_ge20",
                   "reads_ge10_units", "cpg_calls", "frac_meth", "frac_meth_mapq_ge20", "read_level_frac_meth_median"]) + NL)
for smp in SAMPLES:
    MQ = {}
    r = subprocess.run("samtools view -F 0x900 %s | cut -f1,5" % os.path.join(Z, "bam", smp + ".g12.5S.bam"),
                       shell=True, capture_output=True, text=True)
    for L in r.stdout.splitlines():
        a, q = L.split(TAB)
        MQ[a] = int(q)
    RU = collections.defaultdict(lambda: [0, 0, 0.0])      # (read, array, unit) -> n, n_meth, sum_prob
    f = io.open(os.path.join(Z, "meth", smp + ".5S.reads.tsv"), encoding="utf8")
    ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    for L in f:
        x = L.rstrip(NL).split(TAB)
        if x[ix["mod_code"]] != "m":
            continue
        i, u = locate(x[ix["chrom"]], int(x[ix["ref_position"]]) + 1)
        if i is None or u is None:
            continue
        p = float(x[ix["mod_qual"]])
        e = RU[(x[ix["read_id"]], i, u)]
        e[0] += 1
        e[1] += (p >= 0.5)
        e[2] += p
    O = io.open(os.path.join(Z, smp + ".read_unit.tsv"), "w", encoding="utf8", newline=NL)
    O.write(TAB.join(["read_id", "mapq", "array_idx", "chrom", "unit_idx", "n_cpg", "n_meth", "frac_meth", "mean_prob"]) + NL)
    UA = collections.defaultdict(lambda: [set(), 0, 0, 0.0])
    AR = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    for (rid, i, u), (n, m, sp) in sorted(RU.items(), key=lambda t: (t[0][1], t[0][2], t[0][0])):
        O.write(TAB.join([rid, str(MQ.get(rid, -1)), str(i), UNITS[i][0], str(u), str(n), str(m), "%.3f" % (m / float(n)), "%.3f" % (sp / n)]) + NL)
        a = UA[(i, u)]
        a[0].add(rid)
        a[1] += n
        a[2] += m
        a[3] += sp
        AR[i][rid][0] += n
        AR[i][rid][1] += m
    O.close()
    for (i, u), (rs, n, m, sp) in sorted(UA.items()):
        c, s, e, st, strand = UNITS[i]
        US.write(TAB.join([smp, str(i), c, str(s), str(e), strand, str(u), str(st[u - 1]), str(st[u] - 1), str(len(rs)), str(n),
                           "%.3f" % (m / float(n)), "%.3f" % (sp / n)]) + NL)
    UC = collections.Counter((rid, i) for (rid, i, u) in RU)
    for i, (c, s, e, st, strand) in UNITS.items():
        rd = AR.get(i, {})
        tot = sum(v[0] for v in rd.values())
        me = sum(v[1] for v in rd.values())
        hq = [v for rid, v in rd.items() if MQ.get(rid, 0) >= 20]
        th, mh = sum(v[0] for v in hq), sum(v[1] for v in hq)
        fr = sorted(v[1] / float(v[0]) for v in rd.values() if v[0] >= 5)
        AS.write(TAB.join([smp, str(i), c, str(s), str(e), str(max(0, len(st) - 1)), str(len(rd)), str(len(hq)),
                           str(sum(1 for (rid, j), k in UC.items() if j == i and k >= 10)), str(tot),
                           "%.3f" % (me / float(tot)) if tot else "NA", "%.3f" % (mh / float(th)) if th else "NA",
                           "%.3f" % fr[len(fr) // 2] if fr else "NA"]) + NL)
    print("%s : read-유닛 %d 행 · read %d" % (smp, len(RU), len(set(k[0] for k in RU))))
US.close()
AS.close()
