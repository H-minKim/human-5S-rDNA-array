# -*- coding: utf-8 -*-
# 2 층 메틸 맵 (집계 + read 별 원본)
#   단계 : 배치·메틸
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 하플로타입 참조에 array read 를 붙여 메틸 맵을 만든다 (2 층 구조)
#
#  참조 = 그 개체의 하플로타입별 array 서열
#     tier1 : 관통 read 의 array 구간
#     tier2 : hifiasm 재조립 array   (HG002 에서 T2T Q100 과 100.0000% 일치 검증)
#
#  두 하플로타입을 한 참조에 함께 넣고 붙인다.
#    그래야 read 가 자기 하플로타입으로 간다 (일차 정렬이 곧 배정이다).
#    한쪽씩 따로 붙이면 반대쪽 read 도 억지로 붙는다.
#
#  MAPQ 로 거르지 않는다. 탠덤 array 에서 -q 20 은 array 창을 통째로 지운다
#    (2026-08 실측: read 63% 가 NO_ARRAY 로 오판됐다). 일차 정렬만 쓴다(-F 0x900).
#
#  MM/ML 은 pysam 의 modified_bases 로 읽는다. 역방향 정렬일 때 좌표를 스스로
#    바로잡아 준다. 직접 파싱하면 가닥을 틀리기 쉽다(이 프로젝트에서 2 회 사고).
#
#  산출 2 층
#    ① 집계층 methmap_agg  : 하플로타입 · 유닛마다 read 평균 β 와 표준편차
#    ② 개별층 methmap_read : read 하나하나의 유닛별 β (원본 메틸을 그대로 본다)
import io
import os
import sys
import gzip
import subprocess
import collections
import math

W = os.path.expanduser("~/methmap")
os.makedirs(W, exist_ok=True)
GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
NL = "\n"
import pysam


def read_fa(p):
    nm = None
    buf = []
    out = {}
    for L in io.open(p, encoding="utf8"):
        if L[0] == ">":
            if nm:
                out[nm] = "".join(buf)
            nm = L[1:].strip()
            buf = []
        else:
            buf.append(L.strip())
    if nm:
        out[nm] = "".join(buf)
    return out


def unit_bounds(seq, tag):
    """참조 array 를 유닛으로 자른다 -> [(start,end), ...]"""
    fa = os.path.join(W, "b_%s.fa" % tag)
    io.open(fa, "w", newline=NL).write(">x" + NL + seq + NL)
    db = os.path.join(W, "bd_%s" % tag)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                        "-evalue", "1e-5", "-max_target_seqs", "100000000",
                        "-num_threads", TH], capture_output=True, text=True, check=True)
    st = set()
    for L in r.stdout.splitlines():
        f = L.split("\t")
        if len(f) >= 2:
            st.add(min(int(f[0]), int(f[1])) - 1)
    K = sorted(st)
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    return [(K[i], K[i + 1]) for i in range(len(K) - 1)
            if 2100 <= K[i + 1] - K[i] <= 2400]


def find_reads(smp):
    for p in (os.path.expanduser("~/s5s_ont/hits/%s.hits.fastq.gz" % smp),
              os.path.expanduser("~/kgp_ont/out/%s.array.fastq.gz" % smp),
              os.path.expanduser("~/giab_hits/%s.hits.fastq.gz" % smp),
              os.path.expanduser("~/cancer5s/hits/%s.hits.fastq.gz" % smp)):
        if os.path.exists(p):
            return p
    return None


AGG = io.open(os.path.join(W, "methmap_agg_%s.tsv" % TAG), "a",
              encoding="utf8", newline=NL)
RDL = io.open(os.path.join(W, "methmap_read_%s.tsv" % TAG), "a",
              encoding="utf8", newline=NL)

for spec in sys.argv[1:]:
    # spec = "sample:ref1.fa[:ref2.fa]"  참조는 하플로타입별 array 서열
    parts = spec.split(":")
    smp = parts[0]
    refs = [p for p in parts[1:] if os.path.exists(p)]
    if not refs:
        sys.stderr.write("  참조 없음 %s\n" % smp)
        continue
    q = find_reads(smp)
    if not q:
        sys.stderr.write("  read 없음 %s\n" % smp)
        continue

    # 두 하플로타입을 한 참조로 합친다
    ref = os.path.join(W, "ref_%s.fa" % smp)
    UB = {}
    with io.open(ref, "w", newline=NL) as g:
        for hi, rp in enumerate(refs):
            d = read_fa(rp)
            nm0, sq = list(d.items())[0]
            name = "%s_hap%d" % (smp, hi + 1)
            g.write(u">%s\n%s\n" % (name, sq))
            UB[name] = (unit_bounds(sq, "%s_%d" % (smp, hi + 1)), len(sq), nm0)
    for name, (ub, L, nm0) in UB.items():
        sys.stderr.write("  참조 %s : %d bp · 유닛 %d  [%s]\n" % (name, L, len(ub), nm0[:50]))

    bam = os.path.join(W, "%s.bam" % smp)
    # -y 로 fastq 주석의 MM/ML 을 SAM 태그로 옮긴다 (minimap2 2.26 확인)
    p1 = subprocess.Popen(["minimap2", "-y", "-ax", "map-ont", "-t", TH,
                           "--secondary=no", ref, q], stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["samtools", "sort", "-@", "2", "-o", bam],
                          stdin=p1.stdout, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p2.wait()
    subprocess.run(["samtools", "index", bam], check=True, capture_output=True)

    bf = pysam.AlignmentFile(bam, "rb")
    nr = collections.Counter()
    per = collections.defaultdict(list)   # (hap, unit) -> [beta...]
    for a in bf.fetch(until_eof=False):
        if a.is_unmapped or a.is_secondary or a.is_supplementary:
            continue
        hap = a.reference_name
        nr[hap] += 1
        mb = None
        try:
            mb = a.modified_bases
        except Exception:
            mb = None
        if not mb:
            continue
        # 5mC 만 (C+m). 값은 0-255 확률.
        calls = []
        for k, v in mb.items():
            if len(k) >= 3 and (k[2] == "m" or k[2] == 77):
                calls.extend(v)
        if not calls:
            continue
        # read 좌표 -> 참조 좌표
        r2t = dict(a.get_aligned_pairs(matches_only=True))
        ub = UB[hap][0]
        bypos = collections.defaultdict(list)
        for rp, ql in calls:
            tp = r2t.get(rp)
            if tp is None:
                continue
            bypos[tp].append(ql / 255.0)
        if not bypos:
            continue
        # 유닛별로 모은다
        import bisect
        starts = [x[0] for x in ub]
        acc = collections.defaultdict(list)
        for tp, vs in bypos.items():
            i = bisect.bisect_right(starts, tp) - 1
            if i < 0 or i >= len(ub) or tp >= ub[i][1]:
                continue
            acc[i].extend(vs)
        for ui, vs in acc.items():
            b = sum(vs) / len(vs)
            per[(hap, ui)].append(b)
            RDL.write("%s\t%s\t%s\t%d\t%d\t%d\t%.4f\t%d%s"
                      % (smp, hap, a.query_name, ui + 1, ub[ui][0], ub[ui][1],
                         b, len(vs), NL))
    bf.close()
    RDL.flush()
    for hap in sorted(UB):
        ub, L, nm0 = UB[hap]
        for ui in range(len(ub)):
            vs = per.get((hap, ui), [])
            if not vs:
                AGG.write("%s\t%s\t%d\t%d\t%d\t0\t-\t-%s"
                          % (smp, hap, ui + 1, ub[ui][0], ub[ui][1], NL))
                continue
            m = sum(vs) / len(vs)
            sd = math.sqrt(sum((x - m) ** 2 for x in vs) / len(vs)) if len(vs) > 1 else 0.0
            AGG.write("%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f%s"
                      % (smp, hap, ui + 1, ub[ui][0], ub[ui][1], len(vs), m, sd, NL))
        sys.stderr.write("    %s : 붙은 read %d · 메틸 있는 유닛 %d/%d\n"
                         % (hap, nr[hap],
                            sum(1 for ui in range(len(ub)) if per.get((hap, ui))),
                            len(ub)))
    AGG.flush()
    sys.stderr.flush()
    os.remove(ref)
AGG.close()
RDL.close()
print("done")
