# -*- coding: utf-8 -*-
# 관통 read 에서 tier1 참조 서열
#   단계 : 참조 확정
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# tier1 참조 만들기 — 관통 read 의 array 구간을 하플로타입별 참조로 쓴다
#
#  하플로타입 묶기: (근위hap, 원위hap) 짝으로 묶는다.
#    근위 번호와 원위 번호는 서로 독립이라 한쪽만 쓰면 절반이 섞인다
#    (실측: 근위=원위인 read 가 46.9% 뿐). 관통 read 는 양쪽을 다 가지므로
#    짝이 하플로타입을 확정한다. 개체당 짝은 최대 2 종이어야 하고 실제로 97.2% 가 그렇다.
#
#  같은 하플로타입에 관통 read 가 여럿이면 array bp 가 가장 긴 것을 쓴다.
#    유닛 수 최대가 아니라 bp 최대다: 실측상 4 유닛 이상 어긋난 쌍은 bp 가 거의 같았다
#    (비 0.10~0.28) — 즉 서열은 멀쩡하고 blastn 이 유전자 HSP 를 놓친 계수 오류다.
#    서열을 참조로 쓸 것이므로 판단 기준은 실제 덮은 길이여야 한다.
#
#  array 경계는 CHM13 junction 바깥 50 kb 탐침 두 개로 잡는다 (다른 산출물과 동일).
import io
import os
import sys
import gzip
import subprocess
import collections

W = os.path.expanduser("~/tier1ref")
os.makedirs(W, exist_ok=True)
PROBE = os.path.expanduser("~/asmunits3/probe.fa")
GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
NL = "\n"

# 하플로타입 표 (side 를 키에 포함해야 한다)
H = collections.defaultdict(dict)
for hp in ("~/v5k/hap_reads_v5.tsv", "~/giabv5/hap_reads.tsv",
           "~/hg008v5/hap_reads.tsv"):
    p = os.path.expanduser(hp)
    if not os.path.exists(p):
        continue
    for i, L in enumerate(io.open(p, encoding="utf8")):
        if i == 0:
            continue
        f = L.rstrip(NL).split("\t")
        if len(f) >= 4:
            q = f[2].split("|")
            if len(q) >= 2:
                H[(f[0], q[1])][f[1]] = int(f[3])

SPAN = collections.defaultdict(list)
for L in io.open(os.path.expanduser("~/v5u/units_all.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) >= 7 and f[2] == "SPAN":
        SPAN[f[0]].append(f[1])


def find_reads(smp):
    for p in ("~/s5s_ont/hits/%s.hits.fastq.gz",
              "~/kgp_ont/out/%s.array.fastq.gz",
              "~/giab_hits/%s.hits.fastq.gz",
              "~/cancer5s/hits/%s.hits.fastq.gz"):
        q = os.path.expanduser(p % smp)
        if os.path.exists(q):
            return q
    return None


OUT = io.open(os.path.join(W, "tier1_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
for smp in sys.argv[1:]:
    if smp not in SPAN:
        continue
    q = find_reads(smp)
    if not q:
        continue
    want = set(SPAN[smp])
    seqs = {}
    with gzip.open(q, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline()
            fh.readline()
            rid = h[1:].rstrip(NL).replace(" ", "\t").split("\t")[0]
            if rid in want:
                seqs[rid] = s
    if not seqs:
        continue
    fa = os.path.join(W, "s_%s.fa" % smp)
    with io.open(fa, "w", newline=NL) as g:
        for r, s in seqs.items():
            g.write(u">%s\n%s\n" % (r, s))
    # 탐침은 -N 20 -p 0.1 로. --secondary=no 면 탐침이 read 하나에만 붙는다.
    paf = subprocess.run(["minimap2", "-x", "map-ont", "-t", TH, "-N", "20", "-p", "0.1",
                          fa, PROBE], capture_output=True, text=True)
    FL = collections.defaultdict(dict)
    for L in paf.stdout.splitlines():
        f = L.split("\t")
        if len(f) < 12:
            continue
        pq, t, ts, te, nm = f[0], f[5], int(f[7]), int(f[8]), int(f[9])
        if nm < 3000:
            continue
        cur = FL[t].get(pq)
        if cur is None or nm > cur[2]:
            FL[t][pq] = (ts, te, nm)
    os.remove(fa)

    cand = collections.defaultdict(list)   # (ph,dh) -> [(bp, read, s, e)]
    for r in seqs:
        d = FL.get(r, {})
        if "PROX50k" not in d or "DIST50k" not in d:
            continue
        hh = H.get((smp, r), {})
        if "prox" not in hh or "dist" not in hh:
            continue
        pm = (d["PROX50k"][0] + d["PROX50k"][1]) // 2
        dm = (d["DIST50k"][0] + d["DIST50k"][1]) // 2
        a = sorted([d["PROX50k"][0], d["PROX50k"][1],
                    d["DIST50k"][0], d["DIST50k"][1]])
        s0, e0 = a[1], a[2]
        if e0 - s0 < 4000:
            continue
        side = "PROX_FIRST" if pm < dm else "DIST_FIRST"
        cand[(hh["prox"], hh["dist"])].append((e0 - s0, r, s0, e0, side))
    for k in sorted(cand):
        cand[k].sort(reverse=True)
        bp, r, s0, e0, side = cand[k][0]
        seq = seqs[r][s0:e0]
        # 참조는 근위->원위 방향으로 통일한다 (다른 산출물과 좌표계를 맞춘다)
        if side == "DIST_FIRST":
            seq = seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]
        nm = "%s_p%d_d%d" % (smp, k[0], k[1])
        io.open(os.path.join(W, "t1_%s.fa" % nm), "w", newline=NL).write(
            ">%s|tier1|%s|%d%s%s%s" % (nm, r, bp, NL, seq, NL))
        OUT.write("%s\t%d\t%d\t%s\t%d\t%d\t%s%s"
                  % (smp, k[0], k[1], r, bp, len(cand[k]), side, NL))
        OUT.flush()
        sys.stderr.write("  %s (근위h%d,원위h%d) read %s · array %d bp · 후보 %d · %s\n"
                         % (smp, k[0], k[1], r[:12], bp, len(cand[k]), side))
        sys.stderr.flush()
OUT.close()
print("done")
