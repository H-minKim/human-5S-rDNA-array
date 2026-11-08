# -*- coding: utf-8 -*-
# 5S 메틸 추출
#   입력 : 배치 BAM (junction ~/placed_bam · ARRAY_ONLY ~/placed_bam_g8) 의 MM/ML — read 원본 그대로
#          유닛 경계 : 지도 CONSENSUS_MAP (canonV8) 의 unit_start/unit_end
#   좌표 규칙 (파일럿에서 서열로 검증 · CpG 적중률 0.986)
#     · MM 은 원본 read 방향 기준이다. 역정렬 기록은 SAM SEQ 가 뒤집혀 있으므로
#       원본 좌표 op -> SEQ 좌표 sp = len(SEQ)-1-op 로 되돌린 뒤 CIGAR 로 참조 좌표에 옮긴다.
#     · 역방향 read 의 C 는 정본 forward 의 G 자리에 떨어진다 -> 1 을 빼서 CpG 의 C 로 맞춘다.
#       (보정 없이 쓰면 적중률 0.65 로 떨어진다 — 과거 64% 오좌표 사고와 같은 자리)
#     · junction BAM 은 [flank5000][정본][flank5000] 프레임이라 참조좌표 -5000, ARRAY_ONLY 는 그대로.
#   산출 (gzip)
#     METH_READUNIT_<노드>.tsv.gz : 개체·출처·설정·슬롯·유닛번호·read·가닥·CpG수·메틸수·beta·calls
#         calls = "유닛내좌표:값" 목록 (값 0-255, ML 그대로) -> 각 read 의 각 부위를 바로 볼 수 있다
#     METH_SITE_<노드>.tsv.gz     : 개체·출처·슬롯·유닛번호·유닛내좌표 별 집계 (n · n_meth · beta)
#   문턱 : 메틸 판정 ML >= 128 (확률 0.5)
#   인자 : 출력폴더 [개체목록파일]
import io
import os
import re
import sys
import gzip
import bisect
import subprocess
import collections

#  노드의 집 폴더. `ubuntu` 는 AMI 의 기본 계정 이름이고 사람 이름이 아니다.
H = os.environ.get("NODE_HOME", "/home/ubuntu")
NL = chr(10)
TAB = chr(9)
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
MLCUT = 128
OUT = sys.argv[1]
LIST = sys.argv[2] if len(sys.argv) > 2 else None
os.makedirs(OUT, exist_ok=True)
NODE = subprocess.run(["hostname"], capture_output=True, text=True).stdout.strip()
MAPF = os.path.join(H, "g8jxmap_v10", "CONSENSUS_MAP_G8_full.tsv")


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


# ── 유닛 경계 (정본 좌표 · 1-based 포함)
U = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] != "UNIT":
        continue
    U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]), int(r[ix["array_pos"]])))
for k in U:
    U[k].sort()
STARTS = dict((k, [x[0] for x in v]) for k, v in U.items())

SMPS = None
if LIST:
    SMPS = set(x.strip() for x in io.open(LIST, encoding="utf8") if x.strip())


def bams():
    for d, off, kind in ((os.path.join(H, "placed_bam"), 5000, "junction"),
                         (os.path.join(H, "NONE_placed_bam_g8"), 0, "array_only")):
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".bam"):
                continue
            p = fn.split(".")
            smp, cfg, src = p[0], p[1], p[2]
            if SMPS and smp not in SMPS:
                continue
            yield os.path.join(d, fn), smp, cfg, src, off, kind


def calls_of(seq, flag, mm, ml):
    """C+m 확률 (SEQ 좌표, 값). MM 은 원본 read 방향 기준."""
    orig = rc(seq) if (flag & 16) else seq
    vals = [int(x) for x in ml.split(",") if x != ""]
    out = []
    base = 0
    idx = None
    for blk in [p for p in mm.split(";") if p]:
        f2 = blk.split(",")
        deltas = [int(x) for x in f2[1:]]
        if f2[0].startswith("C+m"):
            if idx is None:
                idx = [i for i, c in enumerate(orig) if c == "C"]
            j = -1
            for t, d in enumerate(deltas):
                j += d + 1
                if j >= len(idx):
                    break
                op = idx[j]
                sp = (len(seq) - 1 - op) if (flag & 16) else op
                out.append((sp, vals[base + t] if base + t < len(vals) else 0))
        base += len(deltas)
    return out


RU = gzip.open(os.path.join(OUT, "METH_READUNIT_%s.tsv.gz" % NODE), "wt", encoding="utf8", newline=NL)
RU.write(TAB.join(["sample", "data_source", "cfg", "kind", "slot", "array_pos", "read_id",
                   "strand", "n_cpg", "n_meth", "beta", "calls"]) + NL)
SITE = collections.Counter()
SITEM = collections.Counter()
T = collections.Counter()

for b, smp, cfg, src, off, kind in bams():
    pr = subprocess.Popen(["samtools", "view", b], stdout=subprocess.PIPE, text=True)
    for L in pr.stdout:
        r = L.rstrip(NL).split(TAB)
        if len(r) < 11:
            continue
        flag = int(r[1])
        if flag & 0x900 or r[2] == "*":
            continue
        sn, pos, cig, seq = r[2], int(r[3]), r[5], r[9]
        sl = sn.split("_h")[-1]
        key = (smp, sl)
        if key not in U:
            T["슬롯 지도 없음"] += 1
            continue
        mm = ml = None
        for x in r[11:]:
            if x.startswith("MM:Z:"):
                mm = x[5:]
            elif x.startswith("ML:B:C,"):
                ml = x[7:]
        if not mm or not ml:
            T["MM/ML 없음"] += 1
            continue
        m = {}
        rp = pos - 1
        q = 0
        for g in CIG.finditer(cig):
            k, op = int(g.group(1)), g.group(2)
            if op == "S":
                q += k
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k):
                    m[q + t] = rp + t
                rp += k
                q += k
            elif op == "I":
                q += k
            elif op in "DN":
                rp += k
        st = "-" if (flag & 16) else "+"
        per = collections.defaultdict(list)
        for sp, v in calls_of(seq, flag, mm, ml):
            if sp not in m:
                continue
            cp = m[sp] - off                      # 정본 0-based
            if flag & 16:
                cp -= 1                           # 역방향 C -> CpG 의 C
            p1 = cp + 1                           # 정본 1-based
            if p1 < 1:
                continue
            S = STARTS[key]
            i = bisect.bisect_right(S, p1) - 1
            if i < 0:
                continue
            s0, e0, ap = U[key][i]
            if p1 > e0:
                continue
            per[ap].append((p1 - s0 + 1, v))
        T["read"] += 1
        for ap, lst in per.items():
            lst.sort()
            nm = sum(1 for _, v in lst if v >= MLCUT)
            RU.write(TAB.join([smp, src, cfg, kind, sl, str(ap), r[0], st, str(len(lst)), str(nm),
                               "%.4f" % (nm / float(len(lst))),
                               ",".join("%d:%d" % (a, v) for a, v in lst)]) + NL)
            T["read-유닛"] += 1
            for a, v in lst:
                SITE[(smp, src, sl, ap, a)] += 1
                if v >= MLCUT:
                    SITEM[(smp, src, sl, ap, a)] += 1
    pr.stdout.close()
    pr.wait()
RU.close()

SI = gzip.open(os.path.join(OUT, "METH_SITE_%s.tsv.gz" % NODE), "wt", encoding="utf8", newline=NL)
SI.write(TAB.join(["sample", "data_source", "slot", "array_pos", "unit_pos", "n", "n_meth", "beta"]) + NL)
for k in sorted(SITE):
    n = SITE[k]
    nm = SITEM.get(k, 0)
    SI.write(TAB.join([k[0], k[1], k[2], str(k[3]), str(k[4]), str(n), str(nm),
                       "%.4f" % (nm / float(n))]) + NL)
SI.close()
for k, v in sorted(T.items()):
    print("  %-16s %d" % (k, v))
print("  자리 %d" % len(SITE))
