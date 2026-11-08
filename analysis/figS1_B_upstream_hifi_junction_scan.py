# -*- coding: utf-8 -*-
# minimap2 → JXHIFI_ALL.tsv.gz (대학 SLURM 클러스터)
#   단계 : S1 B 상류
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 검정 2 — HiFi flank(junction) 배치 -> HiFi-합의 지도 -> 변이 목록 재현 검정
#   서울대 클러스터 실행판 (2026-09-07)
#
#  무엇을 재나
#    지도(hifiasm/공개조립본)의 유닛별 변이 목록을, HiFi read 로 만든 합의가
#    되부르는가. 척도 = 재현율(지도 변이 중 되부른 비율) · 정밀도.
#    단순 서열 동일성이 아니라 변이 목록 재현이다.
#
#  왜 HiFi 가 심판이 되나
#    지도의 65%는 공개조립본(HiFi 포함), 23%는 우리 hifiasm(ONT 전용)이다.
#    HiFi read 는 ONT 와 오류 특성이 독립이므로, 두 계열을 따로 검정하면
#    "ONT 전용 조립본이 HiFi 로도 지지되는가" 를 가를 수 있다.
#
#  위치는 junction 앵커 — flank 서열로 고정. 지도를 참조하지 않는다.
#  하플로타입은 정렬이 정한다 (HiFi 는 정확도가 높고, 참조 자체가 조립본이라
#    ONT 에서 문제였던 라벨 뒤바뀜이 원리상 생기지 않는다).
import io
import os
import re
import sys
import subprocess
import collections

W = os.environ["WD"]
R = os.environ["REF"]
CAND = os.environ["CAND"]
ARR = os.path.join(W, "arrays")
TMP = os.environ.get("TMPD", os.path.join(W, "tmp"))
os.makedirs(TMP, exist_ok=True)
NL = "\n"
TH = os.environ.get("TH", "8")
TAG = os.environ.get("TAG", "x")
FL = 5000
MIND = int(os.environ.get("MIND", "3"))
CANON = os.path.join(R, "5S_canonical_unit.fa")
C = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r"([:*+-])([A-Za-z0-9]+)")


def rc(s):
    return s.translate(C)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


FP = rd(os.path.join(R, "FLANK_proximal_50000.fa"))
FD = rd(os.path.join(R, "FLANK_distal_50000.fa"))
SEQF = {}
SRC = {}
for i, x in enumerate(io.open(os.path.join(W, "hap_seq_node.tsv"), encoding="utf8")):
    f = x.rstrip(NL).split("\t")
    if i == 0 or len(f) < 3:
        continue
    SEQF[(f[0], f[1])] = f[2]
    SRC[(f[0], f[1])] = f[2].split("__")[0]
RCF = {}
for i, x in enumerate(io.open(os.path.join(W, "CONSENSUS_SUM_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
UB = collections.defaultdict(list)
VAR = {}
f = io.open(os.path.join(W, "CONSENSUS_MAP_G6.tsv"), encoding="utf8")
h = f.readline().rstrip(NL).split("\t")
ix = dict((n, i) for i, n in enumerate(h))
for L in f:
    r = L.rstrip(NL).split("\t")
    if len(r) < 15 or r[ix["seg_type"]].startswith("PARTIAL"):
        continue
    k = (r[0], r[1])
    u = int(r[ix["array_pos"]])
    UB[k].append((u, int(r[ix["unit_start"]]) - 1, int(r[ix["unit_end"]]), int(r[ix["N"]])))
    VAR[(r[0], r[1], u)] = set() if r[ix["varlist"]] == "-" else set(r[ix["varlist"]].split(","))


def varof(seqs, tag):
    q = os.path.join(TMP, "u_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm5", "-c", "--cs", "-t", TH,
                        "--secondary=no", CANON, q], capture_output=True, text=True)
    V = {}
    for L in r.stdout.splitlines():
        f = L.split("\t")
        if len(f) < 12:
            continue
        i = int(f[0][1:])
        if i in V:
            continue
        cs = ""
        for x in f[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        tp = int(f[7])
        out = set()
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == ":":
                tp += int(v)
            elif op == "*":
                out.add("%d%s" % (tp + 1, v[1].upper()))
                tp += 1
            elif op == "-":
                out.add("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                out.add("%dI%d" % (tp + 1, len(v)))
        V[i] = out
    os.remove(q)
    return V


DONE = set()
_dp = os.path.join(W, "out", "done_%s.txt" % TAG)
if os.path.exists(_dp):
    DONE = set(x.strip() for x in io.open(_dp, encoding="utf8") if x.strip())
OUT = io.open(os.path.join(W, "out", "jh_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)

for smp in [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]:
    if smp in DONE:
        continue
    p = os.path.join(CAND, "%s.fa" % smp)
    slots = sorted(s for (a, s) in UB if a == smp)
    if not os.path.exists(p) or not slots:
        continue
    ref = os.path.join(TMP, "%s.ref.fa" % smp)
    AR = {}
    with io.open(ref, "w", newline=NL) as g:
        for sl in slots:
            q = os.path.join(ARR, SEQF.get((smp, sl), "___"))
            if not os.path.exists(q):
                continue
            a = rd(q)
            if RCF.get((smp, sl)):
                a = rc(a)
            AR[sl] = a
            g.write(">%s_h%s%s%s%s" % (smp, sl, NL, rc(FD[:FL]) + a + rc(FP[-FL:]), NL))
    if not AR:
        continue
    pr = subprocess.Popen(["minimap2", "-a", "-x", "map-hifi", "-t", TH, "-K", "100M",
                           "--secondary=no", ref, p],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    ANC = collections.defaultdict(list)
    nanc = 0
    for L in pr.stdout:
        if L.startswith("@"):
            continue
        f = L.split("\t")
        if len(f) < 10 or f[2] == "*":
            continue
        if int(f[1]) & 0x900 or int(f[4]) < 20:
            continue
        sl = f[2].rsplit("_h", 1)[1]
        a = AR.get(sl)
        if a is None:
            continue
        seq = f[9]
        rpos = int(f[3]) - 1
        q2 = 0
        m2 = {}
        for m in CIG.finditer(f[5]):
            k, op = int(m.group(1)), m.group(2)
            if op == "S":
                q2 += k
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k):
                    m2[rpos + t] = q2 + t
                rpos += k
                q2 += k
            elif op == "I":
                q2 += k
            elif op in "DN":
                rpos += k
        if not m2:
            continue
        lo, hi = min(m2), max(m2)
        Ld, Lp = FL, FL + len(a)
        if not ((lo < Ld - 300 and hi > Ld + 300) or (lo < Lp - 300 and hi > Lp + 300)):
            continue
        nanc += 1
        for u, s0, e0, N in UB[(smp, sl)]:
            A, B = s0 + Ld, e0 + Ld - 1
            if A in m2 and B in m2:
                ANC[(sl, u)].append(seq[m2[A]:m2[B] + 1])
    pr.stdout.close()
    if pr.wait() != 0:
        sys.stderr.write("  MM_FAIL %s%s" % (smp, NL))
        continue
    keys = [k for k, v in ANC.items() if len(v) >= MIND]
    for (sl, u) in sorted(keys):
        lst = ANC[(sl, u)]
        V = varof(lst, "%s_%s_%s_%d" % (TAG, smp, sl, u))
        c = collections.Counter()
        for i in range(len(lst)):
            for x in V.get(i, set()):
                c[x] += 1
        tok = ";".join("%s:%d" % (x, c[x]) for x in sorted(c))
        M = VAR.get((smp, sl, u))
        N = UB[(smp, sl)][0][3]
        OUT.write("%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s%s"
                  % (smp, sl, N + 1 - u, u, len(lst), SRC.get((smp, sl), "?"),
                     tok or "-", (",".join(sorted(M)) or "-") if M is not None else "NA", NL))
    OUT.flush()
    try:
        os.remove(ref)
    except OSError:
        pass
    io.open(_dp, "a", encoding="utf8", newline=NL).write(smp + NL)
    sys.stderr.write("  %-10s 앵커 %4d · 자리 %d%s" % (smp, nanc, len(keys), NL))
    sys.stderr.flush()
OUT.close()
print("done")
