# -*- coding: utf-8 -*-
# Fiber-seq read 배치 — hifiasm 합의지도(CONSENSUS_MAP_G6)에 사슬 점수로 농는다 (2026-09-11)
#
#  왜 다시 하나
#    앞서는 관통 read 만 썼는데, Fiber-seq read 중 관통(양쪽 flank>=300bp)은
#    3.6% 뿐이다(실측). 나머지를 다 버렸다.
#    사슬 배치는 flank 없는 read 도 놓으므로 arr_span>=12유닛 인 31.6% 를 쓴다.
#
#  게이트는 일반 ONT 에서 검정된 것을 그대로 쓴다
#    K=12 · Δ>=10 · MINARR=26,880
#    짝지음 검정(같은 read 1,067개): K=12 97.52% vs K=16 97.63%,
#    Δ>=10 둘 다 통과한 957개 중 판정이 갈린 것 0건 · 귀무 둘 다 0.00%
#
#  MAXFL 기본값을 무제한으로 둔다
#    Fiber-seq read 는 HAP_READ_FINAL 에 없어 junction 배치 상대가 없다.
#    겹칠 상대가 없으므로 flank 있는 read 를 버릴 이유가 없다.
#
import io
import os
import re
import sys
import gzip
import glob
import random
import subprocess
import collections

W = os.path.expanduser(os.environ.get("W", "~/fsplace"))
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CMPD = os.path.expanduser(os.environ.get("CMPD", "~/cmp4d"))
NL = "\n"
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
FL = 5000
K = int(os.environ.get("K", "16"))
DMIN = int(os.environ.get("DMIN", "10"))
MINARR = int(os.environ.get("MINARR", "36000"))
MAXFL = int(os.environ.get("MAXFL", "999999999"))
CANON = os.path.join(R, "unit_canon.fa")
C = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
random.seed(20260909)


def rc(s):
    return s.translate(C)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


FP = rd(os.path.join(R, "FLANK_proximal_50000.fa"))
FD = rd(os.path.join(R, "FLANK_distal_50000.fa"))
SEQF = {}
HSNF = os.path.expanduser(os.environ.get("HSNF", "~/hap_seq_node.tsv"))
for x in io.open(HSNF, encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
SUMF = os.path.expanduser(os.environ.get("SUMF", "~/CONSENSUS_SUM_G6.tsv"))
for i, x in enumerate(io.open(SUMF, encoding="utf8")):
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
UB = collections.defaultdict(list)
MAPV = collections.defaultdict(dict)
MAPF = os.path.expanduser(os.environ.get("MAPF", "~/CONSENSUS_MAP_G6.tsv"))
f = io.open(MAPF, encoding="utf8")
h = f.readline().rstrip(NL).split("\t")
ix = dict((n, i) for i, n in enumerate(h))
for L in f:
    r = L.rstrip(NL).split("\t")
    if len(r) < 15 or r[ix["seg_type"]].startswith("PARTIAL"):
        continue
    k = (r[0], r[1])
    u = int(r[ix["array_pos"]])
    UB[k].append((u, int(r[ix["unit_start"]]) - 1, int(r[ix["unit_end"]]), int(r[ix["N"]])))
    MAPV[k][u] = set() if r[ix["varlist"]] == "-" else set(r[ix["varlist"]].split(","))
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r"([:*+-])([A-Za-z0-9]+)")


def tuc_targets(smp):
    """Fiber-seq 판정표에서 array 를 충분히 문 read 를 고른다.

    일반 ONT 판과 다른 점 — 판정표 위치만 다르다
      일반 ONT : ~/{s5s_ont_tuc,new44_tuc,kgp_tuc,giab_tuc,s5s_tuc}/out/<sample>.bound.tsv
      Fiber-seq: ~/fs41/out/<sample>.bound.tsv   (열 구성은 동일)

    flank 제한을 걸지 않는다 (MAXFL 기본 무제한)
      Fiber-seq read 에는 junction 배치가 없다 — HAP_READ_FINAL 이 일반 ONT read 로만
      만들어져 있어 Fiber-seq read 는 그 표에 없다. 겹칠 상대가 없으므로
      flank 있는 read 를 버릴 이유가 없다. 실측: 관통(양쪽 flank>=300bp) 은 3.6% 뿐이라
      flank 제한을 걸면 대부분을 잃는다.
    """
    out = set()
    p = os.path.expanduser("~/fs41/out/%s.bound.tsv" % smp)
    if not os.path.exists(p):
        return out
    f = io.open(p, encoding="utf8")
    hh = f.readline().rstrip(NL).split("\t")
    jx = dict((n, i) for i, n in enumerate(hh))
    if "arr_span" not in jx or "read" not in jx:
        sys.stderr.write("  열 없음 %s : %s%s" % (smp, ",".join(hh), NL))
        return out
    for L in f:
        r = L.rstrip(NL).split("\t")
        if len(r) < len(hh):
            continue
        try:
            if (int(r[jx["arr_span"]]) >= MINARR
                    and int(r[jx["flankL"]]) < MAXFL
                    and int(r[jx["flankR"]]) < MAXFL):
                out.add(r[jx["read"]])
        except ValueError:
            pass
    return out


def varof(seqs, tag):
    q = os.path.join(W, "u_%s.fa" % tag)
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


def sc(S, V):
    return len(S & V) - len(S - V)


def place(chain, M):
    keys = sorted(M)
    best = []
    for i0 in keys:
        if not all((i0 + t) in M for t in range(len(chain))):
            continue
        best.append((sum(sc(M[i0 + t], chain[t]) for t in range(len(chain))), i0))
    if len(best) < 2:
        return None
    best.sort(reverse=True)
    return best[0][1], best[0][0] - best[1][0]


DONE = set()
_dp = os.path.join(W, "done.txt")
if os.path.exists(_dp):
    DONE = set(x.strip() for x in io.open(_dp, encoding="utf8") if x.strip())
OUT = io.open(os.path.join(W, "fsjx_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
# 유닛별 변이 목록 — read-합의 지도 재료 (자리만으로는 서열 비교 불가)
VOU = io.open(os.path.join(W, "fsjv_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
# others 는 개체마다 새로 만든다(아래 루프 안) — 모듈 수준 캐시는 결함이었다.
for smp in [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]:
    if smp in DONE:
        continue
    want = tuc_targets(smp)
    slots = sorted(s for (a, s) in UB if a == smp)
    if not want or len(slots) < 2:
        continue
    # Fiber-seq read 는 압축되지 않은 .hit.fq 다 (헤더에 MM/ML 메틸 태그가 붙어 있다).
    #   read id 는 공백 분리 첫 토큰이라 태그가 있어도 그대로 읽힌다.
    p = os.path.expanduser("~/fs41/out/%s.hit.fq" % smp)
    if not (os.path.exists(p) and os.path.getsize(p) > 1000):
        continue
    ref = os.path.join(W, "%s.ref.fa" % smp)
    ARR = {}
    with io.open(ref, "w", newline=NL) as g:
        for sl in slots:
            q = os.path.join(CMPD, SEQF.get((smp, sl), "_"))
            if not os.path.exists(q):
                continue
            a = rd(q)
            if RCF.get((smp, sl)):
                a = rc(a)
            ARR[sl] = a
            g.write(">%s_h%s%s%s%s" % (smp, sl, NL, rc(FD[:FL]) + a + rc(FP[-FL:]), NL))
    if len(ARR) < 2:
        continue
    sub = os.path.join(W, "%s.fa" % smp)
    nw = 0
    with io.open(sub, "w", newline=NL) as g:
        with io.open(p, encoding="utf8") as fh:
            j = 0
            rid = None
            for L in fh:
                m = j & 3
                j += 1
                if m == 0:
                    rid = L[1:].split()[0]
                elif m == 1 and rid in want:
                    g.write(">%s%s%s%s" % (rid, NL, L.strip(), NL))
                    nw += 1
    if nw == 0:
        for q in (ref, sub):
            try:
                os.remove(q)
            except OSError:
                pass
        continue
    pr = subprocess.Popen(["minimap2", "-a", "-x", "map-ont", "-t", TH, "-K", "100M",
                           "--secondary=no", ref, sub],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    REC = []
    for L in pr.stdout:
        if L.startswith("@"):
            continue
        f = L.split("\t")
        if len(f) < 10 or f[2] == "*":
            continue
        if int(f[1]) & 0x900 or int(f[4]) < 20:
            continue
        ctg = f[2].rsplit("_h", 1)[1]
        a = ARR.get(ctg)
        if a is None:
            continue
        seq = f[9]
        rpos = int(f[3]) - 1
        q2 = 0
        m2 = {}
        for m in CIG.finditer(f[5]):
            k2, op = int(m.group(1)), m.group(2)
            if op == "S":
                q2 += k2
            elif op == "H":
                pass
            elif op in "M=X":
                for t in range(k2):
                    m2[rpos + t] = q2 + t
                rpos += k2
                q2 += k2
            elif op == "I":
                q2 += k2
            elif op in "DN":
                rpos += k2
        if not m2:
            continue
        seg = []
        for u, s0, e0, N in UB[(smp, ctg)]:
            A, B = s0 + FL, e0 + FL - 1
            if A in m2 and B in m2:
                seg.append((u, seq[m2[A]:m2[B] + 1]))
        if len(seg) >= K:
            REC.append((f[0], ctg, seg, int(f[4])))
    pr.stdout.close()
    pr.wait()
    for q in (ref, sub):
        try:
            os.remove(q)
        except OSError:
            pass
    if not REC:
        io.open(_dp, "a", encoding="utf8", newline=NL).write(smp + NL)
        continue
    # 개체마다 재계산 (v11 aoplace_v3fix 와 같은 수정)
    others = sorted(set(a for (a, s) in UB if a != smp))
    flat = [s for _r, _c, seg, _q in REC for _u, s in seg]
    V = varof(flat, "%s_%s" % (TAG, smp))
    j = 0
    npass = nfail = 0
    for rid, ctg, seg, mq in REC:
        ch = []
        for t in range(len(seg)):
            ch.append((seg[t][0], V.get(j)))
            j += 1
        if any(x[1] is None for x in ch):
            continue
        vs = [x[1] for x in ch[:K]]
        r0 = place(vs, MAPV[(smp, ctg)])
        if not r0:
            continue
        # 귀무 대조를 배치마다 건다 (사후가 아니라 실시간)
        o = random.Random(rid).choice(others)
        os_ = sorted(s for (a, s) in UB if a == o)
        nm = 0
        if os_:
            M2 = MAPV.get((o, os_[0]))
            if M2:
                r1 = place(vs, M2)
                if r1:
                    nm = r1[1]
        dlt = r0[1] - nm
        if dlt < DMIN:
            nfail += 1
            continue
        npass += 1
        # 배치 성공 — 이 read 가 덮은 유닛들을 기록
        N = UB[(smp, ctg)][0][3]
        # 유닛 기록을 K 에서 끊지 않는다 (2026-09-09 수정)
        #   자리는 앞 K 유닛으로 이미 확정됐고, 나머지 유닛은 연속이므로
        #   위치가 r0[0]+t 로 자동으로 따라온다. 그런데 원래는 버렸다.
        #   실측: read 가 실제로 걸친 유닛 중앙 23.6 개인데 16 개만 적어 35.9% 손실.
        #   단, K 번째 이후는 배치 점수에 쓰이지 않은 구간이므로
        #     scored 열로 구분해 표시한다 (검정 전까지 섞어 쓰지 않기 위해).
        for t, (u, _s) in enumerate(seg):
            # scored = 1 이면 배치 점수에 쓰인 앞 K 유닛, 0 이면 그 이후
            sc_ = 1 if t < K else 0
            OUT.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d%s"
                      % (smp, rid[:16], ctg, r0[0] + t, N + 1 - (r0[0] + t),
                         r0[1], nm, dlt, mq, sc_, NL))
            # read-합의 지도를 만들려면 변이 목록이 있어야 한다.
            #   자리(array_pos)만으로는 서열을 비교할 수 없다.
            vl = ",".join(sorted(ch[t][1])) if ch[t][1] else "-"
            VOU.write("%s\t%s\t%s\t%d\t%s\t%d%s"
                      % (smp, rid[:16], ctg, r0[0] + t, vl, sc_, NL))
    OUT.flush()
    VOU.flush()
    io.open(_dp, "a", encoding="utf8", newline=NL).write(smp + NL)
    sys.stderr.write("  %-10s TUC대상 %-5d 정렬 %-5d · 채택 %-5d 탈락 %d%s"
                     % (smp, len(want), len(REC), npass, nfail, NL))
    sys.stderr.flush()
OUT.close()
print("done")
