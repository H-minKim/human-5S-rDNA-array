# -*- coding: utf-8 -*-
# PacBio Fiber-seq read 배치 — pbontplace.py(ONT 판) 를 그대로 옮긴 판 (2026-09-12)
#
#  pbontplace.py 에서 바꾼 곳은 플랫폼에 따른 것뿐이다
#    ① 프리셋   minimap2 -x map-hifi   (ONT 는 map-ont)
#    ② 게이트   K=4 · Δ>=5 · MINARR=8,960 (4 유닛)
#               = HiFi 검정값 (hifitest/hifigate 2026-09-11, 211개체 read 3,824:
#                 정답 99.08% · 남의지도우위 0.00% · K=4 vs K=6/8/10 판정 갈림 0)
#               ONT 게이트(K=12 · MINARR 26,880)를 PacBio 에 쓰면 read 대부분이 빠진다
#                 (PacBio Fiber-seq read 중앙 ~17kb ≈ 7 유닛)
#    ③ read 이름  전체 이름을 적는다 (PacBio 는 movie/zmw/ccs 라 앞 16자가 서로 같다)
#    ④ 출력     ~/pbfsplace/pfjx_<TAG>.tsv · pfjv_<TAG>.tsv
#               pbu_<TAG>.tsv  = 재호출용 AO 17 열 표 (유닛의 read 좌표 f0·f1 포함)
#               <SAMD>/<sample>.<PBJOB>.sam = 그 read 들 (flag 0 · SEQ 는 참조 방향)
#               → recall_pad2.py AO pbu_*.tsv 출력.tsv  로 PAD=100 CTX=200 재호출한다.
#                 pfjv 의 변이 토큰은 배치용 지문일 뿐이다(문맥 없이 맨 정본유닛에 asm5).
#                 합본에 넣을 변이는 반드시 재호출본을 쓸 것.
#    ⑤ 지도 경로를 환경변수로 받는다 (기본값은 정본 그대로 — 안 주면 동작 동일)
#         MAPF  CONSENSUS_MAP   (hg002 는 슬롯2 를 T2T 부계로 고친 지도를 줘야 한다)
#         SUMF  CONSENSUS_SUM   (revcomp 는 열 이름으로 찾는다 — G6 는 6번째,
#                                G8v4 는 seq_file 열이 없어 5번째다. 2026-09-12 고침)
#         HSNF  hap_seq_node    · CMPD  array .fa 폴더
#  입력 ~/pbfs/out/<label>.hit.fq · <label>.bound.tsv   (label = <지도 개체>__<자료원>)
#  그 밖(유닛 자르기 = 정렬 좌표 · 변이 부르기 asm5 · 점수 |S∩V|−|S−V| · 귀무 실시간 · 유닛 전량 기록 scored)은 같다.
import io
import os
import re
import sys
import random
import subprocess
import collections

W = os.path.expanduser(os.environ.get("W", "~/pbfsplace"))
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
D = os.path.expanduser(os.environ.get("READD", "~/pbfs/out"))
NL = "\n"
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
FL = 5000
K = int(os.environ.get("K", "4"))
DMIN = int(os.environ.get("DMIN", "5"))
MINARR = int(os.environ.get("MINARR", "8960"))
MAXFL = int(os.environ.get("MAXFL", "999999999"))
PRESET = "map-hifi"
MAPF = os.path.expanduser(os.environ.get("MAPF", "~/CONSENSUS_MAP_G6.tsv"))
SUMF = os.path.expanduser(os.environ.get("SUMF", "~/CONSENSUS_SUM_G6.tsv"))
HSNF = os.path.expanduser(os.environ.get("HSNF", "~/hap_seq_node.tsv"))
CMPD = os.path.expanduser(os.environ.get("CMPD", "~/cmp4d"))
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
for x in io.open(HSNF, encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
RCI = None
for i, x in enumerate(io.open(SUMF, encoding="utf8")):
    f = x.rstrip(NL).split("\t")
    if i == 0:
        # revcomp 는 열 이름으로 찾는다. G6 SUM 은 6번째지만 G8v4 SUM 은
        #   seq_file 열이 없어 5번째다. 번호로 읽으면 off_prox 를 revcomp 로
        #   읽어 모든 array 가 조용히 정방향 처리된다.
        if "revcomp" not in f:
            raise SystemExit("SUMF 에 revcomp 열이 없다: %s" % SUMF)
        RCI = f.index("revcomp")
        continue
    if len(f) > RCI:
        RCF[(f[0], f[1])] = (f[RCI] == "Y")
UB = collections.defaultdict(list)
MAPV = collections.defaultdict(dict)
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


def tuc_targets(lab):
    out = set()
    p = os.path.join(D, "%s.bound.tsv" % lab)
    if not os.path.exists(p):
        return out
    f = io.open(p, encoding="utf8")
    hh = f.readline().rstrip(NL).split("\t")
    jx = dict((n, i) for i, n in enumerate(hh))
    if "arr_span" not in jx or "read" not in jx:
        sys.stderr.write("  열 없음 %s : %s%s" % (lab, ",".join(hh), NL))
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
OUT = io.open(os.path.join(W, "pfjx_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
VOU = io.open(os.path.join(W, "pfjv_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
# 재호출(recall_pad2.py)용 산출 — 유닛의 read 좌표 표 + 그 read 의 SAM
UOU = io.open(os.path.join(W, "pbu_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
PBJOB = os.environ.get("PBJOB", "PBFS")
DSRC = os.environ.get("DSRC", "PacBio_FiberSeq")
SAMD = os.path.expanduser(os.environ.get("SAMD", "~/aog8/sam"))
os.makedirs(SAMD, exist_ok=True)
SAMR = collections.defaultdict(list)
# others 는 개체마다 새로 만든다(아래 루프 안) — 모듈 수준 캐시는 결함이었다.
for lab in [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]:
    if lab in DONE:
        continue
    smp = lab.split("__")[0]
    want = tuc_targets(lab)
    slots = sorted(s for (a, s) in UB if a == smp)
    p = os.path.join(D, "%s.hit.fq" % lab)
    if not want or len(slots) < 2 or not (os.path.exists(p) and os.path.getsize(p) > 1000):
        sys.stderr.write("  %-18s 건너뜀 (TUC대상 %d · 지도 %s 슬롯 %d · read %s)%s"
                         % (lab, len(want), smp, len(slots), os.path.exists(p), NL))
        continue
    ref = os.path.join(W, "%s.ref.fa" % lab)
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
        sys.stderr.write("  %-18s 건너뜀 (array 파일 %d 개)%s" % (lab, len(ARR), NL))
        continue
    sub = os.path.join(W, "%s.fa" % lab)
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
    pr = subprocess.Popen(["minimap2", "-a", "-x", PRESET, "-t", TH, "-K", "100M",
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
                # f0·f1 은 SAM SEQ(참조 방향) 좌표 · f1 은 배타적(다음 유닛의 f0 과 같다).
                #   recall_pad2.py 에 이 좌표를 그대로 넘기고 SAM 은 flag 0 으로 적으므로
                #   그쪽 rev 분기를 타지 않는다 = 방향 뒤집기 사고가 원천적으로 없다.
                seg.append((u, seq[m2[A]:m2[B] + 1], m2[A], m2[B] + 1))
        if len(seg) >= K:
            REC.append((f[0], ctg, seg, int(f[4]), seq))
    pr.stdout.close()
    pr.wait()
    for q in (ref, sub):
        try:
            os.remove(q)
        except OSError:
            pass
    if not REC:
        io.open(_dp, "a", encoding="utf8", newline=NL).write(lab + NL)
        continue
    # 경쟁자 후보는 개체마다 새로 만든다 (2026-09-12 고침).
    #   옛 코드는 모듈 수준 others 를 한 번만 채우고 재사용해, 한 프로세스가 여러 개체를
    #   처리하면 두 번째부터는 후보에서 첫 개체가 빠지고 자기 자신이 들어갔다.
    #   자기 자신이 뽑히면 자기 지도끼리 비교해 Δ=0 이 되어 그 read 는 무조건 탈락한다.
    #   실행 분할에 따라 후보 목록이 달라져 결과가 재현되지 않던 원인이기도 하다.
    others = sorted(set(a for (a, s) in UB if a != smp))
    flat = [s for _r, _c, seg, _q, _sq in REC for _u, s, _a, _b in seg]
    V = varof(flat, "%s_%s" % (TAG, lab))
    j = 0
    npass = nfail = 0
    for rid, ctg, seg, mq, rseq in REC:
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
        # competitors are selected deterministically by read name.
        #   전역 random 을 쓰면 그 프로세스가 앞서 처리한 read 수에 따라 뽑히는 개체가
        #   달라져, 같은 입력이라도 실행 분할이 다르면 채택 목록이 바뀐다(실측 채택 집합
        #   일치율 86.6-90.8%). 어디에 놓이는지는 원래 결정적이었다(슬롯·유닛 불일치 0).
        #   검정 의미(경쟁자 하나와 비교)는 그대로라 K=4·Δ>=5 눈금을 다시 잡을 필요가 없다.
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
        N = UB[(smp, ctg)][0][3]
        SAMR[smp].append((rid, rseq))
        for t, (u, _s, f0_, f1_) in enumerate(seg):
            sc_ = 1 if t < K else 0
            # 재호출용 AO 17 열 표 (recall_pad2.py 입력과 같은 모양)
            #   1 개체 2 read전체 3 read짧은 4 슬롯 5 job 6 출처 7 - 8 array_pos
            #   9 array_pos 10 f0 11 f1 12 - 13 - 14 verify 15 varlist 16 - 17 Y
            #   verify·varlist 는 재호출이 덮어쓴다(여기서는 자리만 잡는다).
            UOU.write("\t".join([
                smp, rid, rid[:16], ctg, PBJOB, DSRC, "-",
                str(r0[0] + t), str(r0[0] + t), str(f0_), str(f1_),
                "RUN", "RUN", "RUN", "-", "-", "Y"]) + NL)
            OUT.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%s%s"
                      % (smp, rid, ctg, r0[0] + t, N + 1 - (r0[0] + t),
                         r0[1], nm, dlt, mq, sc_, lab, NL))
            vl = ",".join(sorted(ch[t][1])) if ch[t][1] else "-"
            VOU.write("%s\t%s\t%s\t%d\t%s\t%d\t%s%s"
                      % (smp, rid, ctg, r0[0] + t, vl, sc_, lab, NL))
    OUT.flush()
    VOU.flush()
    UOU.flush()
    # recall_pad2.py 가 읽을 SAM : flag 0 · SEQ 는 배치에 쓴 것 그대로(참조 방향).
    #   f0·f1 도 같은 좌표계라 그쪽에서 되돌릴 일이 없다.
    if SAMR[smp]:
        with io.open(os.path.join(SAMD, "%s.%s.sam" % (smp, PBJOB)),
                     "w", encoding="utf8", newline=NL) as sg:
            for rid_, sq_ in SAMR[smp]:
                sg.write("\t".join([rid_, "0", "*", "0", "0", "*", "*", "0", "0",
                                    sq_, "*"]) + NL)
        del SAMR[smp]
    io.open(_dp, "a", encoding="utf8", newline=NL).write(lab + NL)
    sys.stderr.write("  %-18s TUC대상 %-5d 정렬 %-5d · 채택 %-5d 탈락 %d%s"
                     % (lab, len(want), len(REC), npass, nfail, NL))
    sys.stderr.flush()
OUT.close()
VOU.close()
print("done")
