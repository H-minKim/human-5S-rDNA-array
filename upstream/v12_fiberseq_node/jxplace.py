# -*- coding: utf-8 -*-
"""junction 배치 — 슬롯이 정해진 read 의 위치(array_pos)를 flank 앵커로 정한다 (2026-09-15)

  들어가기 전 : 슬롯은 slotvote.py 가 flank 변이 투표로 정했다(표차>=0.90 에서 정확도 99.3%).
                그러므로 여기서는 사슬 게이트(K=12)를 쓰지 않는다 — 그 게이트는 슬롯을
                정하려고 있던 것이고, 이미 정해졌으므로 유닛 1 개짜리 read 도 위치를 받는다.

  절차
    ① read 를 자기 슬롯의 flank 붙인 정본에 정렬 (map-ont · --secondary=no)
    ② 정렬 좌표 → 정본 좌표(− flank5) → G6 지도의 유닛 구간을 찾는다
    ③ 유닛마다 read 서열을 잘라 정본유닛(unit_canon.fa)에 정렬해 변이를 부른다
    ④ AOJV 형식 6 열 출력 : 개체 · read(16자) · 슬롯 · array_pos · 변이목록 · scored

  Coordinate convention : the coordinate in the flank-added file = the canonical coordinate + flank5
"""
import io
import os
import sys
import subprocess
import collections
import re

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r"([:*+\-])([A-Za-z0-9]+)")
H = os.path.expanduser
HSN = H(os.environ.get("HSN", "~/hap_seq_node.tsv"))
CMPD = H(os.environ.get("CMPD", "~/cmp4d_flank"))
MAPF = H(os.environ.get("MAPF", "~/CONSENSUS_MAP_G6.tsv"))
READS = H(os.environ.get("READS", "~/fs41/out"))
CANON = H("~/s5s_tuc/ref/unit_canon.fa")
SUMF = H(os.environ.get("SUMF", "~/CONSENSUS_SUM_G6.tsv"))
C_TR = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(x):
    return x.translate(C_TR)[::-1]

TH = os.environ.get("TH", "6")
MARGIN = float(os.environ.get("MARGIN", "0.90"))
# read 이름 길이 : ONT 16 · PacBio 는 0(전체) — 앞 16 자가 movie 라 전건 충돌
RIDLEN = int(os.environ.get("RIDLEN", "16"))
VOTE, OUT = sys.argv[1], sys.argv[2]
NL = "\n"

SEQF = {}
for x in io.open(HSN, encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3:
        SEQF[(f[0], f[1])] = f[2]

# revcomp 표 — 열 이름으로 찾는다(G6 는 6번째, G8 계열은 5번째)
RCF = {}
for i, x in enumerate(io.open(SUMF, encoding="utf8")):
    fz = x.rstrip(NL).split("	")
    if i == 0:
        if "revcomp" not in fz:
            raise SystemExit("SUMF 에 revcomp 열이 없다: %s" % SUMF)
        RCI = fz.index("revcomp")
        continue
    if len(fz) > RCI:
        RCF[(fz[0], fz[1])] = (fz[RCI] == "Y")

# 지도 : (개체,슬롯) -> [(array_pos, start0, end, N)]   ※ 정본 좌표(flank 없음)
UB = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
h = f.readline().rstrip(NL).split("\t")
ix = dict((n, i) for i, n in enumerate(h))
for L in f:
    r = L.rstrip(NL).split("\t")
    if len(r) < 15 or r[ix["seg_type"]].startswith("PARTIAL"):
        continue
    UB[(r[0], r[1])].append((int(r[ix["array_pos"]]), int(r[ix["unit_start"]]) - 1,
                             int(r[ix["unit_end"]]), int(r[ix["N"]])))
for k in UB:
    UB[k].sort()
CANLEN = dict((k, max(e for (_u, _s0, e, _N) in v)) for k, v in UB.items())


def load(p):
    hdr = ""
    buf = []
    for L in io.open(p, encoding="utf8"):
        if L.startswith(">"):
            hdr = L.strip()
        else:
            buf.append(L.strip())
    d = dict(t.split("=", 1) for t in hdr.split("|")[1:] if "=" in t)
    return "".join(buf).upper(), int(d.get("flank5", 0))


# 슬롯 배정 읽기
ASSIGN = collections.defaultdict(dict)      # 개체 -> read -> 슬롯
for i, L in enumerate(io.open(VOTE, encoding="utf8")):
    if i == 0:
        continue
    r = L.rstrip(NL).split("\t")
    if len(r) < 8:
        continue
    if float(r[7]) < MARGIN:
        continue
    ASSIGN[r[0]][r[1]] = "1" if int(r[5]) > int(r[6]) else "2"

o = io.open(OUT, "w", encoding="utf8", newline=NL)
T = collections.Counter()

for smp in sorted(ASSIGN):
    rd2sl = ASSIGN[smp]
    fq = os.path.join(READS, "%s.hit.fq" % smp)
    if not os.path.exists(fq):
        continue
    SEQ = {}
    it = io.open(fq, encoding="utf8")
    while True:
        hh = it.readline()
        if not hh:
            break
        sq = it.readline().strip()
        it.readline()
        it.readline()
        rid = hh[1:].split()[0].split("\t")[0]
        if rid in rd2sl:
            SEQ[rid] = sq
    if not SEQ:
        continue

    for sl in ("1", "2"):
        ids = [r for r in SEQ if rd2sl[r] == sl]
        if not ids:
            continue
        key = (smp, sl)
        if key not in UB or (smp, sl) not in SEQF:
            continue
        p = os.path.join(CMPD, SEQF[(smp, sl)])
        if not os.path.exists(p):
            continue
        s, f5 = load(p)
        # 극성 : 지도의 unit_start/end 는 revcomp 적용 후 좌표계다(fsplace 206-210 행과 같게).
        #   flank 정본은 원래 방향으로 저장돼 있으므로 revcomp=Y 슬롯은 여기서 뒤집어야 한다.
        #   안 뒤집으면 유닛 번호가 N+1-p 로 통째로 뒤집힌다(2026-09-15 실측으로 확인).
        if RCF.get((smp, sl)):
            f3 = len(s) - f5 - CANLEN.get((smp, sl), 0)
            s = rc(s)
            f5 = f3                      # 뒤집으면 5'flank 자리에 3'flank 가 온다
        ref = "/tmp/jx_%s_%s.fa" % (smp, sl)
        io.open(ref, "w", newline=NL).write(">%s_h%s%s%s%s" % (smp, sl, NL, s, NL))
        q = "/tmp/jx_%s_%s.q.fa" % (smp, sl)
        with io.open(q, "w", newline=NL) as g:
            for r in ids:
                g.write(">%s%s%s%s" % (r, NL, SEQ[r], NL))

        # ① read → 자기 슬롯 참조
        pr = subprocess.Popen(["minimap2", "-a", "-x", "map-ont", "-t", TH,
                               "--secondary=no", ref, q],
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        cut = []        # (read, array_pos, 조각서열)
        for L in pr.stdout:
            if L.startswith("@"):
                continue
            fl = L.rstrip(NL).split("\t")
            if len(fl) < 10 or fl[2] == "*":
                continue
            if int(fl[1]) & 0x900:
                continue
            rid = fl[0]
            seq = fl[9]
            rpos = int(fl[3]) - 1        # 참조(= 정본+flank5) 좌표
            qp = 0
            m2 = {}
            for mm in CIG.finditer(fl[5]):
                n_, op = int(mm.group(1)), mm.group(2)
                if op == "S":
                    qp += n_
                elif op == "H":
                    pass
                elif op in "M=X":
                    for t in range(n_):
                        m2[rpos + t] = qp + t
                    rpos += n_
                    qp += n_
                elif op == "I":
                    qp += n_
                elif op in "DN":
                    rpos += n_
            if not m2:
                continue
            # ② 지도 유닛 구간을 참조 좌표로 옮겨(+f5) read 가 다 덮는 유닛만
            for (u, s0, e0, N) in UB[key]:
                A, B = s0 + f5, e0 + f5 - 1
                if A in m2 and B in m2:
                    cut.append((rid, u, seq[m2[A]:m2[B] + 1]))
        pr.stdout.close()
        pr.wait()
        for x in (ref, q):
            try:
                os.remove(x)
            except OSError:
                pass
        if not cut:
            continue

        # ③ 유닛마다 변이 호출
        uq = "/tmp/jx_%s_%s.u.fa" % (smp, sl)
        with io.open(uq, "w", newline=NL) as g:
            for i, (_r, _u, sq) in enumerate(cut):
                g.write(">U%d%s%s%s" % (i, NL, sq, NL))
        rr = subprocess.run(["minimap2", "-x", "asm5", "-c", "--cs", "-t", TH,
                             "--secondary=no", CANON, uq],
                            capture_output=True, text=True)
        V = {}
        for L in rr.stdout.splitlines():
            fl = L.split("\t")
            if len(fl) < 12:
                continue
            i = int(fl[0][1:])
            if i in V:
                continue
            cs = ""
            for x in fl[12:]:
                if x.startswith("cs:Z:"):
                    cs = x[5:]
            tp = int(fl[7])
            out = set()
            for mm in CS.finditer(cs):
                op, v = mm.group(1), mm.group(2)
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
        try:
            os.remove(uq)
        except OSError:
            pass

        for i, (rid, u, _sq) in enumerate(cut):
            vv = V.get(i)
            if vv is None:
                T["변이못부름"] += 1
                continue
            o.write("%s\t%s\t%s\t%d\t%s\t1%s"
                    % (smp, (rid[:RIDLEN] if RIDLEN else rid), sl, u, (",".join(sorted(vv)) or "-"), NL))
            T["유닛"] += 1
        T["read"] += len(set(r for r, _u, _s in cut))
    sys.stderr.write("  %-10s read %d · 유닛 %d%s" % (smp, T["read"], T["유닛"], NL))
o.close()
sys.stderr.write("  합계 " + " · ".join("%s %d" % kv for kv in sorted(T.items())) + NL)
