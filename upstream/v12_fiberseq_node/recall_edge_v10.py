# -*- coding: utf-8 -*-
# read 유닛 가장자리 결실 규칙(B) — NO_FULL 유닛만 다시 부른다 (노드에서 실행 · 다시 자르거나 정렬하지 않는다)
#   read 유닛이 정본 1-2240 을 다 못 덮은 원인은 read 가장자리 1~3bp 결실 오류가 대부분 (실측)
#   규칙 B : 시작·끝 모자람이 각각 10bp 이하면 모자란 만큼 결실 토큰으로 적고 OK_EDGE 로 살린다
#             시작 모자람 s -> "1Ds" · 끝 모자람 e -> "(2241-e)De" (G8 의 끝 채움 규칙이 이미 이 토큰을 적는다)
#             10bp 를 넘으면 그대로 NO_FULL (버림)
#   유닛 서열 = 저장된 정렬 기록(SAM)의 SEQ 에서 저장된 read 좌표로 꺼낸다 (원 read 좌표 · 방향 반영)
#   검산 : 개체마다 OK 유닛 20 개를 같은 방법으로 다시 불러 저장된 변이목록과 같은지 대조 (서열 꺼내기가 맞는지)
#   인자 : JX|AO · 유닛표 · 출력표
import io
import os
import re
import sys
import random
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
KIND, INP, OUTP = sys.argv[1], sys.argv[2], sys.argv[3]
EDGE = 10
CANON = H("~/s5s_tuc/ref/unit_canon.fa")
CS = re.compile(r"([:*+\-=])([A-Za-z0-9]+)")
CANSEQ = "".join(x.strip() for x in io.open(CANON, encoding="utf8").read().split(NL)[1:]).upper()
random.seed(20260911)


def call(seqs, tag):
    q = "/tmp/re_%s.fa" % tag
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "4", "--secondary=no", CANON, q],
                       capture_output=True, text=True)
    OUT = {}
    for L in r.stdout.splitlines():
        f2 = L.split(TAB)
        if len(f2) < 12:
            continue
        i = int(f2[0][1:])
        if i in OUT:
            continue
        cs = ""
        for x in f2[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        qs, qe, ts, te = int(f2[2]), int(f2[3]), int(f2[7]), int(f2[8])
        toks, reb, tp = [], [], ts
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == "=":
                reb.append(v.upper())
                tp += len(v)
            elif op == ":":
                n = int(v)
                reb.append(CANSEQ[tp:tp + n])
                tp += n
            elif op == "*":
                toks.append("%d%s" % (tp + 1, v[1].upper()))
                reb.append(v[1].upper())
                tp += 1
            elif op == "-":
                toks.append("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                toks.append("%dI%d" % (tp + 1, len(v)))
                reb.append(v.upper())
        qseq = seqs[i]
        kt = 0
        while te + kt < 2240 and qe + kt < len(qseq):
            if qseq[qe + kt] != CANSEQ[te + kt]:
                toks.append("%d%s" % (te + kt + 1, qseq[qe + kt]))
            kt += 1
        if te + kt < 2240:
            toks.append("%dD%d" % (te + kt + 1, 2240 - (te + kt)))
        kh = 0
        while ts - 1 - kh >= 0 and qs - 1 - kh >= 0:
            if qseq[qs - 1 - kh] != CANSEQ[ts - 1 - kh]:
                toks.append("%d%s" % (ts - kh, qseq[qs - 1 - kh]))
            kh += 1
        OUT[i] = dict(toks=toks, qs=qs, qe=qe, reb="".join(reb), fill=kh + kt, cov_e=te + kt, cov_s=ts - kh)
    os.remove(q)
    return OUT


def judge(seq, res, edge):
    if res is None:
        return "NOALN", None
    if res["reb"] != seq[res["qs"]:res["qe"]]:
        return "SEQ_MISMATCH", None
    toks = list(res["toks"])
    sh, eh = res["cov_s"], 2240 - res["cov_e"]
    if sh == 0 and eh == 0:
        return ("OK_FILLED" if res["fill"] else "OK"), sorted(set(toks))
    if edge and sh <= EDGE and eh <= EDGE:
        if sh > 0:
            toks.append("1D%d" % sh)
        return "OK_EDGE", sorted(set(toks))
    return "NO_FULL", None


# 표 열 위치
if KIND == "JX":
    C_SMP, C_RID, C_JOB, C_STR, C_F0, C_F1, C_VD, C_VL = 0, 1, 3, 5, 7, 8, 11, 12
else:
    C_SMP, C_RID, C_JOB, C_STR, C_F0, C_F1, C_VD, C_VL = 0, 1, 4, 6, 9, 10, 13, 14
rows = [L.rstrip(NL).split(TAB) for L in io.open(INP, encoding="utf8") if L.strip()]
by = collections.defaultdict(list)
for i, r in enumerate(rows):
    by[(r[C_SMP], r[C_JOB])].append(i)


def samfile(smp, job):
    if KIND == "JX":
        return H("~/rjrepro/%s.placed.sam" % smp) if job == "JX" else H("~/rjrepro_c2/%s.placed.sam" % smp)
    return H("~/aog8v10/sam/%s.%s.sam" % (smp, job))


T = collections.Counter()
for (smp, job), idx in sorted(by.items()):
    need_nf = [i for i in idx if rows[i][C_VD] == "NO_FULL"]
    oks = [i for i in idx if rows[i][C_VD] in ("OK", "OK_FILLED")]
    chk = random.sample(oks, min(20, len(oks)))
    todo = need_nf + chk
    if not todo:
        continue
    want = set(rows[i][C_RID] for i in todo)
    SEQ = {}
    sp = samfile(smp, job)
    if not os.path.exists(sp):
        T["SAM 없음"] += len(todo)
        continue
    for L in io.open(sp, encoding="utf8"):
        if L.startswith("@"):
            continue
        f = L.split(TAB, 11)
        if f[0] in want and f[0] not in SEQ:
            SEQ[f[0]] = (f[9], int(f[1]) & 16)
    seqs = []
    for i in todo:
        r = rows[i]
        s, rev = SEQ[r[C_RID]]
        L0 = len(s)
        f0, f1 = int(r[C_F0]), int(r[C_F1])
        seqs.append(s[L0 - f1:L0 - f0] if rev else s[f0:f1])
    V = call(seqs, "%s_%s_%s" % (KIND, smp, job))
    for k, i in enumerate(todo):
        if k < len(need_nf):
            vd, toks = judge(seqs[k], V.get(k), True)
            T["NO_FULL -> " + vd] += 1
            rows[i][C_VD] = vd
            rows[i][C_VL] = (",".join(toks) or "-") if toks is not None else "-"
        else:
            vd, toks = judge(seqs[k], V.get(k), False)
            old = rows[i][C_VL]
            new = (",".join(toks) or "-") if toks is not None else "-"
            T["검산 OK 재호출 같음" if (vd == rows[i][C_VD] and new == old) else "검산 OK 재호출 다름"] += 1
with io.open(OUTP, "w", encoding="utf8", newline=NL) as o:
    for r in rows:
        o.write(TAB.join(r) + NL)
print("%s %s · %s" % (KIND, os.path.basename(INP), " · ".join("%s %d" % kv for kv in sorted(T.items()))))
