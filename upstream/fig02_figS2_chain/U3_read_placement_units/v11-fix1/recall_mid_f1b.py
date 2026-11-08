# -*- coding: utf-8 -*-
# handling read units short by more than 10bp (run on the node)
#   "if it is not at the end of the read and there is a unit next to it (a unit inside the array), mark it simply as a deletion even if it is 10bp or more"
#   (1) read middle unit (not the minimum or maximum of the unit numbers covered by that read, with units before and after it in coordinates)
#        -> write the shortfall as a deletion token and mark OK_MIDDEL  (start shortfall s -> "1Ds", end shortfall e -> "(2241-e)De" , same tokens as rule B)
#   (2) units at either end of the read -> NO_FULL_READEND (discarded, as now)
#   Also measured (user question "what is the state of the neighboring unit on that side?"): for each neighboring unit on the side short by more than 10bp,
#        neighbor check status, read gap between the two units (overlap is negative), extra bases the neighbor has beyond the canonical unit on the facing side, and the neighbor's own shortfall
#   unit sequence = taken from the stored SAM SEQ using the stored read coordinates (same as recall_edge.py, nothing re-cut or re-aligned)
#   Arguments: JX|AO, input table (.edge.tsv), output table (.mid.tsv)   ->  output table plus output_table.nb.tsv (neighbor unit diagnostics)
import io
import os
import re
import sys
import subprocess
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
KIND, INP, OUTP = sys.argv[1], sys.argv[2], sys.argv[3]
EDGE = 10
CANON = H("~/s5s_tuc/ref/unit_canon.fa")
CS = re.compile(r"([-:*+=])([A-Za-z0-9]+)")
CANSEQ = "".join(x.strip() for x in io.open(CANON, encoding="utf8").read().split(NL)[1:]).upper()


def call(seqs, tag):
    q = "/tmp/rm_%s.fa" % tag
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
        OUT[i] = dict(toks=toks, qs=qs, qe=qe, reb="".join(reb), sh=ts - kh, eh=2240 - (te + kt),
                      xs=qs - kh, xe=len(qseq) - (qe + kt))
    os.remove(q)
    return OUT


def sz(n):
    return "11-50" if n <= 50 else "51-200" if n <= 200 else "201-1000" if n <= 1000 else ">1000"


if KIND == "JX":
    C_SMP, C_RID, C_JOB, C_POS, C_F0, C_F1, C_VD, C_VL = 0, 1, 3, 6, 7, 8, 11, 12
else:
    C_SMP, C_RID, C_JOB, C_POS, C_F0, C_F1, C_VD, C_VL = 0, 1, 4, 7, 9, 10, 13, 14
rows = [L.rstrip(NL).split(TAB) for L in io.open(INP, encoding="utf8") if L.strip()]
NOTE = ["-"] * len(rows)
byread = collections.defaultdict(list)
need = collections.defaultdict(set)
for i, r in enumerate(rows):
    byread[(r[C_SMP], r[C_JOB], r[C_RID])].append(i)
    if r[C_VD] == "NO_FULL":
        need[(r[C_SMP], r[C_JOB])].add(r[C_RID])


def samfile(smp, job):
    if KIND == "JX":
        return H("~/rjrepro/%s.placed.sam" % smp) if job == "JX" else H("~/rjrepro_c2/%s.placed.sam" % smp)
    return H("~/fix1/aog8f1b/sam/%s.%s.sam" % (smp, job))


T = collections.Counter()
NB = io.open(OUTP + ".nb.tsv", "w", encoding="utf8", newline=NL)
NB.write(TAB.join(["sample", "read", "job", "array_pos", "position", "side", "short_bp", "nb_array_pos", "nb_verify",
                   "gap_bp", "nb_extra_facing", "nb_short_facing", "nb_class"]) + NL)
for (smp, job), rids in sorted(need.items()):
    SEQ = {}
    for L in io.open(samfile(smp, job), encoding="utf8"):
        if L.startswith("@"):
            continue
        f = L.split(TAB, 11)
        if f[0] in rids and f[0] not in SEQ:
            SEQ[f[0]] = (f[9], int(f[1]) & 16)
    plan, todo, seqs, crd = [], {}, [], {}
    for rid in sorted(rids):
        s, rev = SEQ[rid]
        L0 = len(s)
        U = []
        for i in byread[(smp, job, rid)]:
            f0, f1 = int(rows[i][C_F0]), int(rows[i][C_F1])
            s0, s1 = (L0 - f1, L0 - f0) if rev else (f0, f1)
            U.append((s0, s1, i))
            crd[i] = (s0, s1)
        U.sort()
        ps = [int(rows[i][C_POS]) for _, _, i in U]
        lo, hi = min(ps), max(ps)
        for k, (s0, s1, i) in enumerate(U):
            if rows[i][C_VD] != "NO_FULL":
                continue
            prv = U[k - 1][2] if k > 0 else None
            nxt = U[k + 1][2] if k + 1 < len(U) else None
            mid_idx = lo < int(rows[i][C_POS]) < hi
            mid_crd = prv is not None and nxt is not None
            if mid_idx != mid_crd:
                T["중간판정 번호/좌표 엇갈림"] += 1
            plan.append((i, prv, nxt, mid_idx and mid_crd))
            for j in (i, prv, nxt):
                if j is not None and j not in todo:
                    todo[j] = len(seqs)
                    a, b = crd[j]
                    seqs.append(s[a:b])
    V = call(seqs, "%s_%s_%s" % (KIND, smp, job))
    for i, prv, nxt, mid in plan:
        r = rows[i]
        res = V.get(todo[i])
        sq = seqs[todo[i]]
        if res is None or res["reb"] != sq[res["qs"]:res["qe"]]:
            T["재호출 실패(그대로 NO_FULL)"] += 1
            continue
        sh, eh = res["sh"], res["eh"]
        pos = "MID" if mid else "READEND"
        for side, short, nb in (("start", sh, prv), ("end", eh, nxt)):
            if short <= EDGE:
                continue
            T["%s · %s 모자람 >10bp" % (pos, side)] += 1
            T["%s · %s 모자람 크기 %s" % (pos, side, sz(short))] += 1
            if nb is None:
                cls, row = "옆 유닛 없음(read 끝)", ["-", "-", "-", "-", "-"]
            else:
                rn = V.get(todo[nb])
                if side == "start":
                    gap = crd[i][0] - crd[nb][1]
                    ex, nsh = (rn["xe"], rn["eh"]) if rn else (None, None)
                else:
                    gap = crd[nb][0] - crd[i][1]
                    ex, nsh = (rn["xs"], rn["sh"]) if rn else (None, None)
                if rn is None:
                    cls = "옆 유닛 정렬 안 됨"
                elif nsh > EDGE:
                    cls = "옆 유닛도 마주보는 쪽이 모자람(>10bp)"
                elif ex >= max(5, short // 2):
                    cls = "옆 유닛이 정본 밖 염기를 더 가짐(경계 밀림)"
                elif gap >= max(5, short // 2):
                    cls = "두 유닛 사이 read 염기가 비어 있음(경계 틈)"
                else:
                    cls = "옆 유닛 정상 · 경계 붙어 있음(read 에 그 서열 없음=결실)"
                row = [rows[nb][C_POS], rows[nb][C_VD], str(gap), str(ex), str(nsh)]
                T["%s · %s 옆 유닛 상태 %s" % (pos, side, rows[nb][C_VD])] += 1
            T["%s · %s → %s" % (pos, side, cls)] += 1
            NB.write(TAB.join([smp, r[C_RID], job, r[C_POS], pos, side, str(short)] + row + [cls]) + NL)
        NOTE[i] = "short_start=%d;short_end=%d;pos=%s" % (sh, eh, pos)
        if mid:
            toks = list(res["toks"])
            if sh > 0:
                toks.append("1D%d" % sh)
            r[C_VD] = "OK_MIDDEL"
            r[C_VL] = ",".join(sorted(set(toks))) or "-"
            T["NO_FULL -> OK_MIDDEL (결실로 표시)"] += 1
        else:
            r[C_VD] = "NO_FULL_READEND"
            T["NO_FULL -> NO_FULL_READEND (버림)"] += 1
NB.close()
with io.open(OUTP, "w", encoding="utf8", newline=NL) as o:
    for i, r in enumerate(rows):
        o.write(TAB.join(r + [NOTE[i]]) + NL)
for k, v in sorted(T.items()):
    print("CNT%s%s%s%s%s%d" % (TAB, KIND, TAB, k, TAB, v))
