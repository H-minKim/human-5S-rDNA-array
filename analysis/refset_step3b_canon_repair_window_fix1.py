# -*- coding: utf-8 -*-
# 수술 창 (앵커 직접 지정판)
#   단계 : ③b 정본 수리 (수정1 · 2 자리)
# Reason: no machine-specific path is left in the code. A default would leave such a path on this line, so none is set.
#       to run : export NODE_HOME=/home/<user>

# -*- coding: utf-8 -*-
# [v11-수정1 사본 2026-09-16 : 지도 g8canonv10 · 정본 cmp4d_v11 · 앵커 직접 지정 · 입력 ~/fix1/pol · 출력 ~/fix1/splice. 그 밖 불변]
# 수술식 교체 v2 — 앵커를 위치로 고정한다 (2026-09-13)
#
#   v1 의 실패 : 앵커 유닛을 합의 전체에서 '최고 정합률' 로 찾았더니 탠덤 안 다른 사본에 붙었다
#                (NA20282 : 합의 3유닛인데 뒤 앵커가 가운데 유닛에 붙어 진짜 유닛 하나를 삼킴 -> 2유닛)
#   고침      : 합의는 정의상 [앞앵커 ... 뒤앵커] 구간이다. 앞앵커는 합의의 맨 앞 HSP,
#               뒤앵커는 맨 뒤 HSP 로 고정한다. 그리고 '합의에서 앵커 밖으로 남는 bp' 를 찍어
#               정말 끝에 붙었는지 확인한다(0 에 가까워야 한다).
#   검정      : ① 남는 bp ② 앵커 정합률(조립본 기준 100) ③ 끊김 ④ 유닛 수 = 합의와 같아야
import io
import os
import subprocess
import collections

H = os.environ["NODE_HOME"]   # the compute node's home directory
NL = chr(10)
TAB = chr(9)
W = os.path.join(H, "fix1", "splice")
os.makedirs(W, exist_ok=True)
GENE = os.path.join(H, "s5s_tuc", "ref", "5S_gene_query.fa")
CMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
LO, HI = 2100, 2400


def rc(s):
    return "".join(CMP.get(c, "N") for c in reversed(s))


def rdfa(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]).upper()


GLEN = len(rdfa(GENE))
SEQF = {}
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split(TAB)
    if len(f) >= 3:
        SEQF[(f[0], f[1])] = f[2]
FLIP = {}
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_SUM_G8_full.tsv"), encoding="utf8")
f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB)
    FLIP[(r[0], r[1])] = (r[4] == "Y")
U = collections.defaultdict(list)
f = io.open(os.path.join(H, "g8canonv10", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
for L in f:
    r = L.rstrip(NL).split(TAB)
    if r[ix["seg_type"]] == "UNIT":
        U[(r[0], r[1])].append((int(r[ix["array_pos"]]), int(r[ix["unit_start"]]), int(r[ix["unit_end"]])))
for k in U:
    U[k].sort()


def canon(smp, sl):
    a = rdfa(os.path.join(H, "cmp4d_v11", SEQF[(smp, sl)]))
    return rc(a) if FLIP.get((smp, sl)) else a


def hsps(q, t, tag):
    """덮음 >= 98% 인 HSP 전부 (pident, sstart, send)"""
    qf = os.path.join(W, "q_%s.fa" % tag)
    io.open(qf, "w", newline=NL).write(">q" + NL + q + NL)
    tf = os.path.join(W, "t_%s.fa" % tag)
    io.open(tf, "w", newline=NL).write(">t" + NL + t + NL)
    subprocess.run(["makeblastdb", "-in", tf, "-dbtype", "nucl", "-out", tf + ".db"], capture_output=True)
    r = subprocess.run(["blastn", "-task", "megablast", "-query", qf, "-db", tf + ".db",
                        "-outfmt", "6 pident length qstart qend sstart send", "-dust", "no",
                        "-evalue", "1e-20", "-max_target_seqs", "200"], capture_output=True, text=True)
    out = []
    for L in r.stdout.splitlines():
        g = L.split(TAB)
        if len(g) < 6:
            continue
        qs, qe = int(g[2]), int(g[3])
        if 100.0 * (qe - qs + 1) / len(q) < 98:
            continue
        out.append((float(g[0]), min(int(g[4]), int(g[5])), max(int(g[4]), int(g[5]))))
    return out


def gene_ends(seq, tag):
    fa = os.path.join(W, "g_%s.fa" % tag)
    db = os.path.join(W, "gd_%s" % tag)
    io.open(fa, "w", newline=NL).write(">x" + NL + seq + NL)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db], capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send length pident qstart qend", "-dust", "no",
                        "-word_size", "9", "-evalue", "1e-5", "-max_target_seqs", "100000000"],
                       capture_output=True, text=True)
    hits = []
    for L in r.stdout.splitlines():
        g = L.split(TAB)
        if len(g) < 6:
            continue
        a, b, qs, qe = int(g[0]), int(g[1]), int(g[4]), int(g[5])
        if a <= b:
            s_, e_ = a - (qs - 1), b + (GLEN - qe)
        else:
            s_, e_ = b - (GLEN - qe), a + (qs - 1)
        hits.append((max(1, s_), min(len(seq), e_), int(g[2]), (qs - 1) + (GLEN - qe)))
    keep = []
    for h in sorted(hits):
        if keep and h[0] <= keep[-1][1]:
            if h[2] > keep[-1][2]:
                keep[-1] = h
        else:
            keep.append(h)
    keep.sort(key=lambda x: x[1])
    i = 0
    while i + 1 < len(keep):
        if keep[i + 1][1] - keep[i][1] < LO:
            fa_, fb_ = keep[i][3] > 5, keep[i + 1][3] > 5
            if fa_ and not fb_:
                keep.pop(i)
                i = max(0, i - 1)
                continue
            if fb_ and not fa_:
                keep.pop(i + 1)
                i = max(0, i - 1)
                continue
        i += 1
    return [x[1] for x in keep]


def brk(ends):
    return sum(1 for i in range(len(ends) - 1) if not (LO <= ends[i + 1] - ends[i] <= HI))


SITES = [("HG00133", "2", 119, 129), ("HG03270", "1", 40, 47)]
OUT = io.open(os.path.join(W, "SPLICE_MID2_FIX1.tsv"), "w", encoding="utf8", newline=NL)
OUT.write(TAB.join(["sample", "slot", "after_unit", "asm_bp", "cons_bp", "splice_bp",
                    "left_over", "right_over", "mid_bp", "pid5", "pid3",
                    "genes_asm", "genes_cons", "genes_splice",
                    "brk_asm", "brk_cons", "brk_splice", "verdict"]) + NL)
print("%-9s %-3s %-8s %-6s %-6s %-6s %-9s %-7s %-11s %s"
      % ("개체", "h", "앵커", "조립본", "합의", "수술식", "남는bp", "가운데", "유전자수", "끊김(조립/합의/수술)"))
for (smp, sl, ap, apR) in SITES:
    A = canon(smp, sl)
    v = U[(smp, sl)]
    i = [j for j, x in enumerate(v) if x[0] == ap][0]
    iR = [j for j, x in enumerate(v) if x[0] == apR][0]
    a0, a1 = v[i][1], v[iR][2]
    u1 = A[v[i][1] - 1:v[i][2]]
    u2 = A[v[iR][1] - 1:v[iR][2]]
    polp = os.path.join(H, "fix1", "pol", "pol_%s_%s_%d.fa" % (smp, sl, ap))
    cons = rdfa(polp)
    h1 = hsps(u1, cons, "a1_%s%s%d" % (smp, sl, ap))
    h2 = hsps(u2, cons, "a2_%s%s%d" % (smp, sl, ap))
    if not h1 or not h2:
        print("%-9s %-3s 앵커가 합의에 안 붙음 — 건너뜀" % (smp, sl))
        continue
    b1 = min(h1, key=lambda x: x[1])          # 맨 앞 HSP = 앞 앵커
    b2 = max(h2, key=lambda x: x[2])          # 맨 뒤 HSP = 뒤 앵커
    left_over = b1[1] - 1                     # 합의에서 앞 앵커 앞에 남는 bp
    right_over = len(cons) - b2[2]            # 뒤 앵커 뒤에 남는 bp
    mid = cons[b1[2]:b2[1] - 1] if b2[1] - 1 > b1[2] else ""
    new_win = u1 + mid + u2
    old_win = A[a0 - 1:a1]
    ga = gene_ends(old_win, "o_%s%s%d" % (smp, sl, ap))
    gc = gene_ends(cons, "c_%s%s%d" % (smp, sl, ap))
    gn = gene_ends(new_win, "n_%s%s%d" % (smp, sl, ap))
    ba, bc, bn = brk(ga), brk(gc), brk(gn)
    same = (len(gn) == len(gc))
    vd = ("OK" if (bn == 0 and same) else
          ("유닛수 다름" if not same else "끊김 남음"))
    print("%-9s %-3s %-8s %-6d %-6d %-6d %-9s %-7d %-11s %s"
          % (smp, sl, "%d/%d" % (v[i][0], v[iR][0]), len(old_win), len(cons), len(new_win),
             "%d/%d" % (left_over, right_over), len(mid),
             "%d/%d/%d" % (len(ga), len(gc), len(gn)),
             "%d/%d/%d  %s" % (ba, bc, bn, vd)))
    OUT.write(TAB.join([smp, sl, str(ap), str(len(old_win)), str(len(cons)), str(len(new_win)),
                        str(left_over), str(right_over), str(len(mid)),
                        "%.2f" % b1[0], "%.2f" % b2[0],
                        str(len(ga)), str(len(gc)), str(len(gn)),
                        str(ba), str(bc), str(bn), vd]) + NL)
    io.open(os.path.join(W, "spl_%s_%s_%d.fa" % (smp, sl, ap)), "w", newline=NL).write(
        ">spl|%s|h%s|after_unit%d|anchor%d-%d%s%s%s" % (smp, sl, ap, a0, a1, NL, new_win, NL))
OUT.close()
