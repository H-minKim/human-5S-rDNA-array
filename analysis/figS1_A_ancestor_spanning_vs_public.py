# -*- coding: utf-8 -*-
# 관통 read ↔ 공개조립본 염기자리 비교 (mkpairs3 의 앞 판)
#   단계 : S1 A 상류(v1)

# -*- coding: utf-8 -*-
# 관통 read 와 공개조립본이 염기 자리 수준에서 맞는가 (2026-09-05)
#
#  질문 (사용자)
#    "유닛 단위로 맞는 것 말고, 특징적인 변이 등의 위치가 잘 맞아?"
#
#  이미 아는 것
#    minimap2 전체 동일성 : 공개조립본 vs 우리 hifiasm = 100.0000% (중앙)
#                          공개조립본 vs 관통 read     =  99.4282% (중앙)
#    ⇒ 관통은 ~175 bp 마다 1 개씩 어긋난다. 그 어긋남의 성질을 봐야 한다.
#
#  방법
#    같은 하플로타입(유닛 ±1)으로 짝지어지는 129 쌍에 대해
#    ① 두 서열을 각각 유닛으로 분할하고 (5S 유전자 사슬, start-to-start 2,100-2,400)
#    ② 유닛마다 표준유닛(2,240 bp) 대비 변이 목록을 만들고
#    ③ 같은 유닛 번호끼리 변이 자리를 맞대본다.
#
#    유닛 번호는 양쪽 다 근위 junction 기준이라 위상이 같다
#    (실측: off_prox 중앙 1,516 bp · 616 하플로타입 중 98.9% 가 ±5).
#
#  치환과 indel 을 나눠서 본다.
#    ONT 의 주된 오류는 indel 이므로 섞어 세면 치환 일치율이 가려진다.
#  빈도로 층화한다.
#    한 하플로타입 안에서 대부분 유닛에 있는 변이(고빈도)는 진짜 계통 표지이고,
#    한두 유닛에만 있는 것은 오류일 수 있다. 둘을 나눠 봐야 답이 나온다.
import io
import os
import re
import sys
import subprocess
import collections

W = os.path.expanduser("~/psw")
os.makedirs(W, exist_ok=True)
D = os.path.expanduser("~/cmp4d")
R = os.path.expanduser("~/s5s_tuc/ref")
GENE = os.path.join(R, "5S_gene_query.fa")
CANON = os.path.join(R, "unit_canon.fa")
NL = "\n"
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(CMP)[::-1]


def units_of(p, tag):
    """서열을 유닛으로 분할한다 (방향 통일 포함)."""
    sq = "".join(x.strip() for x in
                 io.open(p, encoding="utf8").read().split(NL)[1:])
    for _ in range(2):
        fa = os.path.join(W, "g_%s.fa" % tag)
        db = os.path.join(W, "gd_%s" % tag)
        io.open(fa, "w", newline=NL).write(">x" + NL + sq + NL)
        subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
        r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                            "-outfmt", "6 sstart send", "-dust", "no",
                            "-word_size", "9", "-evalue", "1e-5",
                            "-max_target_seqs", "100000000", "-num_threads", TH],
                           capture_output=True, text=True)
        plus = minus = 0
        st = set()
        for L in r.stdout.splitlines():
            g = L.split("\t")
            if len(g) < 2:
                continue
            a, b = int(g[0]), int(g[1])
            if a <= b:
                plus += 1
            else:
                minus += 1
            st.add(min(a, b) - 1)
        for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
            try:
                os.remove(db + e)
            except OSError:
                pass
        os.remove(fa)
        if minus > plus:
            sq = rc(sq)
            continue
        break
    K = sorted(st)
    U = []
    for i in range(len(K) - 1):
        if 2100 <= K[i + 1] - K[i] <= 2400:
            U.append((len(U) + 1, sq[K[i]:K[i + 1]]))
    if K:
        last = min(K[-1] + 2240, len(sq))
        if last - K[-1] >= 2000:
            U.append((len(U) + 1, sq[K[-1]:last]))
    return U, (K[0] if K else -1)


def variants(U, tag):
    """유닛마다 표준유닛 대비 변이 — 치환과 indel 을 나눠서."""
    q = os.path.join(W, "u_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in U:
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
        sub, ind = set(), set()
        for m in re.finditer(r"([:*+-])([A-Za-z0-9]+)", cs):
            op, v = m.group(1), m.group(2)
            if op == ":":
                tp += int(v)
            elif op == "*":
                sub.add("%d%s" % (tp + 1, v[1].upper()))
                tp += 1
            elif op == "-":
                ind.add("%dD%d" % (tp + 1, len(v)))
                tp += len(v)
            elif op == "+":
                ind.add("%dI%d" % (tp + 1, len(v)))
        V[i] = (sub, ind)
    os.remove(q)
    return V


OUT = io.open(os.path.join(W, "ps_%s.tsv" % TAG), "w", encoding="utf8", newline=NL)
OUT.write("\t".join(["sample", "units", "n_unit_cmp", "class", "kind",
                     "n_pub", "n_span", "n_both", "recall", "precision"]) + NL)
SITE = io.open(os.path.join(W, "pssite_%s.tsv" % TAG), "w", encoding="utf8", newline=NL)
SITE.write("sample\tsite\tkind\tclass\tpub_units\tspan_units\tboth_units\n")

rows = [L.rstrip(NL).split("\t") for L in
        io.open(os.path.expanduser("~/pubspan_pairs.tsv"), encoding="utf8") if L.strip()]
mine = [r for i, r in enumerate(rows) if i % int(sys.argv[1]) == int(sys.argv[2])]

for n, r in enumerate(mine):
    smp, nu, fpub, fspan = r[0], int(r[1]), r[2], r[3]
    pp, ps = os.path.join(D, fpub), os.path.join(D, fspan)
    if not (os.path.exists(pp) and os.path.exists(ps)):
        continue
    tag = "%s_%s_%d" % (TAG, smp, n)
    UA, ka = units_of(pp, tag + "a")
    UB, kb = units_of(ps, tag + "b")
    if not UA or not UB:
        continue
    VA = variants(UA, tag + "a")
    VB = variants(UB, tag + "b")
    m = min(len(UA), len(UB))
    if m < 5:
        continue
    # 빈도 층화 — 그 하플로타입 안에서 몇 개 유닛에 나오는가
    cntA = collections.Counter()
    cntB = collections.Counter()
    for i in range(1, m + 1):
        for x in VA.get(i, (set(), set()))[0] | VA.get(i, (set(), set()))[1]:
            cntA[x] += 1
        for x in VB.get(i, (set(), set()))[0] | VB.get(i, (set(), set()))[1]:
            cntB[x] += 1
    for kind, gi in (("치환", 0), ("indel", 1)):
        for cls, lo, hi in (("고빈도(>=50% 유닛)", 0.5, 1.01),
                            ("중간(10-50%)", 0.1, 0.5),
                            ("저빈도(<10%)", 0.0, 0.1)):
            A = set()
            B = set()
            for i in range(1, m + 1):
                for x in VA.get(i, (set(), set()))[gi]:
                    if lo <= cntA[x] / float(m) < hi:
                        A.add((i, x))
                for x in VB.get(i, (set(), set()))[gi]:
                    if lo <= cntB[x] / float(m) < hi:
                        B.add((i, x))
            both = A & B
            OUT.write("\t".join([smp, str(nu), str(m), cls, kind,
                                 str(len(A)), str(len(B)), str(len(both)),
                                 "%.4f" % (len(both) / float(len(A))) if A else "-",
                                 "%.4f" % (len(both) / float(len(B))) if B else "-"]) + NL)
    # 자리 단위 — 그 자리가 두 서열에서 몇 개 유닛에 나오는가
    for x in set(cntA) | set(cntB):
        a = cntA.get(x, 0)
        b = cntB.get(x, 0)
        bo = sum(1 for i in range(1, m + 1)
                 if x in (VA.get(i, (set(), set()))[0] | VA.get(i, (set(), set()))[1])
                 and x in (VB.get(i, (set(), set()))[0] | VB.get(i, (set(), set()))[1]))
        fr = max(a, b) / float(m)
        cls = ("고빈도" if fr >= 0.5 else "중간" if fr >= 0.1 else "저빈도")
        SITE.write("%s\t%s\t%s\t%s\t%d\t%d\t%d%s" %
                   (smp, x, "indel" if ("D" in x or "I" in x) else "치환",
                    cls, a, b, bo, NL))
    OUT.flush()
    SITE.flush()
    sys.stderr.write("  %d/%d %s 유닛 %d\n" % (n + 1, len(mine), smp, m))
    sys.stderr.flush()
OUT.close()
SITE.close()
print("done")
