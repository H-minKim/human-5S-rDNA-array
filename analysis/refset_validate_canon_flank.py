# -*- coding: utf-8 -*-
# 정본 flank 검정 (조립본 flank ↔ 정본 양끝이 맞는가)
#   단계 : 참조세트 검정
#   포함 판단 : Fig S1 이 주장하는 참조세트 품질의 근거를 내는 검정이라 넣는다.

# -*- coding: utf-8 -*-
# flank read 만으로 만든 consensus 와 지도 서열을 대조한다 (2026-09-06)
#
#  묻는 것 (사용자)
#    "배정된 flank read 의 서열만으로 consensus map 을 만들어,
#     우리가 아는 read 로 만든 서열과 지도 서열을 비교하자"
#
#  왜 의미가 있나
#    지도 서열은 조립본에서 왔다 (공개조립본 64.3% · 우리 hifiasm 33.2% · 관통 2.5%).
#    read 는 독립 자료다. 둘이 맞으면 지도가 그 개체의 진짜 분자를 담고 있다는 뜻이고,
#    어긋나면 지도(=조립본)나 read 배정 둘 중 하나가 틀린 것이다.
#
#  방법
#    ① HAP_READ_FINAL.tsv 로 read 를 슬롯에 배정 (flank read 만, 이미 확정된 것)
#    ② 각 read 를 유닛으로 자르고 유닛별 변이집합을 낸다
#    ③ 사슬 배치로 read 가 지도의 몇 번 유닛부터인지 정한다
#       (chainplace.py 규칙: 점수 |S∩V|−|S−V| · 연속 K 유닛 · LOO k=8 에서 100%)
#    ④ 유닛마다 그 자리를 덮은 read 들의 변이를 모아 다수결 consensus 를 만든다
#    ⑤ 지도의 varlist 와 대조 — 재현율·정밀도
#
#  ONT read 는 ~99% 라 개별 변이엔 잡음이 있다. 그래서 빈도로 층화해서 본다
#    (앞서 실측: 고빈도 변이는 정밀도 99.3%, 저빈도는 23.5%).
#  배치가 불확실한 read 는 뺀다 (margin 하한).
import io
import os
import re
import sys
import gzip
import subprocess
import collections

W = os.path.expanduser("~/consflank")
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
NL = "\n"
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
K = int(os.environ.get("K", "8"))
MINMG = int(os.environ.get("MINMG", "1"))     # 배치 최소 마진
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(CMP)[::-1]


def score(S, V):
    return len(S & V) - len(S - V)


MAP = collections.defaultdict(dict)
for i, L in enumerate(io.open(os.path.expanduser("~/CONSENSUS_MAP_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 15 or f[4] == "PARTIAL_PROX":
        continue
    MAP[(f[0], f[1])][int(f[3])] = set() if f[14] == "-" else set(f[14].split(","))

ASG = collections.defaultdict(list)          # smp -> [(read, slot)]
seen = set()
for i, L in enumerate(io.open(os.path.expanduser("~/HAP_READ_FINAL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 7 or f[4] == "-":
        continue
    if (f[0], f[1]) in seen:
        continue
    seen.add((f[0], f[1]))
    ASG[f[0]].append((f[1], f[4]))


def find_reads(s):
    for p in ("~/s5s_ont/hits/%s.hits.fastq.gz" % s,
              "~/s5s_new/hits/%s.hits.fastq.gz" % s,
              "~/kgp_ont/out/%s.array.fastq.gz" % s):
        q = os.path.expanduser(p)
        if os.path.exists(q) and os.path.getsize(q) > 1000:
            return q
    return None


def units_of(seq, tag):
    fa = os.path.join(W, "r_%s.fa" % tag)
    io.open(fa, "w", newline=NL).write(">r" + NL + seq + NL)
    db = os.path.join(W, "rd_%s" % tag)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                        "-evalue", "1e-5", "-max_target_seqs", "100000000",
                        "-num_threads", TH], capture_output=True, text=True)
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
        return None, True
    Kp = sorted(st)
    U = []
    for i in range(len(Kp) - 1):
        if 2100 <= Kp[i + 1] - Kp[i] <= 2400:
            U.append(seq[Kp[i]:Kp[i + 1]])
    return U, False


def variants(units, tag):
    if not units:
        return []
    q = os.path.join(W, "u_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(units):
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
        for m in re.finditer(r"([:*+-])([A-Za-z0-9]+)", cs):
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
    return [V.get(i, set()) for i in range(len(units))]


def place(chain, smp, slot):
    M = MAP.get((smp, slot))
    if not M or len(chain) < K:
        return None
    keys = sorted(M)
    best = []
    for i0 in keys:
        if not all((i0 + t) in M for t in range(K)):
            continue
        n = min(len(chain), max(keys) - i0 + 1)
        if n < K:
            continue
        s = sum(score(M[i0 + t], chain[t]) for t in range(n) if (i0 + t) in M)
        best.append((s, i0))
    if len(best) < 2:
        return None
    best.sort(reverse=True)
    return best[0][1], best[0][0] - best[1][0]


OUT = io.open(os.path.join(W, "cf_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
smps = [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]
DONE = os.path.join(W, "done")
os.makedirs(DONE, exist_ok=True)

for smp in smps:
    mk = os.path.join(DONE, smp + ".done")
    if os.path.exists(mk):
        continue
    q = find_reads(smp)
    if not q or smp not in ASG:
        continue
    want = dict(ASG[smp])
    seqs = {}
    with gzip.open(q, "rt") as fh:
        j = 0
        rid = None
        for L in fh:
            m = j & 3
            j += 1
            if m == 0:
                # HAP_READ_FINAL.tsv 의 read ID 는 전체 UUID(36자)다.
                #   자르면 매칭이 통째로 실패한다 (실측: 배치 0 건).
                rid = L[1:].split()[0]
            elif m == 1 and rid in want:
                seqs[rid] = L.strip()
    # 유닛 자리마다 read 변이를 모은다
    COV = collections.defaultdict(list)      # (slot, unit) -> [변이집합]
    nplace = nskip = 0
    for rid, slot in sorted(want.items()):
        sq = seqs.get(rid)
        if not sq:
            continue
        tag = "%s_%s_%s" % (TAG, smp, rid[:8])
        U, isrc = units_of(sq, tag)
        if isrc:
            U, _ = units_of(rc(sq), tag)
        if not U or len(U) < K:
            nskip += 1
            continue
        ch = variants(U, tag)
        p = place(ch, smp, slot)
        if not p or p[1] < MINMG:
            nskip += 1
            continue
        i0, mg = p
        nplace += 1
        for t, v in enumerate(ch):
            COV[(slot, i0 + t)].append(v)
    # 유닛마다 다수결 consensus 를 만들어 지도와 대조
    for (slot, u), lst in sorted(COV.items()):
        M = MAP.get((smp, slot), {}).get(u)
        if M is None:
            continue
        n = len(lst)
        c = collections.Counter()
        for v in lst:
            for x in v:
                c[x] += 1
        for th, lab in ((0.5, "maj"), (0.8, "strict")):
            cons = set(x for x, k in c.items() if k / float(n) >= th)
            OUT.write("%s\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d%s" %
                      (smp, slot, u, lab, n, len(M), len(cons),
                       len(M & cons), len(cons - M), NL))
    OUT.flush()
    io.open(mk, "w").write("ok")
    sys.stderr.write("  %s 배치 %d · 제외 %d · 유닛자리 %d\n" %
                     (smp, nplace, nskip, len(COV)))
    sys.stderr.flush()
OUT.close()
print("done")
