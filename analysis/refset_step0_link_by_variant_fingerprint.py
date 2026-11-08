# -*- coding: utf-8 -*-
# 남은 72 개체를 변이 지문 사슬로 잇는다
#   단계 : ① 하플로타입 확정

# -*- coding: utf-8 -*-
# 남은 72 개체의 prox <-> dist 연결 — 변이 지문 사슬 (2026-09-06)
#
#  문제는 작다.
#    72 개체도 side 별 하플로타입 분리는 이미 끝나 있다 (hap_reads_v5.tsv,
#    네 군 prox1/prox2/dist1/dist2 가 72/72 전부 존재, 군 read 수 균형 중앙 0.80).
#    없는 것은 오직 prox 무리와 dist 무리를 잇는 것 —
#    개체마다 (prox1-dist1) 인가 (prox1-dist2) 인가, 두 경우 중 하나를 고르는 문제다.
#    read 를 다시 배치하는 문제가 아니다.
#
#  방법은 이 프로젝트가 이미 검증한 것을 그대로 쓴다 (chainplace.py / loo2.py)
#    점수 = |S∩V| − |S−V|   ("일치−결측")
#    S = 지도의 그 유닛 변이집합 · V = read 의 그 유닛 변이집합
#    연속 K 유닛 사슬. LOO 실측: k=1 49% · k=2 70% · k=4 94% · k=8 100%
#    minimap2 는 쓰지 않는다 — 실측상 하플로타입 배정 도구가 아니다
#      (관통 read 로 잰 정확도 MAPQ 1-19 에서 42.1%, 근위 flank read 의 44%가 MAPQ<20)
#
#  지도 varlist 를 지문으로 쓴다 (chainplace.py 와 동일).
#    맵은 여러 read 의 합의라 '공통 변이' 역할을 한다.
#
#  공통 지문은 정보를 버린다(사용자 지적). 그래서 2 단계로 간다:
#    1 단계 — 군의 공통 변이(그 군 read 의 다수가 갖는 것)로 판정
#    2 단계 — 1 단계에서 margin 이 안 나오면 read 개별 점수를 다 합산해 판정
#             (개별 정보를 버리지 않는다)
import io
import os
import sys
import gzip
import subprocess
import collections

W = os.path.expanduser("~/linkhap")
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
NL = "\n"
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
K = int(os.environ.get("K", "8"))
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(CMP)[::-1]


def score(S, V):
    return len(S & V) - len(S - V)


# ---- 지도: (개체,슬롯) -> {유닛번호: 변이집합} ----
MAP = collections.defaultdict(dict)
for i, L in enumerate(io.open(os.path.expanduser("~/CONSENSUS_MAP_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 15 or f[4] == "PARTIAL_PROX":
        continue
    vs = set()
    if f[14] != "-":
        for v in f[14].split(","):
            vs.add(v)
    MAP[(f[0], f[1])][int(f[3])] = vs

# ---- read 의 side 라벨 ----
SIDE = collections.defaultdict(dict)
for i, L in enumerate(io.open(os.path.expanduser("~/hap_reads_v5.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 4:
        continue
    p = f[2].split("|")
    rid = p[1] if len(p) > 1 else f[2]
    SIDE[f[0]][(rid, f[1])] = f[3]


def find_reads(s):
    for p in ("~/s5s_ont/hits/%s.hits.fastq.gz" % s,
              "~/s5s_new/hits/%s.hits.fastq.gz" % s,
              "~/kgp_ont/out/%s.array.fastq.gz" % s):
        q = os.path.expanduser(p)
        if os.path.exists(q) and os.path.getsize(q) > 1000:
            return q
    return None


def unit_chain(seq, tag):
    """read 를 유닛으로 자르고 유닛마다 표준유닛 대비 변이집합을 낸다."""
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
    import re
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
    """사슬을 그 슬롯 지도에 붙여 (최고점, 마진) 반환."""
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
    return best[0][0], best[0][0] - best[1][0]


OUT = io.open(os.path.join(W, "link_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
smps = [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]
DONE = os.path.join(W, "done")
os.makedirs(DONE, exist_ok=True)

for smp in smps:
    mk = os.path.join(DONE, smp + ".done")
    if os.path.exists(mk):
        continue
    q = find_reads(smp)
    slots = sorted(s for (a, s) in MAP if a == smp)
    if not q or len(slots) < 2:
        sys.stderr.write("  건너뜀 %s (read %s · 슬롯 %d)\n" % (smp, bool(q), len(slots)))
        continue
    want = set(r for (r, sd) in SIDE.get(smp, {}))
    seqs = {}
    with gzip.open(q, "rt") as fh:
        j = 0
        rid = None
        for L in fh:
            m = j & 3
            j += 1
            if m == 0:
                rid = L[1:].split()[0]
            elif m == 1 and rid in want:
                seqs[rid] = L.strip()
    # 군별 점수 합
    SC = collections.defaultdict(lambda: collections.Counter())
    NR = collections.Counter()
    for (rid, sd), h in sorted(SIDE.get(smp, {}).items()):
        sq = seqs.get(rid)
        if not sq:
            continue
        tag = "%s_%s_%s" % (TAG, smp, rid[:8])
        U, isrc = unit_chain(sq, tag)
        if isrc:
            U, _ = unit_chain(rc(sq), tag)
        if not U or len(U) < K:
            continue
        ch = variants(U, tag)
        NR[(sd, h)] += 1
        for slot in slots:
            p = place(ch, smp, slot)
            if p:
                SC[(sd, h)][slot] += p[0]
    # 군 -> 슬롯
    call = {}
    for k in SC:
        c = SC[k]
        if len(c) < 2:
            continue
        a = c.most_common()
        call[k] = (a[0][0], a[0][1] - a[1][1], NR[k])
    for (sd, h), (slot, mg, nr) in sorted(call.items()):
        OUT.write("%s\t%s\t%s\t%s\t%d\t%d\n" % (smp, sd, h, slot, mg, nr))
    OUT.flush()
    io.open(mk, "w").write("ok")
    sys.stderr.write("  %s 군 %d 판정 (read %d)\n" % (smp, len(call), sum(NR.values())))
    sys.stderr.flush()
OUT.close()
print("done")
