# -*- coding: utf-8 -*-
# 그림이 쓰는 자를 만든 글 — g8canonv11f1/CONSENSUS_MAP_G8_full.tsv (bb68fe3a2a23) · CONSENSUS_SUM_G8_full.tsv (4c8252462824) · 616 슬롯 전체 재절단 · 앞서 예치했던 consmapG8.py (9960cb01b156) 는 두 세대 전이라 갈아끼웠다 (옛 판은 $BACKUP_ROOT/_old_versions/2026-10-24_v11_generation_fix/)
#   단계 : ④⑤ 방향 통일 · 유닛 자르기 · 변이

# -*- coding: utf-8 -*-
# G8v3 (2026-09-11) — G8v2 에서 유전자 경계 한 곳만 바꿨다
#   blastn 유전자 HSP 는 유전자 끝 염기가 변이면 한두 염기 앞에서 멈춘다 (616 정본 실측 57,878 HSP 중 50 개)
#   G8v2 는 HSP 의 정본 좌표만 써서 그 경계가 밀렸다 -> 질의 좌표(qstart/qend)로 유전자 전체 길이까지 늘려 경계를 정한다
#     + 가닥 : 시작 = sstart - (qstart - 1) · 끝 = send + (119 - qend)
#     - 가닥 : 시작 = send - (119 - qend)   · 끝 = sstart + (qstart - 1)
#   작업 폴더 ~/g8v3 · 그 밖(방향 · 자르기 규칙 · 변이 부르기 · 검산)은 G8v2 그대로
# G8v2 (2026-09-11) — METHOD_CONSENSUS_MAP_G8_REBUILD.md 대로 최신 정본(합본 v7 이 쓴 정본)에서 다시 자른다
#   바꾼 곳은 METHOD §5-5 표 네 곳 + 기록용 두 곳뿐. 자르기·변이·검산 방법은 G8 그대로.
#   (1) SEQF  : 환경변수 HSN 의 새 hap_seq_node (없으면 중단 — 옛 파일로 조용히 돌지 않게)
#   (2) .fa   : 새 폴더 CMPD (기본 ~/cmp4d_v7) — 옛 ~/cmp4d 는 그대로 둔다
#   (3) RCF   : 삭제 — 가닥 다수결만 (METHOD §4). SUM 의 revcomp 열 = 가닥 다수결로 뒤집었는가
#   (4) W     : ~/g8v2 — 옛 ~/g8 를 덮어쓰지 않는다
#   (5) slot_space 열 추가 (= assembly, METHOD §3-3) · 슬롯 번호는 옛 hap_seq_node 그대로
#   (6) 정본 없음 · 파일 없음 · 뒤집은 뒤에도 - 가닥 우세 를 조용히 넘기지 않고 센다
# consensus map G8 — 조립본 array 를 처음부터 다시 자른다 (2026-09-10)
#
#  G6 의 결함과 G7 의 임시조치
#    G6 : 유닛을 유전자 시작점에서 잘랐다 -> [유전자 119bp][AluY..NTS2]
#         정본유닛은 [AluY..NTS2][유전자] 이므로 119bp 회전 어긋남.
#         minimap2 는 감아 돌지 못해 정본 2122-2240 을 영영 안 덮었다.
#         결과: 지도 변이가 2101 에서 끊기고 5S 유전자 변이 0 건.
#    G7 : 잘린 유닛을 사후에 119bp 회전시켜 때웠다. 동작은 했으나
#         잘못 자른 것을 덧붙여 고치는 방식이라 정본이 될 수 없다.
#
#  G8 = 자를 때부터 정본 위상으로 자른다
#    정본유닛 = [AluY 1-318][NTS1 319-1178][GArich 1179-1551][NTS2 1552-2121][유전자 2122-2240]
#    ⇒ 정본은 유전자가 끝난 바로 다음부터 시작한다.
#    ⇒ 그러므로 유닛 경계를 유전자 HSP 의 end 로 잡는다.
#         유닛 i = (gene_end[i] + 1) .. gene_end[i+1]
#       이러면 회전이 필요 없다. 자른 유닛이 곧 정본 위상이다.
#
#  flank 는 이미 처리돼 있다 (공유)
#    cmp4d/*.fa 는 asm_units3.py 가 50kb flank 탐침으로 잘라낸 flank-to-flank 서열.
#    .fa 의 경계가 곧 junction 이다. 여기서 flank 를 다시 다루지 않는다.
#
#  근위 1번(PARTIAL_PROX)은 빠진다
#    junction 에서 첫 유전자 끝까지는 온전한 유닛이 아니다 (실측 중앙 1,516 bp).
#    array 마지막의 남는 조각도 마찬가지다. 둘 다 PARTIAL 로 표시하고
#    변이를 부르지 않는다.
#
#  검산은 서열로 한다
#    ① 변이 목록으로 유닛 서열을 복원해 원본과 글자 단위로 대조
#    ② 정렬기가 안 덮은 말단은 서열로 직접 비교해 채운다
#       (minimap2 는 말단이 불일치면 안 늘린다 -> 그 변이를 잘라먹는다.
#        실측 400 건 중 237 건이 진짜 변이였다)
#    ③ 유닛 경계는 유전자 HSP 로 확인한다 — bp/2240 환산 금지
import io
import os
import re
import sys
import subprocess
import collections

NL = "\n"
W = os.path.expanduser(os.environ.get("WDIR", "~/g8v4"))
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
PRESET = os.environ.get("PRESET", "asm20")
LO = int(os.environ.get("LO", "2100"))          # 유닛 간격 허용 하한
HI = int(os.environ.get("HI", "2400"))          # 상한
FRAG = int(os.environ.get("FRAG", "5"))         # 확장량이 이보다 크면 조각 HSP 로 본다
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CS = re.compile(r"([:*+\-=])([A-Za-z0-9]+)")
NUM = re.compile(r"\d+")


def rc(s):
    return s.translate(CMP)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


CANSEQ = rd(CANON).upper()
PAD = int(os.environ.get("PAD", "100"))           # 유닛 자를 때 앞뒤 여유
CTX = int(os.environ.get("CTX", "200"))           # 정본유닛 양옆 문맥 (탠덤이라 같은 유닛의 끝·앞)
REFSEQ = CANSEQ[-CTX:] + CANSEQ + CANSEQ[:CTX]
CANREF = os.path.join(W, "canon_ctx.fa")
io.open(CANREF, "w", newline=NL).write(">canon_ctx" + NL + REFSEQ + NL)
if len(CANSEQ) != 2240:
    sys.stderr.write("정본유닛 길이가 2240 이 아니다: %d%s" % (len(CANSEQ), NL))
    raise SystemExit(1)
GENESEQ = rd(GENE).upper()
if len(GENESEQ) != 119:
    sys.stderr.write("유전자 질의 길이가 119 가 아니다: %d%s" % (len(GENESEQ), NL))
    raise SystemExit(1)
# 정본 안에서 유전자가 2122..2240 인지 서열로 확인한다 (가정하지 않는다)
if CANSEQ[2121:2240] != GENESEQ:
    sys.stderr.write("정본 2122-2240 이 유전자와 다르다 — 중단%s" % NL)
    raise SystemExit(1)

SEQF = {}
for x in io.open(os.path.expanduser(os.environ["HSN"]), encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
for i, x in enumerate([]):   # RCF 사용 안 함 (METHOD §4)
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
FLIP = set()


def gene_hits(seq, tag):
    """유전자 HSP 의 (시작, 끝) 과 우세 가닥. 좌표는 1 기반 포함."""
    fa = os.path.join(W, "g_%s.fa" % tag)
    db = os.path.join(W, "gd_%s" % tag)
    io.open(fa, "w", newline=NL).write(">x" + NL + seq + NL)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send length pident qstart qend", "-dust", "no",
                        "-word_size", "9", "-evalue", "1e-5",
                        "-max_target_seqs", "100000000", "-num_threads", TH],
                       capture_output=True, text=True)
    plus = minus = 0
    hits = []
    for L in r.stdout.splitlines():
        g = L.split("\t")
        if len(g) < 6:
            continue
        a, b = int(g[0]), int(g[1])
        qs_, qe_ = int(g[4]), int(g[5])
        glen = len(GENESEQ)
        if a <= b:
            plus += 1
            s_, e_ = a - (qs_ - 1), b + (glen - qe_)
        else:
            minus += 1
            s_, e_ = b - (glen - qe_), a + (qs_ - 1)
        hits.append((max(1, s_), min(len(seq), e_), int(g[2]),
                     (qs_ - 1) + (glen - qe_)))   # 네 번째 = 확장량
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    # F 방식 (2026-09-13 감사) — 1kb 바구니를 구간 겹침 + 조각 판정으로 바꿨다
    #   바구니 결함 : 유령 끝 145074 와 진짜 끝 145147 이 73bp 차이인데
    #                 바구니 144/145 로 갈려 둘 다 살아남아 가짜 끊김 + 경계 73bp 밀린 유닛
    #   확장 규칙 자체는 유지한다 (끝 염기 변이로 한두 염기 모자란 진짜 끝이 있다)
    keep = []
    for h in sorted(hits):
        if keep and h[0] <= keep[-1][1]:          # ① 구간 겹침 -> 긴 것만
            if h[2] > keep[-1][2]:
                keep[-1] = h
        else:
            keep.append(h)
    keep.sort(key=lambda x: x[1])
    i = 0
    while i + 1 < len(keep):                      # ② 간격 < LO 인 이웃 : 조각 쪽만 버린다
        if keep[i + 1][1] - keep[i][1] < LO:
            fa_, fb_ = keep[i][3] > FRAG, keep[i + 1][3] > FRAG
            if fa_ and not fb_:
                keep.pop(i)
                i = max(0, i - 1)
                continue
            if fb_ and not fa_:
                keep.pop(i + 1)
                i = max(0, i - 1)
                continue
        i += 1
    return [(a_, b_, c_) for (a_, b_, c_, d_) in keep], (minus > plus)


def call_variants(seqs, tag):
    q = os.path.join(W, "q_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", PRESET, "-c", "--cs=long", "-t", TH,
                        "--secondary=no", CANREF, q], capture_output=True, text=True)
    OUT = {}
    for L in r.stdout.splitlines():
        f2 = L.split("\t")
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
        toks = []
        reb = []
        tp = ts
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == "=":
                reb.append(v.upper())
                tp += len(v)
            elif op == ":":
                n = int(v)
                reb.append(REFSEQ[tp:tp + n])
                tp += n
            elif op == "*":
                if CTX <= tp < CTX + 2240:
                    toks.append("%d%s" % (tp - CTX + 1, v[1].upper()))
                reb.append(v[1].upper())
                tp += 1
            elif op == "-":
                if CTX <= tp < CTX + 2240:
                    toks.append("%dD%d" % (tp - CTX + 1, len(v)))
                tp += len(v)
            elif op == "+":
                if CTX <= tp < CTX + 2240:
                    toks.append("%dI%d" % (tp - CTX + 1, len(v)))
                reb.append(v.upper())
        # 말단은 정렬기를 믿지 않고 서열로 직접 채운다
        # 채움 없음 (G8v4) : 문맥 붙인 참조에 여유째 정렬했으므로 정렬 결과만 쓴다
        OUT[i] = dict(toks=sorted(set(toks)), qs=qs, qe=qe, ts=ts, te=te,
                      reb="".join(reb), fill=0, cov_e=(2240 if te >= CTX + 2240 else 0), cov_s=(0 if ts <= CTX else 1))
    os.remove(q)
    return OUT


def verify(seq, res):
    if res is None:
        return "NOALN", "정렬 없음"
    if res["reb"] != seq[res["qs"]:res["qe"]]:
        return "SEQ_MISMATCH", "복원 %d vs 원본 %d" % (len(res["reb"]), res["qe"] - res["qs"])
    if res["cov_e"] < 2240 or res["cov_s"] > 0:
        return "NO_FULL", "정본 %d..%d (+채움 %d)" % (res["ts"], res["te"], res["fill"])
    return ("OK_FILLED" if res["fill"] else "OK"), (
        "말단 %d bp 서열로 채움" % res["fill"] if res["fill"] else "")


HDR = ["sample", "slot", "N", "array_pos", "seg_type", "unit_start", "unit_end",
       "unit_len", "gap_to_next", "chain_ok", "n_var", "n_var_gene",
       "varlist", "verify", "verify_note", "slot_space"]
MAP = io.open(os.path.join(W, "CONSENSUS_MAP_G8_%s.tsv" % TAG), "w",
              encoding="utf8", newline=NL)
MAP.write("\t".join(HDR) + NL)
SUM = io.open(os.path.join(W, "CONSENSUS_SUM_G8_%s.tsv" % TAG), "w",
              encoding="utf8", newline=NL)
SUM.write("sample\tslot\tbp\tN_units\trevcomp\toff_prox\toff_dist\tn_break\n")
STAT = collections.Counter()

for line in io.open(sys.argv[1], encoding="utf8"):
    f = line.strip().split("\t")
    if len(f) < 2:
        continue
    smp, sl = f[0], f[1]
    fn = SEQF.get((smp, sl))
    if not fn:
        STAT["정본 없음"] += 1
        continue
    p = os.path.join(os.path.expanduser(os.environ.get("CMPD", "~/cmp4d_v7")), fn)
    if not os.path.exists(p):
        STAT["파일 없음"] += 1
        sys.stderr.write("  파일 없음 %s_%s %s%s" % (smp, sl, p, NL))
        continue
    a = rd(p).upper()
    if RCF.get((smp, sl)):
        a = rc(a)
    hits, isrc = gene_hits(a, "%s_%s_%s" % (TAG, smp, sl))
    if isrc:
        # 방향이 뒤집혀 있으면 서열을 뒤집고 다시 찾는다 (좌표를 손으로 뒤집지 않는다)
        a = rc(a)
        FLIP.add((smp, sl))
        hits, isrc = gene_hits(a, "%s_%s_%s_r" % (TAG, smp, sl))
        if isrc:
            STAT["뒤집은 뒤에도 - 가닥 우세"] += 1
    if len(hits) < 3:
        STAT["유전자 HSP 부족"] += 1
        continue
    # 유닛 경계 = 유전자 끝 -> 다음 유전자 끝  (정본 위상)
    ends = [e for (s, e, ln) in hits]
    bounds = []
    nbreak = 0
    for i in range(len(ends) - 1):
        g = ends[i + 1] - ends[i]
        if LO <= g <= HI:
            bounds.append((ends[i] + 1, ends[i + 1]))     # 1 기반 포함
        else:
            nbreak += 1
    if not bounds:
        STAT["사슬 없음"] += 1
        continue
    N = len(bounds)
    seqs = [a[max(0, s - 1 - PAD):min(len(a), e + PAD)] for (s, e) in bounds]   # 앞뒤 여유
    V = call_variants(seqs, "%s_%s_%s" % (TAG, smp, sl))
    off_prox = ends[0]                       # junction ~ 첫 유전자 끝
    off_dist = len(a) - ends[-1]             # 마지막 유전자 끝 ~ junction
    SUM.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n"
              % (smp, sl, len(a), N, "Y" if (smp, sl) in FLIP else "N",
                 off_prox, off_dist, nbreak))
    # 근위 조각 (junction ~ 첫 유전자 끝) — 변이 부르지 않음
    MAP.write("\t".join([smp, sl, str(N), "0", "PARTIAL_PROX", "1", str(ends[0]),
                         str(ends[0]), "0", "-", "0", "0", "-", "SKIP", "온전유닛 아님", "assembly"]) + NL)
    STAT["PARTIAL_PROX"] += 1
    for i, (s, e) in enumerate(bounds):
        res = V.get(i)
        vd, note = verify(seqs[i], res)
        STAT[vd] += 1
        gap = (bounds[i + 1][0] - s) if i + 1 < N else 0
        if vd in ("OK", "OK_FILLED"):
            toks = res["toks"]
            ng = sum(1 for t in toks if int(NUM.match(t).group()) > 2121)
            STAT["유전자변이"] += ng
            vl = ",".join(toks) if toks else "-"
            MAP.write("\t".join([smp, sl, str(N), str(i + 1), "UNIT", str(s), str(e),
                                 str(e - s + 1), str(gap), "Y", str(len(toks)),
                                 str(ng), vl, vd, note, "assembly"]) + NL)
        else:
            MAP.write("\t".join([smp, sl, str(N), str(i + 1), "UNIT", str(s), str(e),
                                 str(e - s + 1), str(gap), "Y", "0", "0", "-",
                                 vd, note, "assembly"]) + NL)
    # 원위 조각 (마지막 유전자 끝 ~ junction)
    if off_dist > 0:
        MAP.write("\t".join([smp, sl, str(N), str(N + 1), "PARTIAL_DIST",
                             str(ends[-1] + 1), str(len(a)), str(off_dist), "0", "-",
                             "0", "0", "-", "SKIP", "온전유닛 아님", "assembly"]) + NL)
        STAT["PARTIAL_DIST"] += 1
    sys.stderr.write("  %-10s 슬롯%s 유닛 %d · 끊김 %d%s" % (smp, sl, N, nbreak, NL))
    sys.stderr.flush()
MAP.close()
SUM.close()
print("consensus map G8")
for k, v in sorted(STAT.items()):
    print("  %-18s %d" % (k, v))
