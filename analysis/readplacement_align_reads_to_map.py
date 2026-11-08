# -*- coding: utf-8 -*-
# consensus map 에 read 배치 → ra_*.tsv · fr_*.tsv
#   단계 : 배치

# -*- coding: utf-8 -*-
# B-2 v2. consensus map 에 read 배치 + flank read 직접 검출 (2026-09-05)
#
#  v1 대비 바뀐 것
#
#  ① 참조에 flank 50 kb 두 개를 함께 넣는다.
#     참조 = [하플로타입1][하플로타입2][PROX50k][DIST50k]
#     - flank 를 가진 read 는 flank 참조에 붙으므로 정렬로 직접 분류된다.
#       (옛 units_all.tsv 는 264 개체만 담고 있어 나머지가 통째로 빠졌다.
#        사용자 지적: "모든 개체는 flank 가 있는 read 가 있어" — 맞다.
#        표에 없는 것이지 read 가 없는 게 아니다.)
#     - flank 부분이 붙을 곳이 생기므로 정렬 비율을 그대로 해석할 수 있다
#       (v1 에서는 flank 가 clip 되어 정렬 비율이 31-89% 로 낮게 나왔다).
#
#  ② read 분류 (정렬 결과로만 판정)
#       SPAN      : 근위·원위 flank 를 둘 다 건드림 + array 도 건드림
#       SINGLE_P  : 근위 flank 만
#       SINGLE_D  : 원위 flank 만
#       ARRAY_ONLY: flank 없음
#
#  ③ 산출이 둘이다
#       ra_*.tsv : 유닛별 read 수·메틸 (지도용)
#       fr_*.tsv : read 별 배정·염기 (flank read 점검용)
#
#  두 하플로타입을 반드시 함께 넣는다 (한쪽만 넣으면 반대쪽 read 까지 붙는다).
#  방향은 CONSENSUS_SUM_G6.tsv 의 revcomp 플래그를 그대로 쓴다 (재계산 금지).
#  minimap2 -y 로 MM/ML 을 옮긴다 (빼면 beta 가 전부 "-" 가 된다).
#  MAPQ 를 함께 낸다 — 두 하플로타입이 ~99.6% 같아 못 가르면 MAPQ 0 이 된다.
import io
import os
import re
import sys
import gzip
import bisect
import subprocess
import collections
import pysam

W = os.path.expanduser("~/ralign")
os.makedirs(W, exist_ok=True)
D = os.path.expanduser("~/cmp4d")
R = os.path.expanduser("~/s5s_tuc/ref")
NL = "\n"
TH = os.environ.get("TH", "8")
TAG = os.environ.get("TAG", "x")
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
PROX = "PROX50k"
DIST = "DIST50k"
# flank 보유 판정 — 길이만으로 보면 안 된다.
#   실측: 길이 300 bp 만 요구했더니 COLO829-BL 의 read 4,229 중 4,041(96%)이
#   flank 보유로 잡혔다. flank 50 kb 안의 Alu·L1 같은 반복요소에 붙은 위양성이다.
#   ⇒ junction 까지의 거리를 함께 본다.
#     FLANK_proximal_50000 : junction 이 끝    (CHM13 227,694,885-227,744,884)
#     FLANK_distal_50000   : junction 이 시작  (CHM13 228,025,900-228,075,899)
#   진짜 flank read 는 array 와 이어져 있으므로 junction 쪽 끝에 붙어야 한다.
#   거리 경계는 실측으로 정했다 (시범 13 개체 · 정본 units_all.tsv 와 대조).
#     junction 거리 분포가 뚜렷한 이봉이다:
#       근위 0 에 327 개 · 꼬리 44 까지  ||  그다음이 168 (9 개) · 23,641-33,533
#       원위 0 에 307 개 · 꼬리 44 까지  ||  그다음이 5,652-7,503
#     정본이 SPAN 이라 한 18 개는 전부 양쪽 jx <= 17 이다.
#     168 짜리는 근위 flank 안 ~5,922 bp 반복블록에 붙은 것으로, 8 개가
#     길이·거리까지 똑같았다 (5,916-5,924 bp · jx 정확히 168) — 명백한 위양성.
#     ⇒ jx<=2000 이면 SPAN 을 27 개로 과대호출, jx<=50 이면 19 개이고
#       정본 SPAN 회수는 18/18 로 같다.
MINFL = 1000         # flank 참조에 붙은 최소 길이
JXWIN = 50           # junction 에서 이 거리 안에 닿아야 한다 (실측 이봉 경계)


def rc(s):
    return s.translate(CMP)[::-1]


def rdfa(p):
    return "".join(x.strip() for x in
                   io.open(p, encoding="utf8").read().split(NL)[1:])


FL_P = rdfa(os.path.join(R, "FLANK_proximal_50000.fa"))
FL_D = rdfa(os.path.join(R, "FLANK_distal_50000.fa"))


def find_reads(s):
    """read 자료원 — ONT fastq 를 먼저, 없으면 TUC 의 HiFi FASTA.

    둘은 같은 자료가 아니다. 반드시 구분해 기록한다:
      ONT   `*.hits.fastq.gz` / `*.array.fastq.gz`
            - 길다(중앙 수십~수백 kb) -> flank-to-flank 관통이 가능
            - 헤더에 MM/ML 이 있어 메틸을 낼 수 있다
      HIFI  `~/s5s_tuc/cand/<sample>.fa`  (TUC 파이프라인 산출)
            - PacBio HiFi, 중앙 18.5 kb (~8 유닛) -> 관통은 원리상 거의 불가
            - FASTA 라 MM/ML 이 없다 -> beta 는 전부 "-" 가 된다
            - 대신 염기 정확도가 높아 서열 대조에는 오히려 낫다
    이걸 놓쳐서 "24 개체 재추출에 10 시간" 이라고 잘못 말한 적이 있다.
      실제로 자료가 아예 없는 개체는 7 개뿐이다.
    """
    for p in ("~/s5s_ont/hits/%s.hits.fastq.gz" % s,
              "~/s5s_new/hits/%s.hits.fastq.gz" % s,
              "~/kgp_ont/out/%s.array.fastq.gz" % s,
              "~/giab_hits/%s.hits.fastq.gz" % s,
              "~/cancer5s/hits/%s.hits.fastq.gz" % s):
        q = os.path.expanduser(p)
        if os.path.exists(q) and os.path.getsize(q) > 1000:
            return q, "ONT"
    q = os.path.expanduser("~/s5s_tuc/cand/%s.fa" % s)
    if os.path.exists(q) and os.path.getsize(q) > 1000:
        return q, "HIFI"
    return None, None


def sniff_platform(path, guess):
    """플랫폼은 디렉터리가 아니라 read 이름으로 판별한다.

    실측 사고: `~/cancer5s/hits/*.fastq.gz` 를 디렉터리만 보고 ONT 로 판정해
    PacBio HiFi read(m84039_...)에 map-ont 를 썼다. COLO829 의 정렬 동일성이
    69% 로 떨어졌다 (정상 개체는 97-99%).
    PacBio 는 `m<숫자>_<날짜>_...` 또는 `/ccs` 를 이름에 갖는다.
    """
    try:
        op = gzip.open if path.endswith(".gz") else io.open
        with op(path, "rt") as fh:
            for L in fh:
                L = L.strip()
                if not L or L[0] not in "@>":
                    continue
                nm = L[1:].split()[0]
                if "|" in nm:
                    nm = nm.split("|", 1)[1]
                if re.match(r"^m\d+_\d+", nm) or "/ccs" in nm:
                    return "HIFI"
                return "ONT" if guess == "ONT" else guess
    except Exception:
        pass
    return guess


SEQF = {}
for L in io.open(os.path.expanduser("~/hap_seq_node.tsv"), encoding="utf8"):
    f = L.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
UNI = {}
for i, L in enumerate(io.open(os.path.expanduser("~/CONSENSUS_SUM_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
        UNI[(f[0], f[1])] = int(f[4])

UB = collections.defaultdict(list)
for i, L in enumerate(io.open(os.path.expanduser("~/CONSENSUS_MAP_G6.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip(NL).split("\t")
    if len(f) < 8 or f[4] == "PARTIAL_PROX":
        continue
    UB[(f[0], f[1])].append((int(f[3]), int(f[5]) - 1, int(f[6])))

OUT = io.open(os.path.join(W, "ra_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
FR = io.open(os.path.join(W, "fr_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
smps = [x.strip() for x in io.open(sys.argv[1], encoding="utf8") if x.strip()]
DONE = os.path.join(W, "done")
os.makedirs(DONE, exist_ok=True)

for smp in smps:
    mk = os.path.join(DONE, smp + ".done")
    if os.path.exists(mk):
        continue
    q, rsrc = find_reads(smp)
    if q:
        rsrc = sniff_platform(q, rsrc)      # 디렉터리가 아니라 read 이름으로
    slots = sorted(k[1] for k in SEQF if k[0] == smp)
    if not q or not slots:
        sys.stderr.write("  건너뜀 %s (read %s · 슬롯 %d)\n" % (smp, bool(q), len(slots)))
        continue
    ref = os.path.join(W, "%s.ref.fa" % smp)
    names = {}
    reflen = {}
    with io.open(ref, "w", newline=NL) as g:
        for sl in slots:
            p = os.path.join(D, SEQF[(smp, sl)])
            if not os.path.exists(p):
                continue
            sq = rdfa(p)
            if RCF.get((smp, sl)):
                sq = rc(sq)
            nm = "%s_h%s" % (smp, sl)
            names[nm] = sl
            reflen[nm] = len(sq)
            g.write(">" + nm + NL + sq + NL)
        # flank 참조 — read 의 flank 보유를 정렬로 직접 판정하기 위해서다
        g.write(">" + PROX + NL + FL_P + NL)
        g.write(">" + DIST + NL + FL_D + NL)
    if not names:
        continue
    bam = os.path.join(W, "%s.bam" % smp)
    # 플랫폼에 맞는 프리셋 — HiFi 에 map-ont 를 쓰면 정확도를 버린다
    preset = "map-ont" if rsrc == "ONT" else "map-hifi"
    p1 = subprocess.Popen(["minimap2", "-ax", preset, "-y", "-t", TH,
                           "--secondary=no", ref, q],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    p2 = subprocess.Popen(["samtools", "sort", "-@", "2", "-o", bam],
                          stdin=p1.stdout, stderr=subprocess.DEVNULL)
    p1.stdout.close()
    p2.wait()
    subprocess.run(["samtools", "index", bam], check=True, capture_output=True)
    bf = pysam.AlignmentFile(bam, "rb")

    # ---- ① read 별 분류 : 어느 참조를 건드렸나 ----
    HIT = collections.defaultdict(dict)      # rid -> refname -> 정렬길이 합
    JX = collections.defaultdict(dict)       # rid -> refname -> junction 까지 최소거리
    BEST = {}                                # rid -> 하플로타입 위 1차정렬
    NSUP = collections.Counter()
    LP, LD = len(FL_P), len(FL_D)
    for a in bf.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary:
            continue
        rid = a.query_name
        rn = a.reference_name
        HIT[rid][rn] = HIT[rid].get(rn, 0) + a.query_alignment_length
        if rn == PROX:                       # junction 이 끝
            d = LP - a.reference_end
        elif rn == DIST:                     # junction 이 시작
            d = a.reference_start
        else:
            d = None
        if d is not None:
            JX[rid][rn] = min(JX[rid].get(rn, 10 ** 9), max(0, d))
        if a.is_supplementary:
            NSUP[rid] += 1
        elif rn in names:
            BEST[rid] = a
    KLS = {}
    for rid, h in HIT.items():
        hp = (h.get(PROX, 0) >= MINFL and JX[rid].get(PROX, 10 ** 9) <= JXWIN)
        hd = (h.get(DIST, 0) >= MINFL and JX[rid].get(DIST, 10 ** 9) <= JXWIN)
        KLS[rid] = ("SPAN" if (hp and hd) else "SINGLE_P" if hp
                    else "SINGLE_D" if hd else "ARRAY_ONLY")
    for rid, a in BEST.items():
        cl = KLS.get(rid, "ARRAY_ONLY")
        nmt = a.get_tag("NM") if a.has_tag("NM") else -1
        al = a.query_alignment_length
        idt = 100.0 * (1.0 - nmt / float(al)) if (al and nmt >= 0) else -1
        rl = a.infer_read_length() or len(a.query_sequence or "")
        # 정렬 비율은 보조정렬까지 합쳐서 본다 (flank 조각이 따로 붙는다)
        tot_al = sum(HIT[rid].values())
        sl = names[a.reference_name]
        rcv = 100.0 * (a.reference_end - a.reference_start) / reflen[a.reference_name]
        jp = JX[rid].get(PROX, -1)
        jd = JX[rid].get(DIST, -1)
        FR.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%.1f\t%.1f\t%.1f\t%d\t%d"
                 "\t%d\t%d\t%d\t%d%s" %
                 (smp, rid[:12], rsrc, cl, sl, UNI.get((smp, sl), 0),
                  a.mapping_quality, idt,
                  100.0 * al / rl if rl else 0,             # 1차정렬 비율
                  100.0 * tot_al / rl if rl else 0,         # 보조 포함 정렬 비율
                  rcv, NSUP[rid], rl,
                  HIT[rid].get(PROX, 0), jp,                # 근위 flank 길이·junction거리
                  HIT[rid].get(DIST, 0), jd, NL))           # 원위 flank 길이·junction거리

    # ---- ② 유닛별 read·메틸 ----
    for nm, sl in names.items():
        ub = UB.get((smp, sl), [])
        if not ub:
            continue
        span = collections.Counter()
        over = collections.Counter()
        spq = collections.Counter()
        spf = collections.Counter()          # flank 보유 read 만
        msum = collections.defaultdict(lambda: [0.0, 0])
        asum = collections.defaultdict(lambda: [0.0, 0])
        qsum = collections.defaultdict(lambda: [0.0, 0])
        ntot = 0
        nq0 = 0
        for a in bf.fetch(nm):
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            ntot += 1
            hiq = (a.mapping_quality >= 20)
            hasfl = KLS.get(a.query_name, "ARRAY_ONLY") != "ARRAY_ONLY"
            if a.mapping_quality == 0:
                nq0 += 1
            lo, hi = a.reference_start, a.reference_end
            try:
                mb = a.modified_bases
            except Exception:
                mb = None
            rmod = []
            if mb:
                q2r = dict(a.get_aligned_pairs(matches_only=True))
                for key, lst in mb.items():
                    if len(key) < 3 or key[2] not in ("m", 27551):
                        continue
                    for qp, ql in lst:
                        rp = q2r.get(qp)
                        if rp is not None:
                            rmod.append((rp, 1.0 if ql >= 128 else 0.0))
                rmod.sort()
            # 유닛마다 rmod 를 전부 훑으면 O(유닛 x 메틸호출) 이라 개체당 5.3 분이 걸렸다
            #   (유닛 154 x 호출 1만 x read 3천 = 46 억 회).
            #   rmod 도 ub 도 좌표순이므로 이분 탐색으로 구간만 잘라 쓴다.
            RP = [x[0] for x in rmod]
            for pos, s0, e0 in ub:
                if hi <= s0 or lo >= e0:
                    continue
                full = (lo <= s0 and hi >= e0)
                over[pos] += 1
                if full:
                    span[pos] += 1
                    if hiq:
                        spq[pos] += 1
                    if hasfl:
                        spf[pos] += 1
                if not RP:
                    continue
                i0 = bisect.bisect_left(RP, s0)
                i1 = bisect.bisect_left(RP, e0)
                if i1 <= i0:
                    continue
                sv = 0.0
                for j in range(i0, i1):
                    sv += rmod[j][1]
                nv = i1 - i0
                asum[pos][0] += sv
                asum[pos][1] += nv
                if full:
                    msum[pos][0] += sv
                    msum[pos][1] += nv
                    if hiq:
                        qsum[pos][0] += sv
                        qsum[pos][1] += nv
        for pos, s0, e0 in ub:
            bb = msum[pos]
            aa = asum[pos]
            qq = qsum[pos]
            OUT.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d%s" %
                      (smp, sl, rsrc, pos, span[pos], over[pos] - span[pos], ntot,
                       ("%.4f" % (bb[0] / bb[1])) if bb[1] else "-",
                       ("%.4f" % (aa[0] / aa[1])) if aa[1] else "-",
                       bb[1], spq[pos],
                       ("%.4f" % (qq[0] / qq[1])) if qq[1] else "-",
                       qq[1], nq0, spf[pos], NL))
    bf.close()
    for x in (bam, bam + ".bai", ref):
        try:
            os.remove(x)
        except OSError:
            pass
    io.open(mk, "w").write("ok")
    OUT.flush()
    FR.flush()
    sys.stderr.write("  %s [%s] 완료 (read %d · flank 보유 %d · SPAN %d)\n" %
                     (smp, rsrc, len(HIT),
                      sum(1 for v in KLS.values() if v != "ARRAY_ONLY"),
                      sum(1 for v in KLS.values() if v == "SPAN")))
    sys.stderr.flush()
OUT.close()
FR.close()
print("done")
