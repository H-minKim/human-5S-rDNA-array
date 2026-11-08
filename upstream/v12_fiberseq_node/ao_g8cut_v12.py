# -*- coding: utf-8 -*-
# ARRAY_ONLY 배치 read 를 G8 틀로 자른다 (2026-09-11) — 로컬 WSL 에서 실행
#   배치는 다시 하지 않는다 : read · 슬롯 · 유닛 번호는 저장된 v7 배치(AOJV) 그대로
# ① source reads : fastq taken from the compute node (header MM/ML, sequence and quality unchanged)
#   ② 정렬 : 자기 슬롯 정본에서 저장된 유닛 구간(정본 안으로 줄임) +-250bp 창 (minimap2 -a -x map-ont --secondary=no -y)
#            (2026-09-11 고침) 옛 창 +-3 유닛은 비슷한 유닛이 이어진 배열에서 read 가 1~3 유닛 미끄러짐(ARRAY_ONLY 11%)
#            저장 배치 대조 : read 단위로 스페이서 변이를 저장 번호 g+d 와 비교해 뚜렷이 어긋나면(같은 번호<0.6 · 옆 번호>=0.8)
#            저장 위치로 옮긴 read 전체 발자국 창(+-300bp)에 다시 정렬하고, 저장 번호와 맞으면 채택 (READS 11열에 기록)
#            -y 로 fastq 머리줄의 MM/ML 이 SAM 에 그대로 따라온다 · FASTQ 로 넣으므로 QUAL 도 그대로
#   ③ 경계 : 정본(G8v3) 유닛 경계를 정렬로 read 에 옮긴다. 유전자 끝(T) 뒤에 종결신호 TTTT 가 있어
#            경계가 T 연속 한가운데에 놓인다 -> T 연속 바깥의 T 아닌 염기 둘을 기준점으로 삼고,
#            read 의 T 연속에서 정본과 같은 수(보통 1 개)의 T 를 앞 유닛 끝에, 나머지를 다음 유닛 머리에 준다.
#            T 연속이 아닌 경계는 그대로 옮기고, 대응 염기가 없으면 가장 가까운 대응 염기(표시 NEAREST)
#   ④ 변이 : G8v3 지도와 같은 방식 (asm20 · --cs=long · 말단 서열 채움 · 복원 대조)
#   ⑤ 번호 : 정렬로 얻은 정본 좌표 -> G6 번호(저장값과 대조) -> G8 번호. 스페이서 변이로 저장 번호와 맞는지 확인
#   인자 : 작업폴더(aog8) · 개체 목록 · 병렬 수
import io
import os
import re
import sys
import gzip
import glob
import bisect
import subprocess
import collections
import multiprocessing

NL = chr(10)
TAB = chr(9)
D = sys.argv[1]
LST = sys.argv[2]
NP = int(sys.argv[3]) if len(sys.argv) > 3 else 6
CMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CS = re.compile(r"([:*+\-=])([A-Za-z0-9]+)")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
NUM = re.compile(r"\d+")
WIN = 3            # (옛 값 · 쓰지 않음) 2026-09-11 창 +-3 유닛에서 read 가 1~3 유닛 미끄러짐(ARRAY_ONLY 8%) 확인 -> PADBP 로 바꿈
PADBP = 250        # 창 = 저장 G6 유닛 범위(최소~최대) 정본 좌표 +- 250bp : 유닛 하나(2,240bp)만큼 미끄러질 여지를 없앤다
SHIFTMAX = 3       # 저장 배치와 read 단위 번호 차이를 볼 범위 (+-유닛)


def rc(s):
    return s.translate(CMP)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


REF = os.path.join(D, "s5s_tuc", "ref")
CANON = os.path.join(REF, "unit_canon.fa")
CANSEQ = rd(CANON).upper()
assert len(CANSEQ) == 2240


def samples(p):
    return set(L.split(TAB)[0] for L in io.open(p, encoding="utf8") if L.strip())


def loadmap(p, types=("UNIT",)):
    U = collections.defaultdict(list)
    f = io.open(p, encoding="utf8")
    ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    for L in f:
        r = L.rstrip(NL).split(TAB)
        if r[ix["seg_type"]] in types:
            U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]), int(r[ix["array_pos"]])))
    for k in U:
        U[k].sort()
    return U


def call_variants(seqs, tag, wd):
    q = os.path.join(wd, "q_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "1", "--secondary=no", CANON, q],
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
        OUT[i] = dict(toks=sorted(set(toks)), qs=qs, qe=qe, reb="".join(reb), fill=kh + kt, cov_e=te + kt, cov_s=ts - kh)
    os.remove(q)
    return OUT


def verify(seq, res):
    if res is None:
        return "NOALN"
    if res["reb"] != seq[res["qs"]:res["qe"]]:
        return "SEQ_MISMATCH"
    if res["cov_e"] < 2240 or res["cov_s"] > 0:
        return "NO_FULL"
    return "OK_FILLED" if res["fill"] else "OK"


def boundary(A, m2, b, side):
    """정본 A(0 기반 문자열)의 유전자 끝 b(1 기반) 뒤 경계를 read 로 옮긴다.
    side='end' 이면 앞 유닛의 read 끝 좌표(포함), 'start' 면 다음 유닛의 read 시작 좌표. (read 좌표 0 기반, SAM SEQ 기준)"""
    i = b - 1                                   # 0 기반 유전자 끝
    if 0 <= i < len(A) and A[i] == "T":
        rs = i
        while rs - 1 >= 0 and A[rs - 1] == "T":
            rs -= 1
        re_ = i
        while re_ + 1 < len(A) and A[re_ + 1] == "T":
            re_ += 1
        a0, a1 = rs - 1, re_ + 1                # T 연속 바깥의 기준점 (정본 0 기반)
        k1 = i - rs + 1                         # 정본에서 앞 유닛에 속하는 T 수
        if a0 in m2 and a1 in m2:
            q0, q1 = m2[a0], m2[a1]
            lr = q1 - q0 - 1                    # read 의 T 연속 길이
            qe = q0 + min(k1, max(lr, 0))       # 앞 유닛 끝 (포함)
            return (qe if side == "end" else qe + 1), "RUN"
    # T 연속이 아니거나 기준점이 안 잡힐 때 : 그대로 옮기고, 안 되면 가장 가까운 대응 염기
    j = i if side == "end" else i + 1
    if j in m2:
        return m2[j], "DIRECT"
    for d in range(1, 50):
        k = j - d if side == "end" else j + d
        if k in m2:
            return m2[k] + (0 if side == "end" else 0), "NEAREST"
    return None, "FAIL"


def align_cut(smp, slot, rid, r16, cfg, lab, rec, A, u6, ws, we, wd):
    """창 [ws, we] (정본 1 기반) 에 read 를 정렬하고 G8 유닛으로 자른다.
    돌려주는 값 : (SAM 한 줄(정본 좌표), [(유닛서열, meta)], (정렬 정본시작0, 정렬 정본끝0(미포함), 왼쪽 clip, 오른쪽 clip)) 또는 None"""
    wfa = os.path.join(wd, "w.fa")
    rfq = os.path.join(wd, "r.fq")
    io.open(wfa, "w", newline=NL).write(">%s_h%s%s%s%s" % (smp, slot, NL, A[ws - 1:we], NL))
    io.open(rfq, "w", newline=NL).write(NL.join(rec) + NL)
    r = subprocess.run(["minimap2", "-a", "-x", "map-ont", "--secondary=no", "-y", "-t", "1", wfa, rfq],
                       capture_output=True, text=True)
    prim = None
    for L in r.stdout.splitlines():
        if L.startswith("@"):
            continue
        f = L.split(TAB)
        if f[2] != "*" and not (int(f[1]) & 0x900):
            prim = f
            break
    if prim is None:
        return None
    # 정본 좌표로 옮긴 SAM (POS += 창 시작 - 1) — 원 read 서열 · 품질 · MM/ML 은 손대지 않는다
    prim2 = list(prim)
    prim2[3] = str(int(prim[3]) + ws - 1)
    rpos = int(prim[3]) - 1 + ws - 1      # 정본 0 기반
    ref0 = rpos
    q2 = 0
    m2 = {}
    cig = [(int(m.group(1)), m.group(2)) for m in CIG.finditer(prim[5])]
    for k, op in cig:
        if op == "S":
            q2 += k
        elif op in "M=X":
            for t in range(k):
                m2[rpos + t] = q2 + t
            rpos += k
            q2 += k
        elif op == "I":
            q2 += k
        elif op in "DN":
            rpos += k
    lclip = cig[0][0] if cig and cig[0][1] in "SH" else 0
    rclip = cig[-1][0] if cig and cig[-1][1] in "SH" else 0
    seq = prim[9]
    flag = int(prim[1])
    L0 = len(seq)
    lo_r, hi_r = min(m2), max(m2)
    units = []
    for s8, e8, p8 in U8[(smp, slot)]:
        # 유닛 p8 = [s8, e8] : 시작 경계 = 앞 유닛 유전자 끝(s8-1) 뒤 · 끝 경계 = 이 유닛 유전자 끝(e8) 뒤
        if not (lo_r < s8 - 2 and e8 < hi_r):
            continue
        qs_, bs = boundary(A, m2, s8 - 1, "start")
        qe_, be = boundary(A, m2, e8, "end")
        if qs_ is None or qe_ is None or qe_ < qs_:
            continue
        useq = seq[qs_:qe_ + 1]
        # 원 read(fastq) 기준 좌표
        if flag & 16:
            f0, f1 = L0 - (qe_ + 1), L0 - qs_
        else:
            f0, f1 = qs_, qe_ + 1
        mid = s8 + 1000                     # 이 유닛 스페이서 한가운데 (정본 1 기반)
        i6 = bisect.bisect_right([x[0] for x in u6], mid) - 1
        g6 = u6[i6][2] if (i6 >= 0 and u6[i6][0] <= mid <= u6[i6][1]) else None
        units.append((useq, [smp, rid, r16, slot, cfg, lab, "-" if not (flag & 16) else "+rc", p8, g6, f0, f1, bs, be]))
    return TAB.join(prim2), units, (ref0, rpos, lclip, rclip)


def spacer(toks):
    return set(t for t in toks if int(NUM.match(t).group()) <= 2121)


def jac(a, b):
    un = a | b
    return len(a & b) / float(len(un)) if un else 1.0


def read_offset(items, P):
    """items = [(g6, 스페이서 변이 집합)] (변이를 부른 유닛만) -> {d: 저장 번호 g6+d 와의 평균 자카드}"""
    S = {}
    for d in range(-SHIFTMAX, SHIFTMAX + 1):
        v = [jac(t, spacer(P[g + d])) for g, t in items if g is not None and g + d in P]
        if len(v) >= max(2, len(items) // 2):
            S[d] = sum(v) / len(v)
    return S


def shift_pos(u6, pos1, d):
    """정본 좌표 pos1(1 기반)을 G6 유닛 d 개만큼 옮긴 좌표 (유닛 안 상대 위치 유지)"""
    st = [x[0] for x in u6]
    i = bisect.bisect_right(st, pos1) - 1
    i = max(i, 0)
    j = i + d
    if j < 0 or j >= len(u6):
        return None
    return pos1 + (u6[j][0] - u6[i][0])


def relabel(res, bd):
    """다시 정렬로도 저장 위치에 안 붙는 read : 자른 read 좌표·변이는 그대로, 유닛 번호만 저장 차이 bd 만큼 옮긴다.
    SAM 기록은 실제 정렬 그대로 두고 ZR 태그로 표시 (그 read 의 POS 는 저장 번호 체계보다 bd 유닛 어긋남)"""
    sam, units, foot = res
    units2 = [(us, mt[:7] + [mt[7] + bd, (mt[8] + bd) if mt[8] is not None else None] + mt[9:]) for us, mt in units]
    return (sam + TAB + "ZR:Z:relabeled_d=%+d" % bd, units2, foot)


PBCFG = ("AOJV_PB",)       # 이 cfg 는 read 이름을 자르지 않는다 (PacBio movie 이름은 앞 16 자가 겹친다)
def rkey(rid, cfg):
    return rid if cfg in PBCFG else rid[:16]


def work(arg):
    smp, cfg = arg                      # v12 : 개체마다 cfg 가 여럿일 수 있다
    wd = os.path.join(D, "tmp", "%s.%s" % (smp, cfg))
    os.makedirs(wd, exist_ok=True)
    sv = SAVED[(cfg, smp)]
    files = sorted(glob.glob(os.path.join(D, "aoreads", "%s.%s.*.fastq.gz" % (smp, cfg))))
    out_u, out_r, out_s = [], [], []
    FQ = {}
    for fp in files:
        lab = fp.split(".")[-3]
        with gzip.open(fp, "rt", encoding="utf8") as fh:
            rec = []
            for x in fh:
                rec.append(x.rstrip(NL))
                if len(rec) == 4:
                    rid = rec[0][1:].split()[0]
                    FQ.setdefault(rkey(rid, cfg), (rid, rec, lab))   # v12
                    rec = []
    RES = {}                  # r16 -> (sam, units, foot)
    NOTE = {}                 # r16 -> 다시 정렬 기록
    # ── 1차 : 창 = 저장 G6 유닛 범위(정본 안으로 줄임) +- PADBP
    for r16, (slot, P) in sorted(sv.items()):
        if r16 not in FQ:
            out_r.append([smp, r16, "-", slot, cfg, "-", "NOT_IN_FASTQ"])
            continue
        rid, rec, lab = FQ[r16]
        A = ARR[(smp, slot)]
        u6 = U6[(smp, slot)]
        p6 = dict((p, (s, e)) for s, e, p in u6)
        inr = sorted(p for p in p6 if min(P) <= p <= max(P))      # 저장 범위 중 정본에 있는 유닛 (배열 끝 너머 저장 번호는 잘림)
        if not inr:
            out_r.append([smp, r16, rid, slot, cfg, lab, "NO_WINDOW"])
            continue
        if max(P) > inr[-1] or min(P) < inr[0]:
            NOTE[r16] = "saved_beyond_canon(%d-%d>%d-%d)" % (min(P), max(P), inr[0], inr[-1])
        ws, we = max(1, p6[inr[0]][0] - PADBP), min(len(A), p6[inr[-1]][1] + PADBP)
        x = align_cut(smp, slot, rid, r16, cfg, lab, rec, A, u6, ws, we, wd)
        if x is None:
            out_r.append([smp, r16, rid, slot, cfg, lab, "NOALN"])
            continue
        RES[r16] = x
    order = [(r16, i) for r16 in sorted(RES) for i in range(len(RES[r16][1]))]
    V = call_variants([RES[r][1][i][0] for r, i in order], smp, wd) if order else {}
    CALL = {}                 # (r16, i) -> (vd, toks)
    for k, (r16, i) in enumerate(order):
        res = V.get(k)
        vd = verify(RES[r16][1][i][0], res)
        CALL[(r16, i)] = (vd, set(res["toks"]) if (res and vd in ("OK", "OK_FILLED")) else None)
    # ── 2차 : 저장 배치와 read 단위 번호가 뚜렷이 어긋나면(같은 번호 < 0.6 · 옆 번호 >= 0.8)
    #          저장 위치로 옮긴 'read 전체 발자국' 창(양끝 +-300bp)에 다시 정렬 — 유닛 하나만큼 미끄러지면 read 2kb 가 남아 불리
    REDO = {}
    RELAB = {}
    for r16, (sam, units, foot) in RES.items():
        slot, P = sv[r16]
        items = [(units[i][1][8], spacer(CALL[(r16, i)][1])) for i in range(len(units)) if CALL[(r16, i)][1] is not None]
        S = read_offset(items, P)
        if 0 not in S:
            continue
        bd = max(S, key=lambda d: (S[d], -abs(d)))
        if bd == 0 or not (S[bd] >= 0.8 and S[0] < 0.6):
            continue
        u6 = U6[(smp, slot)]
        A = ARR[(smp, slot)]
        ref0, ref1, lcl, rcl = foot
        fs, fe = shift_pos(u6, ref0 + 1, bd), shift_pos(u6, ref1, bd)
        if fs is None or fe is None:
            NOTE[r16] = (NOTE.get(r16, "") + ";" if r16 in NOTE else "") + "shift_d=%+d_no_room" % bd
            RELAB[r16] = bd
            continue
        ws2, we2 = max(1, fs - lcl - 300), min(len(A), fe + rcl + 300)
        rid, rec, lab = FQ[r16]
        y = align_cut(smp, slot, rid, r16, cfg, lab, rec, A, u6, ws2, we2, wd)
        if y is None or not y[1]:
            NOTE[r16] = (NOTE.get(r16, "") + ";" if r16 in NOTE else "") + "shift_d=%+d_realign_noaln" % bd
            RELAB[r16] = bd
            continue
        REDO[r16] = (y, bd, S[0], S[bd])
    if REDO:
        order2 = [(r16, i) for r16 in sorted(REDO) for i in range(len(REDO[r16][0][1]))]
        V2 = call_variants([REDO[r][0][1][i][0] for r, i in order2], smp + "_re", wd)
        CALL2 = {}
        for k, (r16, i) in enumerate(order2):
            res = V2.get(k)
            vd = verify(REDO[r16][0][1][i][0], res)
            CALL2[(r16, i)] = (vd, set(res["toks"]) if (res and vd in ("OK", "OK_FILLED")) else None)
        for r16, (y, bd, s0, sb) in REDO.items():
            slot, P = sv[r16]
            units = y[1]
            items = [(units[i][1][8], spacer(CALL2[(r16, i)][1])) for i in range(len(units)) if CALL2[(r16, i)][1] is not None]
            S2 = read_offset(items, P)
            if 0 in S2 and S2[0] >= 0.7 and S2[0] > s0 + 0.1 and max(S2, key=lambda d: (S2[d], -abs(d))) == 0:
                RES[r16] = y
                for i in range(len(units)):
                    CALL[(r16, i)] = CALL2[(r16, i)]
                NOTE[r16] = (NOTE.get(r16, "") + ";" if r16 in NOTE else "") + "realigned_d=%+d(%.3f->%.3f)" % (bd, s0, S2[0])
            else:
                NOTE[r16] = (NOTE.get(r16, "") + ";" if r16 in NOTE else "") + "shift_d=%+d_realign_failed(%.3f,%.3f->%s)" % (bd, s0, sb, "%.3f" % S2[0] if 0 in S2 else "-")
                RELAB[r16] = bd
    for r16, bd in RELAB.items():
        RES[r16] = relabel(RES[r16], bd)
        NOTE[r16] = NOTE[r16] + ";relabeled_d=%+d" % bd
    # ── 출력
    per = collections.defaultdict(list)
    for r16 in sorted(RES):
        sam, units, foot = RES[r16]
        out_s.append(sam)
        slot, P = sv[r16]
        for i, (useq, mt) in enumerate(units):
            vd, toks = CALL[(r16, i)]
            sc = "-"
            if toks is not None and mt[8] in P:
                sc = "%.3f" % jac(spacer(toks), spacer(P[mt[8]]))
                per[r16].append(float(sc))
            out_u.append(mt + [vd, (",".join(sorted(toks)) or "-") if toks is not None else "-", sc,
                               "Y" if mt[8] in P else "N"])
    for r16, (slot, P) in sorted(sv.items()):
        if r16 in FQ and not any(x[1] == r16 for x in out_r):
            s = per.get(r16, [])
            out_r.append([smp, r16, FQ[r16][0], slot, cfg, FQ[r16][2], "OK" if s else "NOUNIT",
                          str(len(s)), "%.3f" % (sum(s) / len(s)) if s else "-", str(len(P)), NOTE.get(r16, "-")])
    for fn in glob.glob(os.path.join(wd, "*")):
        os.remove(fn)
    os.rmdir(wd)
    return smp, cfg, out_u, out_r, out_s


def init():
    global CFG, SAVED, ARR, U6, U8
    AGS = samples(os.path.join(D, "AOJV_GIAB.tsv"))
    ACS = samples(os.path.join(D, "AOJV8.tsv"))
    AFS = samples(os.path.join(D, "AOJV5.tsv"))
    SAVED = collections.defaultdict(dict)
    for cfg, fn in (("AOJV3", "AOJV3.tsv"), ("AOJV5", "AOJV5.tsv"), ("AOJV8", "AOJV8.tsv"), ("AOJVG", "AOJV_GIAB.tsv"),
                    ("AOJV_FS", "AOJV_FS.tsv"), ("AOJV_PB", "AOJV_PB.tsv")):   # v12 : Fiber-seq
        if not os.path.exists(os.path.join(D, fn)):
            continue                   # 없는 자료원은 건너뛴다 (v11 동작 보존)
        for L in io.open(os.path.join(D, fn), encoding="utf8"):
            r = L.rstrip(NL).split(TAB)
            if len(r) < 5:
                continue
            e = SAVED[(cfg, r[0])].setdefault(r[1], [r[2], {}])
            e[1][int(r[3])] = set() if r[4] == "-" else set(r[4].split(","))
    # v12 : all available data sources are kept side by side rather than one being chosen
    #   기존 4 종은 예전대로 우선순위로 하나만 고른다 — 같은 read 가 중복되지 않게 (v11 과 동일)
    CFG = {}
    for s in set(k[1] for k in SAVED):
        base = "AOJVG" if s in AGS else "AOJV8" if s in ACS else "AOJV5" if s in AFS else "AOJV3"
        CFG[s] = [base] + [c for c in ("AOJV_FS", "AOJV_PB") if (c, s) in SAVED]
    SEQF = {}
    for L in io.open(os.path.join(D, "g8in", "hap_seq_node_v7.tsv"), encoding="utf8"):
        f = L.rstrip(NL).split(TAB)
        if len(f) >= 3:
            SEQF[(f[0], f[1])] = f[2]
    FLIP = {}
    f = io.open(os.path.join(D, "SUM_G8.tsv"), encoding="utf8")
    f.readline()
    for L in f:
        r = L.rstrip(NL).split(TAB)
        FLIP[(r[0], r[1])] = (r[4] == "Y")
    ARR = {}
    for k, fn in SEQF.items():
        a = rd(os.path.join(D, "cmp4d_v7", fn)).upper()
        ARR[k] = rc(a) if FLIP.get(k) else a
    U6 = loadmap(os.path.join(D, "MAP_G6C2.tsv"), ("UNIT", "UNIT_LAST"))   # 배열 마지막 온전 유닛(UNIT_LAST)도 창·G6 번호에 쓴다 (빠지면 저장 최대번호가 마지막 유닛인 read 가 NO_WINDOW)
    U8 = loadmap(os.path.join(D, "MAP_G8.tsv"))


if __name__ == "__main__":
    init()
    smps = [x.strip() for x in io.open(LST, encoding="utf8") if x.strip()]
    smps = [s for s in smps if s in CFG]
    OU = io.open(os.path.join(D, "AO_G8_UNITS.tsv"), "a", encoding="utf8", newline=NL)
    OR = io.open(os.path.join(D, "AO_G8_READS.tsv"), "a", encoding="utf8", newline=NL)
    os.makedirs(os.path.join(D, "sam"), exist_ok=True)
    tasks = [(s, c) for s in smps for c in CFG[s]]        # v12
    with multiprocessing.Pool(NP, initializer=init) as pool:
        for smp, cfg, uu, rr, ss in pool.imap_unordered(work, tasks):
            for x in uu:
                OU.write(TAB.join(str(y) for y in x) + NL)
            for x in rr:
                OR.write(TAB.join(str(y) for y in x) + NL)
            io.open(os.path.join(D, "sam", "%s.%s.sam" % (smp, cfg)), "w", encoding="utf8", newline=NL).write(NL.join(ss) + (NL if ss else ""))
            OU.flush()
            OR.flush()
            sys.stderr.write("  %s %s 유닛 %d · read %d%s" % (smp, cfg, len(uu), len(rr), NL))
    OU.close()
    OR.close()
    print("done")
