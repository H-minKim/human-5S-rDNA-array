# -*- coding: utf-8 -*-
# 수리 뒤 read 사건 재측정 (18→0 · 10→4)
#   단계 : ③b 정본 수리 검정 (수정1)

# -*- coding: utf-8 -*-
# read 유닛 삽입/결실 전수 v2 (2026-09-16) — junction 을 flank 로 직접 찾아 앵커로 쓴다
#   v1 과 다른 점
#   (a) TUC array 구간에 기대지 않는다 (TUC 는 40유닛 탠덤 참조 정렬이라 긴 read 에서 구간이 안쪽에서 멈춤)
#   (b) read 첫 유전자 앞 꼬리 · 마지막 유전자 뒤 꼬리(각 ≤20kb)를 그 개체 flank(flank_final2, junction 쪽 20kb)에 붙여
#       원위 꼬리 → 원위 flank(L_), 근위 꼬리 → 근위 flank(R_) 이고 flank 의 junction 끝까지(≤JTOL) 닿으면 junction 확정
#   (c) junction 확정 쪽은 유전자 적중을 junction 안쪽 EDGE 까지로 제한 (퇴화 유사 유전자 −238/+481 제외)
#   (d) junction→첫(마지막) 유전자 끝 거리 − 정본 끝 조각(원위 off_prox · 근위 off_dist) = k×2240 ±KTOL 이면
#       빈 유닛 k 개 + junction 토큰(JD/JP). 토큰 짝 = 점수 JSC · 앵커 ANCH 유닛으로 셈. 토큰은 유닛과 짝 불가
#   (e) 양성 대조 : 사건 없는 read 에서 무작위 위치 + junction 바로 옆에 −1/+1 · −2 → 위치(정본 번호·끝까지 거리) 기록
#   나머지(유전자 분할·빈 유닛·끊김·변이·DP·슬롯 판정·사건 인정·음성 대조)는 v1 과 같다.
#   인자 : 병렬수 [개체 목록(쉼표) — 파일럿]
import io, os, re, sys, gzip, glob, random, socket, subprocess, collections, statistics, multiprocessing
NL = chr(10); TAB = chr(9); H = os.path.expanduser("~"); hn = socket.gethostname()
NP = int(sys.argv[1]) if len(sys.argv) > 1 else 8
ONLY = set(sys.argv[2].split(",")) if len(sys.argv) > 2 and sys.argv[2] else None
OUTD = os.path.join(H, "uindel2_fix1" if ONLY is None else "uindel2_fix1"); os.makedirs(OUTD, exist_ok=True)
REF = os.path.join(H, "s5s_tuc", "ref"); GENE = os.path.join(REF, "5S_gene_query.fa"); CANON = os.path.join(REF, "unit_canon.fa")
FIN = os.path.join(H, "flank_final2") if os.path.isdir(os.path.join(H, "flank_final2")) else os.path.join(H, "flank_in", "flank_final2")
LO, HI, GLEN, U, FRAG, PAD = 2100, 2400, 119, 2240, 5, 100
EDGE, TAIL, JTOL, KTOL, JMIN = 150, 20000, 1000, 250, 500
GO, GE, CLEAR, ANCH, MINSEG, JSC = 0.8, 0.1, 0.05, 3, 6, 3.0
CMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
TOK = re.compile(r"^(\d+)([ACGT]|[DI])(\d*)$"); CS = re.compile(r"([:*+\-=])([A-Za-z0-9]+)")
def rc(s): return s.translate(CMP)[::-1]
def rdfa(p): return "".join(x.strip() for x in io.open(p, encoding="utf8") if not x.startswith(">")).upper()
def keep(t):
    m = TOK.match(t); return bool(m) and (m.group(2) in "ACGT" or int(m.group(3) or 0) > 3)
def jac(a, b):
    u = a | b; return len(a & b) / float(len(u)) if u else 0.0

TGT = collections.defaultdict(dict)
for L in io.open(os.path.join(H, "UINDEL_TARGETS_%s.tsv" % hn), encoding="utf8"):
    s, rid, u, fl, st, sl = L.rstrip(NL).split(TAB)
    if ONLY and s not in ONLY: continue
    TGT[s][rid] = (int(u), fl, st, sl)
FQ = collections.defaultdict(list)
for L in io.open(os.path.join(H, "FASTQ_LOC.tsv"), encoding="utf8"):
    h, s, p, b = L.rstrip(NL).split(TAB)
    if h == hn: FQ[s].append(p)
FN = {}
for L in io.open(os.path.join(H, "hap_seq_node_v8.tsv"), encoding="utf8"):
    r = L.rstrip(NL).split(TAB)
    if len(r) >= 3 and r[0] != "sample": FN[(r[0], r[1])] = r[2]
SUMF = {}
f = io.open(os.path.join(H, "g8canonv11f1", "CONSENSUS_SUM_G8_full.tsv"), encoding="utf8"); f.readline()
for L in f:
    r = L.rstrip(NL).split(TAB); SUMF[(r[0], r[1])] = dict(N=int(r[3]), flip=(r[4] == "Y"), op=int(r[5]), od=int(r[6]))

def load_map(smp):
    M = collections.defaultdict(dict)
    f = io.open(os.path.join(H, "g8canonv11f1", "CONSENSUS_MAP_G8_full.tsv"), encoding="utf8"); ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TAB)))
    for L in f:
        if not L.startswith(smp + TAB): continue
        r = L.rstrip(NL).split(TAB)
        if r[ix["seg_type"]] != "UNIT": continue
        ok = r[ix["verify"]] in ("OK", "OK_FILLED")
        M[r[1]][int(r[ix["array_pos"]])] = set(t for t in r[ix["varlist"]].split(",") if keep(t)) if (ok and r[ix["varlist"]] != "-") else set()
    allu = [v for d in M.values() for v in d.values()]; c = collections.Counter(t for v in allu for t in v)
    FIX = set(t for t, x in c.items() if x >= 0.9 * max(1, len(allu)))
    out = {}
    for sl, d in M.items():
        aps = sorted(d); n_ = SUMF.get((smp, sl), {}).get("N", aps[-1] if aps else 0)
        out[sl] = (["JD"] + [d[a] - FIX for a in aps] + ["JP"], [0] + aps + [n_ + 1])
    return out, FIX

def flank_ref(smp, W):
    """개체 슬롯별 원위(L_)·근위(R_) flank 의 junction 쪽 20kb (지도 방향). 없으면 None"""
    p = os.path.join(W, "flank.fa"); n = 0
    with io.open(p, "w", newline=NL) as g:
        for (s, sl), fn in FN.items():
            if s != smp or (s, sl) not in SUMF: continue
            q = os.path.join(FIN, fn)
            if not os.path.exists(q): continue
            h = io.open(q, encoding="utf8").readline().rstrip(NL)
            m = dict(x.split("=", 1) for x in h[1:].split("|") if "=" in x); l5, l3 = int(m.get("flank5", 0)), int(m.get("flank3", 0))
            full = rdfa(q)
            if SUMF[(s, sl)]["flip"]: full = rc(full); left, right = l3, l5
            else: left, right = l5, l3
            if left >= JMIN: g.write(">L_%s%s%s%s" % (sl, NL, full[:left][-TAIL:], NL)); n += 1
            if right >= JMIN: g.write(">R_%s%s%s%s" % (sl, NL, full[len(full) - right:][:TAIL], NL)); n += 1
    return p if n else None

def gene_hits(seqs, W, tag):
    fa = os.path.join(W, "r_%s.fa" % tag); db = os.path.join(W, "db_%s" % tag)
    with io.open(fa, "w", newline=NL) as g:
        for k, s in seqs.items(): g.write(">%s%s%s%s" % (k, NL, s, NL))
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db], capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db, "-outfmt", "6 sseqid sstart send length pident qstart qend",
                        "-dust", "no", "-word_size", "9", "-evalue", "1e-5", "-max_target_seqs", "100000000", "-num_threads", "2"],
                       capture_output=True, text=True)
    HT = collections.defaultdict(list)
    for L in r.stdout.splitlines():
        g = L.split(TAB); a, b, qs_, qe_ = int(g[1]), int(g[2]), int(g[5]), int(g[6]); ext = (qs_ - 1) + (GLEN - qe_)
        if a <= b: HT[g[0]].append(("+", a - (qs_ - 1), b + (GLEN - qe_), int(g[3]), ext))
        else: HT[g[0]].append(("-", b - (GLEN - qe_), a + (qs_ - 1), int(g[3]), ext))
    for fn in glob.glob(os.path.join(W, "*_%s*" % tag)):
        try: os.remove(fn)
        except OSError: pass
    return HT

def gene_ends(hits, lo, hi):
    hs = sorted([h for h in hits if h[0] == "+" and lo <= h[2] <= hi], key=lambda x: x[2])
    keep_ = []
    for h in hs:
        if keep_ and h[1] <= keep_[-1][2]:
            if h[3] > keep_[-1][3]: keep_[-1] = h
        else: keep_.append(h)
    i = 0
    while i + 1 < len(keep_):
        if keep_[i + 1][2] - keep_[i][2] < LO:
            fa_, fb_ = keep_[i][4] > FRAG, keep_[i + 1][4] > FRAG
            if fa_ and not fb_: keep_.pop(i); i = max(0, i - 1); continue
            if fb_ and not fa_: keep_.pop(i + 1); i = max(0, i - 1); continue
            if keep_[i + 1][3] >= keep_[i][3]: keep_.pop(i)
            else: keep_.pop(i + 1)
            i = max(0, i - 1); continue
        i += 1
    return [h[2] for h in keep_]

def chain_from(ends):
    out = []
    for i in range(len(ends) - 1):
        g = ends[i + 1] - ends[i]
        if LO <= g <= HI: out.append(("U", ends[i] + 1, ends[i + 1]))
        else:
            k = int(round(g / float(U)))
            if k >= 2 and abs(g - k * U) <= 120 * k: out.extend([("P",)] * k)
            else: out.append(("B", g))
    return out

def call_vars(useqs, W, tag):
    q = os.path.join(W, "u_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(useqs): g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", "asm20", "-c", "--cs=long", "-t", "2", "--secondary=no", CANON, q], capture_output=True, text=True)
    V = {}
    for L in r.stdout.splitlines():
        f = L.split(TAB); i = int(f[0][1:])
        if i in V: continue
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), ""); tp = int(f[7]); toks = set()
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == "=": tp += len(v)
            elif op == ":": tp += int(v)
            elif op == "*": toks.add("%d%s" % (tp + 1, v[1].upper())); tp += 1
            elif op == "-": toks.add("%dD%d" % (tp + 1, len(v))); tp += len(v)
            elif op == "+": toks.add("%dI%d" % (tp + 1, len(v)))
        V[i] = set(t for t in toks if keep(t))
    os.remove(q)
    return V

def score(a, b, base):
    ta, tb = isinstance(a, str), isinstance(b, str)
    if ta or tb: return JSC if (ta and tb and a == b) else -50.0
    return (jac(a, b) - base) if (a and b) else 0.0

def dp(R, Mv, base):
    n, m = len(R), len(Mv); NEG = -1e18
    S0 = [[NEG] * (m + 1) for _ in range(n + 1)]; S1 = [[NEG] * (m + 1) for _ in range(n + 1)]; S2 = [[NEG] * (m + 1) for _ in range(n + 1)]
    B0 = [[0] * (m + 1) for _ in range(n + 1)]; B1 = [[0] * (m + 1) for _ in range(n + 1)]; B2 = [[0] * (m + 1) for _ in range(n + 1)]
    for j in range(m + 1): S0[0][j] = 0.0
    for i in range(1, n + 1): S1[i][0] = -(GO + GE * (i - 1)) - (50.0 if isinstance(R[i - 1], str) else 0.0); B1[i][0] = 1 if i > 1 else 0
    for i in range(1, n + 1):
        a = R[i - 1]; r0, r1, r2 = S0[i], S1[i], S2[i]; p0, p1, p2 = S0[i - 1], S1[i - 1], S2[i - 1]
        tokA = 50.0 if isinstance(a, str) else 0.0          # read 토큰은 삽입으로 버릴 수 없음
        for j in range(1, m + 1):
            b = Mv[j - 1]; sc = score(a, b, base)
            tokB = 50.0 if (isinstance(b, str) and 1 < j < m) else 0.0
            x0, x1, x2 = p0[j - 1], p1[j - 1], p2[j - 1]
            if x0 >= x1 and x0 >= x2: r0[j] = x0 + sc; B0[i][j] = 0
            elif x1 >= x2: r0[j] = x1 + sc; B0[i][j] = 1
            else: r0[j] = x2 + sc; B0[i][j] = 2
            y0, y1, y2 = p0[j] - GO - tokA, p1[j] - GE - tokA, p2[j] - GO - tokA
            if y0 >= y1 and y0 >= y2: r1[j] = y0; B1[i][j] = 0
            elif y1 >= y2: r1[j] = y1; B1[i][j] = 1
            else: r1[j] = y2; B1[i][j] = 2
            tb = 50.0 if isinstance(b, str) else 0.0             # 정본 토큰은 결실로 건너뛸 수 없음 (양끝 무료는 시작·끝 열에서)
            z0, z1, z2 = r0[j - 1] - GO - tb, r1[j - 1] - GO - tb, r2[j - 1] - GE - tb
            if z0 >= z1 and z0 >= z2: r2[j] = z0; B2[i][j] = 0
            elif z1 >= z2: r2[j] = z1; B2[i][j] = 1
            else: r2[j] = z2; B2[i][j] = 2
    best, je, st = NEG, 0, 0
    for j in range(m + 1):
        for s_, v in ((0, S0[n][j]), (1, S1[n][j]), (2, S2[n][j])):
            if v > best: best, je, st = v, j, s_
    ops = []; i, j = n, je
    while i > 0:
        if st == 0: ops.append(("M", i - 1, j - 1)); nst = B0[i][j]; i -= 1; j -= 1
        elif st == 1: ops.append(("I", i - 1, j)); nst = B1[i][j]; i -= 1
        else: ops.append(("D", i, j - 1)); nst = B2[i][j]; j -= 1
        st = nst
    return best, ops[::-1]

def events(ops, aps, R):
    """사건 = 연속 I 또는 D. 앵커 = 앞뒤 짝 수 (junction 토큰 짝은 ANCH 로 셈)"""
    w = [(ANCH if isinstance(R[o[1]], str) else 1) if o[0] == "M" else 0 for o in ops]
    tot = sum(w); ev = []; seen = 0; k = 0
    while k < len(ops):
        op = ops[k][0]
        if op == "M": seen += w[k]; k += 1; continue
        k2 = k
        while k2 < len(ops) and ops[k2][0] == op: k2 += 1
        if op == "I":
            prevj = next((ops[x][2] for x in range(k - 1, -1, -1) if ops[x][0] == "M"), None)
            pos = aps[prevj] if prevj is not None else None
        else:
            pos = aps[ops[k][2]]
        ev.append((op, k2 - k, pos, seen, tot - seen))
        k = k2
    return ev

def work(smp):
    random.seed(hash(smp) & 0xffffffff)
    W = os.path.join(OUTD, "tmp_" + smp); os.makedirs(W, exist_ok=True)
    want = TGT[smp]; MAPS, FIX = load_map(smp)
    rows, evrows, ctrl = [], [], []
    if len(MAPS) < 1: return smp, rows, evrows, ctrl, "지도 없음"
    SEQ = {}
    for fq in FQ.get(smp, []):
        with gzip.open(fq, "rt") as fh:
            j = 0; rid = None
            for L in fh:
                m = j & 3; j += 1
                if m == 0: rid = L[1:].split()[0][:16]
                elif m == 1 and rid in want and rid not in SEQ: SEQ[rid] = L.strip().upper()
    if not SEQ: return smp, rows, evrows, ctrl, "read 못 찾음"
    HT = gene_hits(SEQ, W, "a")
    flip = {k for k, hs in HT.items() if sum(1 for h in hs if h[0] == "-") > sum(1 for h in hs if h[0] == "+")}
    if flip:
        HT2 = gene_hits(dict((k, rc(SEQ[k])) for k in flip), W, "b")
        for k in flip: SEQ[k] = rc(SEQ[k]); HT[k] = HT2.get(k, [])
    # ── junction 찾기 : 꼬리를 flank 에
    E0 = {k: gene_ends(HT.get(k, []), -10 ** 9, 10 ** 9) for k in SEQ}
    FR = flank_ref(smp, W); JL, JR = {}, {}
    if FR:
        q = os.path.join(W, "tails.fa")
        with io.open(q, "w", newline=NL) as g:
            for k, e in E0.items():
                if not e: continue
                a0 = max(0, e[0] - GLEN - TAIL); a1 = e[0] - GLEN
                if a1 - a0 >= JMIN: g.write(">%s|L|%d%s%s%s" % (k, a0, NL, SEQ[k][a0:a1], NL))
                b0 = e[-1]; b1 = min(len(SEQ[k]), e[-1] + TAIL)
                if b1 - b0 >= JMIN: g.write(">%s|R|%d%s%s%s" % (k, b0, NL, SEQ[k][b0:b1], NL))
        r = subprocess.run(["minimap2", "-x", "map-ont", "-t", "2", "--secondary=no", FR, q], capture_output=True, text=True)
        best = {}
        for L in r.stdout.splitlines():
            g = L.split(TAB)
            qn, ql, qs, qe, strand, tn, tl, ts, te, nmatch = g[0], int(g[1]), int(g[2]), int(g[3]), g[4], g[5], int(g[6]), int(g[7]), int(g[8]), int(g[9])
            if qn not in best or nmatch > best[qn][-1]: best[qn] = (qs, qe, strand, tn, tl, ts, te, nmatch)
        for qn, (qs, qe, strand, tn, tl, ts, te, nm) in best.items():
            k, side, off0 = qn.split("|"); off0 = int(off0)
            if strand != "+" or qe - qs < JMIN: continue
            if side == "L" and tn.startswith("L_") and tl - te <= JTOL:
                JL[k] = (off0 + qe + (tl - te), tn[2:])
            if side == "R" and tn.startswith("R_") and ts <= JTOL:
                JR[k] = (off0 + qs - ts, tn[2:])
    CH, TOKL, TOKR = {}, {}, {}
    for k in SEQ:
        lo = JL[k][0] - EDGE if k in JL else -10 ** 9
        hi = JR[k][0] + EDGE if k in JR else 10 ** 9
        ends = gene_ends(HT.get(k, []), lo, hi); ch = chain_from(ends)
        pre, post = [], []
        if ends and k in JL:
            sl = JL[k][1]; off = SUMF.get((smp, sl), {}).get("op", 0); d = ends[0] - JL[k][0]; kk = int(round((d - off) / float(U)))
            if kk >= 0 and abs(d - off - kk * U) <= KTOL: pre = [("J", "JD")] + [("P",)] * kk; TOKL[k] = kk
            else: TOKL[k] = "resid%+d" % (d - off - kk * U)
        if ends and k in JR:
            sl = JR[k][1]; off = SUMF.get((smp, sl), {}).get("od", 0); d = JR[k][0] - ends[-1]; kk = int(round((d - off) / float(U)))
            if kk >= 0 and abs(d - off - kk * U) <= KTOL: post = [("P",)] * kk + [("J", "JP")]; TOKR[k] = kk
            else: TOKR[k] = "resid%+d" % (d - off - kk * U)
        CH[k] = pre + ch + post
    useqs, uidx = [], []
    for k, ch in CH.items():
        for i, it in enumerate(ch):
            if it[0] == "U": useqs.append(SEQ[k][max(0, it[1] - 1 - PAD):it[2] + PAD]); uidx.append((k, i))
    VV = call_vars(useqs, W, "u") if useqs else {}
    FP = collections.defaultdict(dict)
    for n_, (k, i) in enumerate(uidx): FP[k][i] = VV.get(n_, set()) - FIX
    pool_r = [v for d in FP.values() for v in d.values() if v]; pool_m = [v for (mv, _a) in MAPS.values() for v in mv if isinstance(v, set) and v]
    base = statistics.median([jac(random.choice(pool_r), random.choice(pool_m)) for _ in range(2000)]) if pool_r and pool_m else 0.1
    keep_n, add_n = [], []; tokpool = list(set(t for v in pool_m for t in v)) or ["0A"]; ctrl_todo = []
    for k in SEQ:
        u_tuc, fl, st, slk = want[k]; ch = CH[k]
        segs = []; cur = []
        for i, it in enumerate(ch):
            if it[0] == "B":
                if cur: segs.append(cur)
                cur = []
            elif it[0] == "J": cur.append(it[1])
            else: cur.append(FP[k].get(i, set()) if it[0] == "U" else None)
        if cur: segs.append(cur)
        nU = sum(1 for it in ch if it[0] == "U"); nP = sum(1 for it in ch if it[0] == "P"); nB = sum(1 for it in ch if it[0] == "B")
        segs = [sg for sg in segs if sum(1 for x in sg if not isinstance(x, str)) >= MINSEG]
        jinfo = "L:%s/R:%s" % (TOKL.get(k, "-") if k in JL else "-", TOKR.get(k, "-") if k in JR else "-")
        if not segs:
            rows.append([smp, k, str(u_tuc), fl, st, slk, str(len(SEQ[k])), str(nU), str(nP), str(nB), "0", "-", "none", "-", "-", "-", "-", "0", "-", "-", jinfo])
            continue
        tot = {sl: 0.0 for sl in MAPS}; EV = {sl: [] for sl in MAPS}; ntot = 0; qlist = {sl: [] for sl in MAPS}; SP = {sl: [] for sl in MAPS}
        for si, sg in enumerate(segs):
            R = [x if (x is not None) else set() for x in sg]; ntot += sum(1 for x in R if not isinstance(x, str))
            for sl, (mv, aps) in MAPS.items():
                sc, ops = dp(R, mv, base); tot[sl] += sc; ev = events(ops, aps, R)
                for e in ev: EV[sl].append((si,) + e)
                mj = [aps[o[2]] for o in ops if o[0] == "M" and not isinstance(R[o[1]], str)]
                if mj: SP[sl].append("%d-%d" % (min(mj), max(mj)))
                for o in ops:
                    if o[0] == "M" and isinstance(R[o[1]], set) and R[o[1]] and isinstance(mv[o[2]], set) and mv[o[2]]: qlist[sl].append(jac(R[o[1]], mv[o[2]]))
                if sl == slk and len(R) >= 12 and not ev:
                    for o in ops:
                        if o[0] == "M" and isinstance(R[o[1]], set) and R[o[1]] and isinstance(mv[o[2]], set) and mv[o[2]]:
                            keep_n.append(len(R[o[1]] & mv[o[2]]) / float(len(mv[o[2]]))); add_n.append(len(R[o[1]] - mv[o[2]]))
                    if random.random() < 0.25: ctrl_todo.append((R, sl, ops))
        sls = sorted(MAPS)
        if slk in MAPS: su, basis = slk, "known"
        else:
            order = sorted(sls, key=lambda s: -tot[s]); su = order[0]
            mg = (tot[order[0]] - tot[order[1]]) / max(1, ntot) if len(order) > 1 else 9.0
            basis = "clear" if mg >= CLEAR else "ambig"
        oth = next((s for s in sls if s != su), None)
        mg_out = (tot[su] - tot[oth]) / max(1, ntot) if oth else 9.0
        acc = 0
        for (si, op, size, pos, aL, aR) in EV[su]:
            reason = []
            if aL < ANCH or aR < ANCH: reason.append("앵커부족")
            if basis == "ambig" and not any(e[1] == op and e[2] == size for e in EV.get(oth, [])): reason.append("모호슬롯·반대슬롯에없음")
            expl = "Y" if (oth and tot[oth] >= tot[su] and not EV.get(oth)) else "N"
            ok = "Y" if not reason else "N"; acc += ok == "Y"
            evrows.append([smp, k, st, su, basis, op, str(size), str(pos), str(aL), str(aR), ok, ";".join(reason) or "-", expl, str(si),
                           "%.3f" % (statistics.mean(qlist[su]) if qlist[su] else 0), jinfo])
        rows.append([smp, k, str(u_tuc), fl, st, slk, str(len(SEQ[k])), str(nU), str(nP), str(nB), str(len(segs)), su, basis,
                     "%.4f" % mg_out, "%.2f" % tot[su], "%.2f" % (tot[oth] if oth else 0), "%.3f" % (statistics.mean(qlist[su]) if qlist[su] else 0),
                     str(acc), "|".join("%s%d@%s" % (e[1], e[2], e[3]) for e in EV[su]) or "-", ";".join(SP[su]) or "-", jinfo])
    # ── 양성 대조 : 무작위 위치 + junction 옆 · 위치 기록
    for R, sl, ops in ctrl_todo[:80]:
        mv, aps = MAPS[sl]; n_ = aps[-1] - 1
        idx_ap = {o[1]: aps[o[2]] for o in ops if o[0] == "M"}
        units = [i for i in range(len(R)) if not isinstance(R[i], str)]
        cand = []
        if len(units) >= 8:
            cand.append(("rand", random.choice(units[3:-3])))
        if R and R[0] == "JD" and len(units) >= 5: cand.append(("jxD", units[0]))
        if R and R[-1] == "JP" and len(units) >= 5: cand.append(("jxP", units[-1]))
        for where, i in cand:
            ap = idx_ap.get(i)
            if ap is None: continue
            dend = min(ap, n_ - ap + 1)
            for kind, R2, want_ in (("del1", R[:i] + R[i + 1:], ("D", 1)), ("del2", R[:i] + R[i + 2:] if i + 2 <= len(R) and not isinstance(R[min(i + 1, len(R) - 1)], str) else None, ("D", 2)),
                                    ("ins1", R[:i] + [R[i]] + R[i:], ("I", 1))):
                if R2 is None: continue
                sc, o2 = dp(R2, mv, base); ev = events(o2, aps, R2)
                hit = any(e[0] == want_[0] and e[1] == want_[1] and e[3] >= ANCH and e[4] >= ANCH for e in ev)
                ctrl.append(["pos_%s_%s" % (kind, where), smp, sl, str(len(R2)), "1" if hit else "0", str(sum(1 for e in ev if e[3] >= ANCH and e[4] >= ANCH)), str(ap), str(dend), "D" if ap <= n_ - ap + 1 else "P"])
    # ── 음성 대조 : 지도에서 뗀 가짜 read (+끝이면 토큰) + 실측 잡음
    pk = statistics.mean(keep_n) if keep_n else 0.6; ad = statistics.mean(add_n) if add_n else 2.0
    for t in range(60):
        sl = random.choice(sorted(MAPS)); mv, aps = MAPS[sl]; core = mv[1:-1]
        if len(core) < 12: continue
        L_ = random.randint(12, min(40, len(core))); st0 = random.choice([0, len(core) - L_] + [random.randint(0, len(core) - L_)] * 2)
        R2 = []
        for v in core[st0:st0 + L_]:
            x = set(tk for tk in v if random.random() < pk)
            for _ in range(sum(1 for _ in range(int(ad * 3)) if random.random() < 1 / 3.0)):
                x.add(random.choice(tokpool) if random.random() < 0.5 else "%d%s" % (random.randint(1, 2240), random.choice("ACGT")))
            R2.append(x)
        if st0 == 0: R2 = ["JD"] + R2
        if st0 + L_ == len(core): R2 = R2 + ["JP"]
        sc, o2 = dp(R2, mv, base); ev = events(o2, aps, R2)
        ctrl.append(["neg_sim%s" % ("_jx" if isinstance(R2[0], str) or isinstance(R2[-1], str) else ""), smp, sl, str(L_), "1" if any(e[3] >= ANCH and e[4] >= ANCH for e in ev) else "0",
                     str(sum(1 for e in ev if e[3] >= ANCH and e[4] >= ANCH)), str(aps[1 + st0]), "-", "-"])
    ctrl.append(["noise_model", smp, "-", "-", "%.3f" % pk, "%.3f" % ad, "-", "-", "-"])
    ctrl.append(["junction", smp, "-", "-", "L%d/R%d" % (len(JL), len(JR)), "tokL%d/tokR%d" % (sum(1 for v in TOKL.values() if isinstance(v, int)), sum(1 for v in TOKR.values() if isinstance(v, int))), "-", "-", "-"])
    for fn in glob.glob(os.path.join(W, "*")):
        try: os.remove(fn)
        except OSError: pass
    try: os.rmdir(W)
    except OSError: pass
    return smp, rows, evrows, ctrl, "ok read %d · junction L%d R%d" % (len(SEQ), len(JL), len(JR))

if __name__ == "__main__":
    RO = io.open(os.path.join(OUTD, "UINDEL_READS_%s.tsv" % hn), "w", encoding="utf8", newline=NL)
    EO = io.open(os.path.join(OUTD, "UINDEL_EVENTS_%s.tsv" % hn), "w", encoding="utf8", newline=NL)
    CO = io.open(os.path.join(OUTD, "UINDEL_CTRL_%s.tsv" % hn), "w", encoding="utf8", newline=NL)
    LG = io.open(os.path.join(OUTD, "UINDEL_LOG_%s.txt" % hn), "w", encoding="utf8", newline=NL)
    RO.write(TAB.join(["sample", "read16", "units_tuc", "flank", "status", "slot_known", "read_len", "n_unit", "n_placeholder", "n_break", "n_seg",
                       "slot_used", "slot_basis", "margin_per_unit", "dp_own", "dp_other", "mean_jac_own", "n_events_accepted", "events_all", "aligned_ap_spans", "junction_k"]) + NL)
    EO.write(TAB.join(["sample", "read16", "status", "slot_used", "slot_basis", "type", "size_units", "canon_ap", "anchor_L", "anchor_R",
                       "accepted", "reject_reason", "other_slot_explains", "seg", "mean_jac_own", "junction_k"]) + NL)
    CO.write(TAB.join(["kind", "sample", "slot", "n_units", "hit", "n_events_anchored", "canon_ap", "dist_to_end", "side"]) + NL)
    with multiprocessing.Pool(NP) as pool:
        for smp, rows, evrows, ctrl, msg in pool.imap_unordered(work, sorted(TGT)):
            for r in rows: RO.write(TAB.join(r) + NL)
            for r in evrows: EO.write(TAB.join(r) + NL)
            for r in ctrl: CO.write(TAB.join(r) + NL)
            LG.write("%s\t%s\n" % (smp, msg)); RO.flush(); EO.flush(); CO.flush(); LG.flush()
    LG.write("DONE\n"); LG.close()
