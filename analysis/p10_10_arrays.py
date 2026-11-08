# -*- coding: utf-8 -*-
# P10 계열별 : [B] 이봉성  [C] array 절대축 지도  [D] 전이의 sub-unit 위치
# 입력 : out/P10/P10_FRAME_ALL.tsv (p10_08, 실측 주기판) · P10_FRAME_SEGS_ALL.tsv · P10_SUBUNIT_FRAMES.json
#
# Rules
#  · 유닛 class = 계열별 2-성분 이항혼합 EM 의 우도비 (0.5 고정 임계 아님)
#  · 결합도 = sub-unit 을 뺀 나머지로 class 를 정하고 그 sub-unit 의 β 차 / 동적범위 (순환 회피)
#  · 경계 프로파일 = 같은 canonical 칸의 LL/HH 대조로 정규화한 score, 지속 교차(200 bp)
#  · DXZ4 는 성별로 나눈다 (Xa/Xi). 극성 CONFLICT 인 DXZ4 조각은 원위/근위를 부르지 않는다.
import io, os, math, json, collections, statistics

def _need(_v):
    """No defaults are set: the script stops if a variable is missing.
    이 글이 읽는 자료의 뿌리는 README.md 의
    "Running this code" 절에 적혀 있다."""
    import os as _o, sys as _s
    _x = _o.environ.get(_v)
    if not _x:
        _s.exit(
            "environment variable %s is not set; it is required by %s"
            % (_v, _o.path.basename(__file__))
            + "\n(a module stops the script that imports it)."
            + "\nThe deposited copy has no default; see README.md"
            + " section \"Running this code\"."
        )
    return _x


W = _need("RDNA_WORK_ROOT")
OUT = os.path.join(W, "out", "P10")
FR = json.load(io.open(os.path.join(OUT, "P10_SUBUNIT_FRAMES.json"), encoding="utf8"))
SEX = {}
for i, L in enumerate(io.open(os.path.join(W, "P10", "out", "P10_SAMPLES.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    SEX[f[0]] = f[2]
BIN, MINC = 40, 20
# 5S 는 canonical 이 CHM13 c 와 반대 방향
CDIR = {"5S": -1, "D4Z4_4q": 1, "D4Z4_10q": 1, "DXZ4": 1, "RNU2": 1, "SST1": 1}


def grp(fam, smp):
    if fam == "DXZ4":
        return "DXZ4_F" if SEX.get(smp) == "XX" else "DXZ4_M"
    return fam


def base(g):
    return "DXZ4" if g.startswith("DXZ4") else g


GROUPS = ["5S", "D4Z4_4q", "D4Z4_10q", "DXZ4_F", "DXZ4_M", "RNU2", "SST1"]

# ---------- 적재 ----------
# 완전 유닛만 : read 가 유닛 중간에서 끝나면 그 read-유닛은 일부 구획만 덮는다.
#   RNU2 에서 부분 유닛 LOW 비율 0.314 vs 완전 유닛 0.004 — 비메틸 구획(U2·CpG섬)만 덮은 조각이
#   가짜 저메틸 유닛을 만들었다. CpG 가 있는 칸의 90% 이상을 덮은 read-유닛만 쓴다.
RAW = []
TOT = collections.defaultdict(collections.Counter)
for i, L in enumerate(io.open(os.path.join(OUT, "P10_FRAME_ALL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    cb = {}
    for t in f[13].split(";"):
        b, n, k = t.split(":")
        cb[int(b)] = (int(n), int(k))
        TOT[f[2]][int(b)] += int(n)
    RAW.append((f, cb))
CPGB = {}
for fam, c in TOT.items():
    mx = max(c.values())
    CPGB[fam] = {b for b, v in c.items() if v >= 0.02 * mx}
COVMIN = 0.9
RU = collections.defaultdict(list)    # group -> [(smp,hap,md,strand,pol,nun,read,idx,ufd,n,k,cbins{b:(n,k)})]
NALL = collections.Counter()
for f, cb in RAW:
    fam = f[2]
    NALL[grp(fam, f[0])] += 1
    if len(set(cb) & CPGB[fam]) / float(len(CPGB[fam])) < COVMIN:
        continue
    RU[grp(fam, f[0])].append((f[0], f[1], f[3], f[4], f[6], int(f[7]), f[8], int(f[9]),
                               None if f[10] == "NA" else int(f[10]), int(f[11]), int(f[12]), cb))
del RAW
print("적재(완전 유닛 / 전체) : " + " · ".join("%s %d/%d" % (g, len(RU[g]), NALL[g]) for g in GROUPS))
# 하플로타입당 조각 수 (조립 끊김 표시)
_segs = set()
for i, L in enumerate(io.open(os.path.join(OUT, "P10_FRAME_SEGS_ALL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    _segs.add((f[0], f[1], f[2], f[3]))
NSEG = collections.Counter((s, h, fam) for s, h, fam, _m in _segs)
print("하플로당 조각 수 분포 : " + " · ".join(
    "%s %s" % (fam, dict(collections.Counter(v for (s, h, f_), v in NSEG.items() if f_ == fam)))
    for fam in ("5S", "D4Z4_4q", "D4Z4_10q", "DXZ4", "RNU2", "SST1")))

O = io.open(os.path.join(OUT, "P10_ARRAYS.tsv"), "w", encoding="utf8", newline="\n")
O.write("group\tkey\tvalue\n")


def w(g, k, v):
    O.write("%s\t%s\t%s\n" % (g, k, v))


# ---------- 이봉성 도구 ----------
def em_binom(dat, p0=0.1, p1=0.85, w0=0.3):
    for _ in range(300):
        s0 = n0 = k0 = n1 = k1 = 0.0
        for n, k in dat:
            a = math.log(max(w0, 1e-9)) + k * math.log(max(p0, 1e-9)) + (n - k) * math.log(max(1 - p0, 1e-9))
            b = math.log(max(1 - w0, 1e-9)) + k * math.log(max(p1, 1e-9)) + (n - k) * math.log(max(1 - p1, 1e-9))
            m = max(a, b)
            r = math.exp(a - m) / (math.exp(a - m) + math.exp(b - m))
            s0 += r
            n0 += r * n
            k0 += r * k
            n1 += (1 - r) * n
            k1 += (1 - r) * k
        w0n, p0n, p1n = s0 / len(dat), (k0 / n0 if n0 else p0), (k1 / n1 if n1 else p1)
        if max(abs(w0n - w0), abs(p0n - p0), abs(p1n - p1)) < 1e-7:
            w0, p0, p1 = w0n, p0n, p1n
            break
        w0, p0, p1 = w0n, p0n, p1n
    return w0, p0, p1


def gmm2(x):
    xs = sorted(x)
    m0, m1 = xs[len(xs) // 5], xs[4 * len(xs) // 5]
    s0 = s1 = 0.1
    w0 = 0.5
    for _ in range(300):
        a0 = a1 = b0 = b1 = c0 = c1 = 0.0
        for v in x:
            p = w0 * math.exp(-(v - m0) ** 2 / (2 * s0 * s0)) / s0
            q = (1 - w0) * math.exp(-(v - m1) ** 2 / (2 * s1 * s1)) / s1
            r = p / (p + q) if p + q > 0 else 0.5
            a0 += r
            b0 += r * v
            c0 += r * v * v
            a1 += 1 - r
            b1 += (1 - r) * v
            c1 += (1 - r) * v * v
        m0n, m1n = b0 / a0, b1 / a1
        s0 = max(0.02, math.sqrt(max(c0 / a0 - m0n ** 2, 1e-6)))
        s1 = max(0.02, math.sqrt(max(c1 / a1 - m1n ** 2, 1e-6)))
        w0n = a0 / len(x)
        if abs(m0n - m0) < 1e-7 and abs(m1n - m1) < 1e-7:
            m0, m1, w0 = m0n, m1n, w0n
            break
        m0, m1, w0 = m0n, m1n, w0n
    if m0 > m1:
        m0, m1, s0, s1, w0 = m1, m0, s1, s0, 1 - w0
    return w0, m0, s0, m1, s1


def bc(x):
    n = len(x)
    m = statistics.mean(x)
    s = statistics.pstdev(x)
    g = sum((v - m) ** 3 for v in x) / n / s ** 3
    k = sum((v - m) ** 4 for v in x) / n / s ** 4 - 3
    return (g * g + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))


FIT = {}
print("\n" + "=" * 100)
print("[B] 유닛 β 이봉성 (배치 read-유닛, CpG>=%d)" % MINC)
print("  %-9s %6s %8s %8s %6s %6s %6s %6s %6s %6s  %s"
      % ("그룹", "유닛", "p_LOW", "p_HIGH", "w_lo", "μ_lo", "μ_hi", "D", "BC", "LOW≤.3", "판정"))
for g in GROUPS:
    rows = [r for r in RU[g] if r[9] >= MINC]
    if len(rows) < 200:
        continue
    dat = [(r[9], r[10]) for r in rows]
    x = [k / float(n) for n, k in dat]
    w0, p0, p1 = em_binom(dat)
    gw, gm0, gs0, gm1, gs1 = gmm2(x)
    D = abs(gm1 - gm0) * math.sqrt(2) / math.sqrt(gs0 ** 2 + gs1 ** 2)
    B = bc(x)
    call = "이봉" if (D > 2 and B > 0.555 and min(gw, 1 - gw) >= 0.05) else "단봉"
    FIT[g] = (p0, p1)
    lo = sum(1 for v in x if v <= 0.3) / float(len(x))
    print("  %-9s %6d %8.3f %8.3f %6.3f %6.3f %6.3f %6.2f %6.3f %6.3f  %s"
          % (g, len(x), p0, p1, gw, gm0, gm1, D, B, lo, call))
    h = [0] * 50
    for v in x:
        h[min(49, int(v * 50))] += 1
    for k, v in (("n_unit", len(x)), ("p_low", p0), ("p_high", p1), ("gmm_w_low", gw),
                 ("gmm_mu_low", gm0), ("gmm_mu_high", gm1), ("ashman_D", D), ("BC", B),
                 ("frac_le03", lo), ("call", call),
                 ("hist50", ",".join("%.5f" % (c / float(len(x))) for c in h))):
        w(g, k, v)
    # 하플로타입별 LOW 비율 (DXZ4 는 Xa/Xi 가 하플로 단위라 여기서 보인다)
    hp = collections.defaultdict(list)
    for r in rows:
        hp[(r[0], r[1])].append(r[10] / float(r[9]))
    hl = sorted((sum(1 for v in xx if v <= 0.3) / float(len(xx)), s, h_, len(xx))
                for (s, h_), xx in hp.items() if len(xx) >= 30)
    w(g, "hap_lowfrac", ";".join("%s_%s:%.3f:%d" % (s, h_, f_, n_) for f_, s, h_, n_ in hl))

# ---------- class 호출 ----------


def cls_of(g, n, k):
    p0, p1 = FIT[g]
    d = k * math.log(p1 / p0) + (n - k) * math.log((1 - p1) / (1 - p0))
    return 1 if d > 0 else 0


# ---------- [C] array 절대축 ----------
print("\n" + "=" * 100)
print("[C] array 절대축 — 원위 끝부터 / 근위 끝부터 (극성 확정 조각만) · 유닛 β 평균")
for g in GROUPS:
    rows = [r for r in RU[g] if r[9] >= MINC]
    if len(rows) < 200:
        continue
    D_, P_, E_ = collections.defaultdict(list), collections.defaultdict(list), collections.defaultdict(list)
    for r in rows:
        b = r[10] / float(r[9])
        nun = r[5]
        if r[8] is not None:
            D_[r[8]].append(b)
            P_[nun - r[8] + 1].append(b)
        # 극성 몰라도 되는 축 : 가까운 끝까지 거리
        e = min(r[7] + 1, nun - r[7])
        E_[e].append(b)
    nd = sum(len(v) for v in D_.values())

    def prof(d, K=40):
        return [(u, statistics.mean(d[u]), len(d[u])) for u in range(1, K + 1) if len(d[u]) >= 20]
    dp, pp, ep = prof(D_), prof(P_), prof(E_)

    def seg(d, a, z):
        v = [x for u in range(a, z + 1) for x in d[u]]
        return (statistics.mean(v), len(v)) if v else (float("nan"), 0)
    inner = [x for u, v in D_.items() for x in v if u > 20 and (u in D_)]
    # 내부 = 양 끝에서 20 유닛 넘게
    inner = []
    for r in rows:
        if r[8] is not None and r[8] > 20 and (r[5] - r[8] + 1) > 20:
            inner.append(r[10] / float(r[9]))
    d15, p15 = seg(D_, 1, 5), seg(P_, 1, 5)
    e15 = seg(E_, 1, 5)
    ein = [x for u, v in E_.items() if u > 20 for x in v]
    print("  %-9s 극성확정 %5d/%5d | 원위1-5 %.3f(%d) 내부(>20) %s 근위1-5 %.3f(%d) | 끝1-5(극성무관) %.3f · 끝>20 %s"
          % (g, nd, len(rows), d15[0], d15[1],
             "%.3f(%d)" % (statistics.mean(inner), len(inner)) if inner else "-",
             p15[0], p15[1], e15[0], "%.3f" % statistics.mean(ein) if ein else "-"))
    w(g, "axis_distal", ";".join("%d:%.4f:%d" % x for x in dp))
    w(g, "axis_proximal", ";".join("%d:%.4f:%d" % x for x in pp))
    w(g, "axis_nearest_end", ";".join("%d:%.4f:%d" % x for x in ep))
    w(g, "distal1_5", "%.4f" % d15[0])
    w(g, "proximal1_5", "%.4f" % p15[0])
    w(g, "interior_gt20", "%.4f" % statistics.mean(inner) if inner else "nan")
    # 민감도 : 조립본이 array 를 끊으면 조각 끝 ≠ array 끝. 하플로당 조각 1 개인 것만
    one = [r for r in rows if NSEG.get((r[0], r[1], base(g)), 0) == 1 and r[8] is not None]
    if one:
        d1 = [r[10] / float(r[9]) for r in one if r[8] <= 5]
        p1_ = [r[10] / float(r[9]) for r in one if r[5] - r[8] + 1 <= 5]
        i1 = [r[10] / float(r[9]) for r in one if r[8] > 20 and r[5] - r[8] + 1 > 20]
        nh = len(set((r[0], r[1]) for r in one))
        f_ = lambda v: "%.3f(%d)" % (statistics.mean(v), len(v)) if v else "-"
        print("  %-9s   └ 하플로당 조각 1개인 %d 하플로만 : 원위1-5 %s · 내부 %s · 근위1-5 %s"
              % ("", nh, f_(d1), f_(i1), f_(p1_)))
        w(g, "onesegment_haps", nh)
        w(g, "onesegment_d_i_p", "%s|%s|%s" % (f_(d1), f_(i1), f_(p1_)))

# ---------- canonical 칸 ----------


def canon_bins(fam, cb):
    """CHM13 c 칸 → canonical 칸 (칸 중심의 canonical 위치로)"""
    fr = FR[fam]
    L = fr["L"]
    out = collections.defaultdict(lambda: [0, 0])
    for b, (n, k) in cb.items():
        c = min(L - 1, b * BIN + BIN // 2)
        q = fr["canon"][c]
        v = out[(q - 1) // BIN]
        v[0] += n
        v[1] += k
    return out


LABQ = {}              # fam -> [라벨] (canonical q 1..maxq, 빈 칸은 가까운 값으로)
for fam, fr in FR.items():
    d = {q: lab for q, lab in zip(fr["canon"], fr["label"])}
    mq = max(fr["canon"])
    arr = [None] * (mq + 2)
    last = None
    for q in range(1, mq + 1):
        if q in d:
            last = d[q]
        arr[q] = last
    nxt = None
    for q in range(mq, 0, -1):
        if arr[q] is None:
            arr[q] = nxt
        else:
            nxt = arr[q]
    LABQ[fam] = arr


def labq(fam, b):
    a = LABQ[fam]
    return a[max(1, min(len(a) - 2, b * BIN + BIN // 2))]

# ---------- [D-1] sub-unit 결합도 (leave-one-subunit-out) ----------
print("\n" + "=" * 100)
print("[D-1] sub-unit 결합도 — 그 sub-unit 을 뺀 나머지로 class 를 정하고, 그 sub-unit 이 따라가나")
print("   결합도 = (HIGH일 때 β − LOW일 때 β) / (그 sub-unit 의 HIGH−LOW 최대폭 : 전 sub-unit 중 최대)")
for g in GROUPS:
    fam = base(g)
    rows = [r for r in RU[g] if r[9] >= MINC]
    if len(rows) < 200 or g not in FIT:
        continue
    order = FR[fam]["order"]
    if len(order) < 2:
        print("  %-9s : 주석 sub-unit 없음 → 칸 단위만 (아래 D-2)" % g)
        continue
    acc = {s: {0: [0, 0], 1: [0, 0]} for s in order}
    for r in rows:
        qb = canon_bins(fam, r[11])
        per = collections.defaultdict(lambda: [0, 0])
        for b, (n, k) in qb.items():
            s = labq(fam, b)
            if s is None:
                continue
            per[s][0] += n
            per[s][1] += k
        tn, tk = r[9], r[10]
        for s in order:
            n, k = per[s]
            if n < 3:
                continue
            no, ko = tn - n, tk - k
            if no < 15:
                continue
            c = cls_of(g, no, ko)
            acc[s][c][0] += n
            acc[s][c][1] += k
    res = []
    for s in order:
        a = acc[s]
        if a[0][0] < 200 or a[1][0] < 200:
            res.append((s, None, None))
            continue
        res.append((s, a[0][1] / float(a[0][0]), a[1][1] / float(a[1][0])))
    ok = [x for x in res if x[1] is not None]
    if not ok:
        continue
    dyn = max(x[2] - x[1] for x in ok)
    print("  %-9s %s" % (g, " | ".join("%s L%.3f H%.3f 결합%.2f" % (s, lo, hi, (hi - lo) / dyn)
                                      for s, lo, hi in ok)))
    w(g, "subunit_coupling", ";".join("%s:%.4f:%.4f:%.4f" % (s, lo, hi, (hi - lo) / dyn) for s, lo, hi in ok))

# ---------- [D-2] canonical 칸 프로파일 (class 별) + 경계 프로파일 ----------
print("\n" + "=" * 100)
print("[D-2] 경계 — 같은 read 에서 canonical 순서로 이웃한 두 유닛의 class 가 다를 때")
print("   새 블록 첫 유닛의 canonical 칸 score (0 = 그 칸의 LOW 내부, 1 = HIGH 내부), 지속 교차 200 bp")
for g in GROUPS:
    fam = base(g)
    rows = [r for r in RU[g] if r[9] >= MINC]
    if len(rows) < 200 or g not in FIT:
        continue
    L = FR[fam]["L"]
    nb = L // BIN + 1
    byread = collections.defaultdict(dict)
    for r in rows:
        byread[(r[0], r[2], r[6])][r[7]] = r
    CT = {0: collections.defaultdict(lambda: [0, 0]), 1: collections.defaultdict(lambda: [0, 0])}
    WIN = {"H2L": {0: collections.defaultdict(lambda: [0, 0]), 1: collections.defaultdict(lambda: [0, 0])},
           "L2H": {0: collections.defaultdict(lambda: [0, 0]), 1: collections.defaultdict(lambda: [0, 0])}}
    NEV = collections.Counter()
    for key, u in byread.items():
        for idx, r in u.items():
            sgn = 1 if r[3] == "+" else -1
            step = sgn * CDIR[fam]              # canonical 다음 유닛
            c = cls_of(g, r[9], r[10])
            prv, nxt = u.get(idx - step), u.get(idx + step)
            if prv and nxt and cls_of(g, prv[9], prv[10]) == c and cls_of(g, nxt[9], nxt[10]) == c:
                for b, (n, k) in canon_bins(fam, r[11]).items():
                    CT[c][b][0] += n
                    CT[c][b][1] += k
            if nxt is None:
                continue
            cn = cls_of(g, nxt[9], nxt[10])
            if cn == c:
                continue
            # 청정 : 양쪽에 한 유닛씩 더 있고 같은 class (없으면 세지 않는다 — v11 은 3 개씩)
            p2, n2 = u.get(idx - step), u.get(idx + 2 * step)
            clean = (p2 is not None and cls_of(g, p2[9], p2[10]) == c) and \
                    (n2 is not None and cls_of(g, n2[9], n2[10]) == cn)
            if not clean:
                continue
            kind = "H2L" if c == 1 else "L2H"
            NEV[kind] += 1
            for off, rr in ((0, r), (1, nxt)):
                for b, (n, k) in canon_bins(fam, rr[11]).items():
                    WIN[kind][off][b][0] += n
                    WIN[kind][off][b][1] += k
    line = []
    for kind in ("H2L", "L2H"):
        if NEV[kind] < 60:
            line.append("%s %d건(부족)" % (kind, NEV[kind]))
            w(g, "%s_n" % kind, NEV[kind])
            continue
        rowsc = {}
        for off in (0, 1):
            rr = []
            for b in range(nb):
                a = WIN[kind][off][b]
                lo, hi = CT[0][b], CT[1][b]
                if a[0] < 40 or lo[0] < 100 or hi[0] < 100:
                    rr.append(None)
                    continue
                bl, bh = lo[1] / float(lo[0]), hi[1] / float(hi[0])
                if bh - bl < 0.1:
                    rr.append(None)
                    continue
                rr.append((a[1] / float(a[0]) - bl) / (bh - bl))
            rowsc[off] = rr
        # 지속 교차 : 새 유닛(off 1) 에서, 옛 상태 쪽 → 새 상태 쪽으로 넘어가 200 bp 유지
        tgt = 0.5
        hf = (kind == "H2L")
        cs = None
        r1 = rowsc[1]
        for b in range(nb):
            if r1[b] is None:
                continue
            cnt, ok = 0, True
            for t in range(b, nb):
                if r1[t] is None:
                    continue
                if (hf and r1[t] <= tgt) or ((not hf) and r1[t] >= tgt):
                    cnt += 1
                else:
                    ok = False
                    break
                if cnt >= 5:
                    break
            if ok and cnt >= 5:
                cs = b
                break
        # 첫 유닛 sub-unit 별 평균 score
        subsc = collections.defaultdict(list)
        for b, v in enumerate(r1):
            if v is None:
                continue
            subsc[labq(fam, b)].append(v)
        sub_txt = " ".join("%s %.2f" % (s, statistics.mean(subsc[s])) for s in FR[fam]["order"] if subsc.get(s))
        if cs is not None:
            q = cs * BIN + BIN // 2
            lab = labq(fam, cs)
            line.append("%s %d건: 교차 q≈%d (%s) | 첫유닛 %s" % (kind, NEV[kind], q, lab, sub_txt))
            w(g, "%s_cross_q" % kind, q)
            w(g, "%s_cross_sub" % kind, lab)
        else:
            line.append("%s %d건: 지속교차 없음 | 첫유닛 %s" % (kind, NEV[kind], sub_txt))
        w(g, "%s_n" % kind, NEV[kind])
        for off in (0, 1):
            w(g, "%s_score_off%d" % (kind, off),
              ",".join("%.4f" % v if v is not None else "" for v in rowsc[off]))
    print("  %-9s %s" % (g, "\n            ".join(line)))
    # class 별 canonical 프로파일 (그림용)
    for c in (0, 1):
        w(g, "ctrl_prof_%d" % c, ",".join("%.4f" % (CT[c][b][1] / float(CT[c][b][0])) if CT[c][b][0] >= 100 else ""
                                         for b in range(nb)))
O.close()
print("\nsaved out/P10/P10_ARRAYS.tsv")
