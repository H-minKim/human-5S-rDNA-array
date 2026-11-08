# -*- coding: utf-8 -*-
# Ask A (2차) — 요소마다 이항 혼합(EM)으로 두 상태를 독립 추정하고, 그 문턱으로 전이 위치를 본다.
#
# 1차 실패: 히스토그램에서 봉우리·골을 찾는 방식은 CpG 가 적은 요소(유전자 9·프로모터 11)에서
#   β 가 이산적이라 β=1.0 막대 옆 빈 칸을 '골'로 잡았다(0.95~0.97). 폐기.
#
# 이번 방법: 요소마다 (n_cpg, n_meth) 로 2성분 이항 혼합을 EM 으로 적합 →
#   p_low, p_high, 혼합비 w. 분류는 균등 사전(우도비)로 해서 위치 곡선이 w 에 끌리지 않게 한다.
#   이분성 판정 = 적합된 혼합이 관측 분포를 설명하는가 + '애매한 유닛'(우도비 1 근처) 비율.
import io, os, gzip, collections, math, statistics

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
V = os.path.join(_need("RDNA_VERSION_ROOT"), "08_versions", "v11")
SU = os.path.join(W, "out", "SUBUNIT.tsv.gz")
SUM = os.path.join(V, "03_maps", "g8canonv10", "CONSENSUS_SUM_G8_full.tsv")
EL = ["AluY", "NTS1", "GA", "NTS2", "prom", "gene"]
MINC = {"AluY": 12, "NTS1": 30, "GA": 12, "NTS2": 20, "prom": 6, "gene": 5}

BRK = set()
for i, L in enumerate(io.open(SUM, encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    if len(f) >= 8 and int(f[7]) > 0:
        BRK.add((f[0], f[1]))

CNT = collections.defaultdict(collections.Counter)       # el -> (n,k) -> 수
ROWS = collections.defaultdict(list)                     # el -> [(axis,idx,n,k)]
for i, L in enumerate(gzip.open(SU, "rt", encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    if len(f) < 20 or (f[0], f[1]) in BRK:
        continue
    pi, di = int(f[3]), int(f[4])
    ax, idx = ("prox", pi) if pi <= di else ("dist", di)
    if idx > 40:
        ax, idx = "int", 0
    for j, e in enumerate(EL):
        nn, bb = f[8 + 2 * j], f[9 + 2 * j]
        if nn and bb and int(nn) >= MINC[e]:
            n_ = int(nn)
            k_ = int(round(float(bb) * n_))
            CNT[e][(n_, k_)] += 1
            ROWS[e].append((ax, idx, n_, k_))
print("요소별 유닛 수:", {e: len(ROWS[e]) for e in EL})


def logbinom_pmf(k, n, p):
    if p <= 0:
        return 0.0 if k == 0 else -1e300
    if p >= 1:
        return 0.0 if k == n else -1e300
    return k * math.log(p) + (n - k) * math.log(1 - p)


def em(cnt, iters=300):
    """2성분 이항 혼합. 조합계수는 두 성분에 공통이라 생략 가능."""
    plo, phi, w = 0.05, 0.85, 0.7
    items = list(cnt.items())
    for _ in range(iters):
        s_lo_k = s_lo_n = s_hi_k = s_hi_n = 0.0
        s_w = s_tot = 0.0
        for (n_, k_), c in items:
            a = math.log(w) + logbinom_pmf(k_, n_, phi)
            b = math.log(1 - w) + logbinom_pmf(k_, n_, plo)
            m = max(a, b)
            ra = math.exp(a - m)
            rb = math.exp(b - m)
            gh = ra / (ra + rb)
            s_hi_k += c * gh * k_
            s_hi_n += c * gh * n_
            s_lo_k += c * (1 - gh) * k_
            s_lo_n += c * (1 - gh) * n_
            s_w += c * gh
            s_tot += c
        np_hi = s_hi_k / s_hi_n if s_hi_n else phi
        np_lo = s_lo_k / s_lo_n if s_lo_n else plo
        nw = s_w / s_tot
        if abs(np_hi - phi) + abs(np_lo - plo) + abs(nw - w) < 1e-9:
            phi, plo, w = np_hi, np_lo, nw
            break
        phi, plo, w = np_hi, np_lo, nw
    return plo, phi, w


FIT = {}
print("\n[A-1] 요소마다 독립으로 적합한 두 상태 (유닛 class 를 전혀 쓰지 않음)")
print("%-6s %9s %8s %9s %9s %8s %10s %10s"
      % ("요소", "유닛", "CpG중앙", "저상태 p", "고상태 p", "고비율 w", "애매유닛", "귀무 애매"))
O = io.open(os.path.join(W, "out", "SUBUNIT_BIMODAL2.tsv"), "w", encoding="utf8", newline="\n")
O.write("kind\telement\tn_unit\tmed_ncpg\tp_low\tp_high\tw_high\tambig_obs\tambig_null\n")
for e in EL:
    plo, phi, w = em(CNT[e])
    FIT[e] = (plo, phi, w)
    ncs = [n_ for (n_, _k), c in CNT[e].items() for _ in range(min(c, 1))]
    med_n = statistics.median([n_ for (n_, _k) in CNT[e] for _ in range(CNT[e][(n_, _k)])][:400000])
    # 애매 유닛 = 우도비가 1/9 ~ 9 사이 (posterior 0.1~0.9, 균등 사전)
    amb = tot = 0
    for (n_, k_), c in CNT[e].items():
        d = logbinom_pmf(k_, n_, phi) - logbinom_pmf(k_, n_, plo)
        tot += c
        if abs(d) < math.log(9):
            amb += c
    # 귀무 : 진짜 두 상태뿐일 때 표집만으로 생기는 애매 비율
    nullamb = nulltot = 0.0
    for (n_, k_), c in CNT[e].items():
        pass
    ns = collections.Counter()
    for (n_, _k), c in CNT[e].items():
        ns[n_] += c
    for n_, c in ns.items():
        for comp, wt in ((phi, w), (plo, 1 - w)):
            # 이 성분에서 k 분포
            lg = [0.0] * (n_ + 1)
            # log C(n,k) 재귀
            lc = 0.0
            for k_ in range(n_ + 1):
                if k_ > 0:
                    lc += math.log((n_ - k_ + 1) / float(k_))
                lg[k_] = lc + logbinom_pmf(k_, n_, comp)
            mx = max(lg)
            ps = [math.exp(x - mx) for x in lg]
            s = sum(ps)
            for k_ in range(n_ + 1):
                d = logbinom_pmf(k_, n_, phi) - logbinom_pmf(k_, n_, plo)
                if abs(d) < math.log(9):
                    nullamb += c * wt * ps[k_] / s
            nulltot += c * wt
    print("%-6s %9d %8d %9.4f %9.4f %8.3f %10.4f %10.4f"
          % (e, tot, med_n, plo, phi, w, amb / tot, nullamb / nulltot))
    O.write("fit\t%s\t%d\t%d\t%.5f\t%.5f\t%.4f\t%.5f\t%.5f\n"
            % (e, tot, med_n, plo, phi, w, amb / tot, nullamb / nulltot))
print("  애매유닛 ≈ 귀무 이면 중간 상태가 없다는 뜻 = 그 요소도 이분")

print("\n[A-2] 요소마다 자기 우도비로 분류한 HIGH 비율 (균등 사전)")
KEYS = [("prox", 1), ("prox", 3), ("prox", 5), ("prox", 10), ("prox", 20), ("int", 0),
        ("dist", 20), ("dist", 16), ("dist", 12), ("dist", 10), ("dist", 5), ("dist", 1)]
print("%-6s %s" % ("요소", " ".join("%5s" % x for x in
                                    ["p1", "p3", "p5", "p10", "p20", "int",
                                     "d20", "d16", "d12", "d10", "d5", "d1"])))
CURVE = {}
for e in EL:
    plo, phi, w = FIT[e]
    acc = collections.defaultdict(lambda: [0, 0])
    for ax, idx, n_, k_ in ROWS[e]:
        d = logbinom_pmf(k_, n_, phi) - logbinom_pmf(k_, n_, plo)
        a = acc[(ax, idx)]
        a[0] += 1
        a[1] += 1 if d > 0 else 0
    C = {k: v[1] / float(v[0]) for k, v in acc.items() if v[0] >= 200}
    CURVE[e] = C
    print("%-6s %s" % (e, " ".join("%5.2f" % C[k] if k in C else "    -" for k in KEYS)))
    for k, v in sorted(C.items()):
        O.write("curve\t%s\t%s\t%d\t%.4f\t\t\t\t\n" % (e, k[0], k[1], v))


def cross(C, ax, t):
    pts = [(i, C[(ax, i)]) for i in range(1, 41) if (ax, i) in C]
    pts.sort()
    for a in range(len(pts) - 1):
        x1, y1 = pts[a]
        x2, y2 = pts[a + 1]
        if (y1 - t) * (y2 - t) <= 0 and y1 != y2:
            return x1 + (t - y1) / (y2 - y1) * (x2 - x1)
    return None


print("\n[A-3] 요소별 교차 지수 — 독립 분류인데도 같은 자리에서 전이하는가")
print("%-6s %9s %9s %9s %9s %9s %9s"
      % ("요소", "근위25%", "근위50%", "근위75%", "원위75%", "원위50%", "원위25%"))
R = {}
for e in EL:
    C = CURVE[e]
    vals = [cross(C, "prox", t) for t in (0.25, 0.5, 0.75)] + \
           [cross(C, "dist", t) for t in (0.75, 0.5, 0.25)]
    R[e] = vals
    print("%-6s %s" % (e, " ".join("%9.2f" % v if v is not None else "        -" for v in vals)))
    O.write("cross\t%s\t%s\t\t\t\t\t\t\n" % (e, ",".join("%.2f" % v if v is not None else ""
                                                        for v in vals)))
for j, lab in enumerate(["근위25%", "근위50%", "근위75%", "원위75%", "원위50%", "원위25%"]):
    v = [R[e][j] for e in EL if R[e][j] is not None]
    if len(v) >= 4:
        print("  %s : 범위 %.2f ~ %.2f (폭 %.2f 유닛) · 중앙 %.2f · 요소 %d"
              % (lab, min(v), max(v), max(v) - min(v), statistics.median(v), len(v)))
O.close()
print("\nsaved: out/SUBUNIT_BIMODAL2.tsv")
