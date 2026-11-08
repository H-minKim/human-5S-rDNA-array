# -*- coding: utf-8 -*-
# 흩어진 5S 자리 : ① 덮는 유닛 부분으로 묶기 ② 자리별 이봉/단봉 ③ 상태별 유닛좌표 프로파일 vs array
import io, os, math, collections, statistics

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


D = os.path.join(_need("RDNA_WORK_ROOT"), "dispersed_vs_array")
R = os.path.join(_need("RDNA_VERSION_ROOT"), "05_checks", "2026-09-17_TUC_dropped_reads_5S_content_methylation")
EL = [("AluY", 1, 318), ("NTS1", 319, 1178), ("GA", 1179, 1551), ("NTS2", 1552, 2031),
      ("prom", 2032, 2121), ("gene", 2122, 2240)]
NB = 45


def sub_of_bin(b):
    c = b * 50 + 25
    return next(n for n, a, z in EL if a <= c <= z)


BSUB = [sub_of_bin(b) for b in range(NB)]
MIN_LOCI, MIN_SMP = 200, 10

# ---------- 적재 ----------
INFO = {}
for i, L in enumerate(io.open(os.path.join(D, "DISP_SITE_INFO_ALL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    cov = {int(a): int(b) for a, b in (x.split(":") for x in f[4].split(";") if x)}
    INFO[f[0]] = dict(n=int(f[1]), groups=f[2], states=f[3], cov=cov)
SITES = {}
for i, L in enumerate(io.open(os.path.join(R, "CHM13_positions", "SITES.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    SITES[f[0]] = dict(chr=f[1], start=int(f[2]), end=int(f[3]), nsmp=int(f[4]), main=f[6], cat=f[18])
BIMO = {}
for i, L in enumerate(io.open(os.path.join(R, "bimodal_intact_units", "SITE_BIMODAL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    f = L.rstrip("\n").split("\t")
    BIMO[f[0]] = dict(dip=f[8], D=f[14], call=f[15])
BINS = collections.defaultdict(lambda: [0, 0])       # (site,state,half,bin)
for i, L in enumerate(io.open(os.path.join(D, "DISP_SITE_BINS_ALL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    s, st, h, b, n, k = L.rstrip("\n").split("\t")
    BINS[(s, st, int(h), int(b))] = [int(n), int(k)]
LB = {}
for i, L in enumerate(io.open(os.path.join(D, "DISP_SITE_LOCUSBETA_ALL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    s, v = L.rstrip("\n").split("\t")
    LB[s] = [float(x) for x in v.split(",") if x]
SB = collections.defaultdict(dict)
for i, L in enumerate(io.open(os.path.join(D, "DISP_SITE_SAMPLE_ALL.tsv"), encoding="utf8")):
    if i == 0:
        continue
    s, smp, n, k = L.rstrip("\n").split("\t")
    if int(n) >= 20:
        SB[s][smp] = int(k) / float(n)
ARR = collections.defaultdict(lambda: [0, 0])        # (pos, state, bin)
for i, L in enumerate(io.open(os.path.join(D, "ARRAY_V11_BINS.tsv"), encoding="utf8")):
    if i == 0:
        continue
    p, st, b, n, k = L.rstrip("\n").split("\t")
    ARR[(p, st, int(b))] = [int(n), int(k)]


def prof(get, st, halves=(0, 1), min_n=30):
    out = []
    for b in range(NB):
        n = sum(get(st, h, b)[0] for h in halves)
        k = sum(get(st, h, b)[1] for h in halves)
        out.append(k / float(n) if n >= min_n else None)
    return out


def site_get(s):
    return lambda st, h, b: BINS.get((s, st, h, b), [0, 0])


def group_get(sites):
    def g(st, h, b):
        n = k = 0
        for s in sites:
            v = BINS.get((s, st, h, b))
            if v:
                n += v[0]
                k += v[1]
        return [n, k]
    return g


def arr_get(pos):
    return lambda st, h, b: ARR.get((pos, st, b), [0, 0]) if h == 0 else [0, 0]


def sub_means(p):
    out = {}
    for sname, a, z in EL:
        v = [p[b] for b in range(NB) if BSUB[b] == sname and p[b] is not None]
        out[sname] = statistics.mean(v) if v else None
    return out


def corr(a, b):
    xy = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
    if len(xy) < 5:
        return None
    xs, ys = [x for x, _ in xy], [y for _, y in xy]
    if statistics.pstdev(xs) < 1e-9 or statistics.pstdev(ys) < 1e-9:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    return sum((x - mx) * (y - my) for x, y in xy) / math.sqrt(
        sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))


# ---------- bimodality (GMM + BC) ----------
def gmm2(x):
    xs = sorted(x)
    m0, m1 = xs[len(xs) // 5], xs[4 * len(xs) // 5]
    if m1 - m0 < 0.05:
        m0, m1 = max(0.0, m0 - 0.1), min(1.0, m1 + 0.1)
    s0 = s1 = 0.1
    w0 = 0.5
    for _ in range(200):
        a0 = a1 = b0 = b1 = c0 = c1 = 0.0
        for v in x:
            p = w0 * math.exp(-(v - m0) ** 2 / (2 * s0 * s0)) / s0
            q = (1 - w0) * math.exp(-(v - m1) ** 2 / (2 * s1 * s1)) / s1
            r = p / (p + q) if p + q > 0 else 0.5
            a0 += r; b0 += r * v; c0 += r * v * v
            a1 += 1 - r; b1 += (1 - r) * v; c1 += (1 - r) * v * v
        if a0 < 1e-6 or a1 < 1e-6:
            break
        m0n, m1n = b0 / a0, b1 / a1
        s0 = max(0.02, math.sqrt(max(c0 / a0 - m0n ** 2, 1e-6)))
        s1 = max(0.02, math.sqrt(max(c1 / a1 - m1n ** 2, 1e-6)))
        w0n = a0 / len(x)
        conv = abs(m0n - m0) < 1e-6 and abs(m1n - m1) < 1e-6
        m0, m1, w0 = m0n, m1n, w0n
        if conv:
            break
    if m0 > m1:
        m0, m1, s0, s1, w0 = m1, m0, s1, s0, 1 - w0
    return w0, m0, s0, m1, s1


def bc(x):
    n = len(x)
    m = statistics.mean(x)
    s = statistics.pstdev(x)
    if s < 1e-9:
        return 0.0
    g = sum((v - m) ** 3 for v in x) / n / s ** 3
    k = sum((v - m) ** 4 for v in x) / n / s ** 4 - 3
    return (g * g + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3)))


def bimodal(x):
    if len(x) < 50:
        return None
    xx = x if len(x) <= 4000 else x[:: max(1, len(x) // 4000)]
    w, m0, s0, m1, s1 = gmm2(xx)
    D_ = abs(m1 - m0) * math.sqrt(2) / math.sqrt(s0 ** 2 + s1 ** 2)
    B = bc(xx)
    call = "이봉" if (D_ > 2 and B > 0.555 and min(w, 1 - w) >= 0.05) else "단봉"
    return dict(w=w, m0=m0, m1=m1, D=D_, BC=B, call=call)


# ---------- ① 덮는 부분으로 묶기 ----------
O = io.open(os.path.join(D, "DISP_SITES_TABLE.tsv"), "w", encoding="utf8", newline="\n")
O.write("site\tchr\tstart\tmain_group\tn_loci\tn_samples\tcover_pattern\t" +
        "\t".join("cov_" + e[0] for e in EL) + "\tmean_beta\tfrac_low\tfrac_high\t"
        "gmm_m0\tgmm_m1\tgmm_w0\tashmanD\tBC\tcall\ttheir_call\tsample_sd\t" +
        "\t".join("LOW_" + e[0] for e in EL) + "\t" + "\t".join("HIGH_" + e[0] for e in EL) + "\n")
KEEP = []
PAT = {}
for s, inf in INFO.items():
    if s == "ARRAY":
        continue
    sm = SITES.get(s, {})
    if inf["n"] < MIN_LOCI or len(SB.get(s, {})) < MIN_SMP:
        continue
    covf = {}
    for sname, a, z in EL:
        bs = [b for b in range(NB) if BSUB[b] == sname]
        covf[sname] = max(inf["cov"].get(b, 0) for b in bs) / float(inf["n"])
    pat = "+".join(e[0] for e in EL if covf[e[0]] >= 0.5) or "(짧음)"
    PAT[s] = pat
    KEEP.append(s)
    lb = LB.get(s, [])
    bm = bimodal(lb) if lb else None
    lo = sub_means(prof(site_get(s), "low"))
    hi = sub_means(prof(site_get(s), "high"))
    sv = list(SB[s].values())
    f_ = lambda v: "%.3f" % v if v is not None else ""
    O.write("%s\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
        s, sm.get("chr", ""), sm.get("start", 0), sm.get("main", ""), inf["n"], len(SB[s]), pat,
        "\t".join("%.2f" % covf[e[0]] for e in EL),
        f_(statistics.mean(lb)) if lb else "",
        f_(sum(1 for v in lb if v < 0.3) / float(len(lb))) if lb else "",
        f_(sum(1 for v in lb if v > 0.7) / float(len(lb))) if lb else "",
        f_(bm["m0"]) if bm else "", f_(bm["m1"]) if bm else "", f_(bm["w"]) if bm else "",
        f_(bm["D"]) if bm else "", f_(bm["BC"]) if bm else "", bm["call"] if bm else "",
        BIMO.get(s, {}).get("call", ""), f_(statistics.pstdev(sv)) if len(sv) > 1 else "",
        "\t".join(f_(lo[e[0]]) for e in EL), "\t".join(f_(hi[e[0]]) for e in EL)))
O.close()
print("자리 %d 개 (좌위>=%d · 개체>=%d)" % (len(KEEP), MIN_LOCI, MIN_SMP))
pc = collections.Counter(PAT[s] for s in KEEP)
print("\n[① 덮는 유닛 부분] 패턴별 자리 수 · 좌위 수 · 주 구조그룹")
PG = collections.defaultdict(list)
for s in KEEP:
    PG[PAT[s]].append(s)
for p_, ss in sorted(PG.items(), key=lambda kv: -sum(INFO[s]["n"] for s in kv[1])):
    mg = collections.Counter(SITES.get(s, {}).get("main", "?") for s in ss)
    print("  %-34s 자리 %3d · 좌위 %7d · 구조 %s" % (p_, len(ss), sum(INFO[s]["n"] for s in ss),
                                               ", ".join("%s %d" % x for x in mg.most_common(4))))

# ---------- ② 이봉 / 단봉 ----------
print("\n[② 자리별 이봉/단봉] (좌위 β, NTS1+NTS2+코딩 CpG>=10 · GMM Ashman D>2 & BC>0.555 & 작은성분>=5%)")
rows = []
for s in KEEP:
    lb = LB.get(s, [])
    bm = bimodal(lb) if lb else None
    if bm:
        rows.append((PAT[s], s, bm, statistics.mean(lb), len(lb)))
cc = collections.Counter((p_, r["call"]) for p_, s, r, m, n in rows)
for p_ in sorted(set(x[0] for x in rows)):
    print("  %-34s 이봉 %3d · 단봉 %3d" % (p_, cc[(p_, "이봉")], cc[(p_, "단봉")]))
print("\n  이봉 자리 (좌위 수 많은 순)")
for p_, s, r, m, n in sorted([x for x in rows if x[2]["call"] == "이봉"], key=lambda x: -x[4])[:20]:
    sm = SITES.get(s, {})
    print("   %s %s:%.2fM %-10s n=%6d mean %.2f | components %.2f / %.2f (w_lo %.2f) D %.2f BC %.2f | source %s"
          % (s, sm.get("chr"), sm.get("start", 0) / 1e6, sm.get("main"), n, m, r["m0"], r["m1"], r["w"],
             r["D"], r["BC"], BIMO.get(s, {}).get("call", "-")))

# ---------- ③ 상태별 프로파일 vs array ----------
A_ALL = {st: prof(arr_get("ALL"), st, halves=(0,)) for st in ("low", "high")}
A_POS = {(p, st): prof(arr_get(p), st, halves=(0,)) for p in ("distal1_10", "interior", "proximal1_10")
         for st in ("low", "high")}
A_CTL = {st: prof(site_get("ARRAY"), st) for st in ("low", "high")}
print("\n[0] do the two array definitions give the same profile when the state is the same")
for st in ("low", "high"):
    print("  %s : v11 placement against the dispersed ARRAY control, column correlation %.3f · subunit v11 %s | control %s"
          % (st, corr(A_ALL[st], A_CTL[st]),
             " ".join("%s %.2f" % (k, v) for k, v in sub_means(A_ALL[st]).items() if v is not None),
             " ".join("%s %.2f" % (k, v) for k, v in sub_means(A_CTL[st]).items() if v is not None)))

OUTP = io.open(os.path.join(D, "DISP_PROFILES.tsv"), "w", encoding="utf8", newline="\n")
OUTP.write("kind\tname\tstate\t" + "\t".join("b%d" % b for b in range(NB)) + "\n")


def dump(kind, name, st, p):
    OUTP.write("%s\t%s\t%s\t%s\n" % (kind, name, st, "\t".join("%.4f" % v if v is not None else "" for v in p)))


for st in ("low", "high"):
    dump("array_v11", "ALL", st, A_ALL[st])
    dump("array_ctl", "ARRAY", st, A_CTL[st])
    for p in ("distal1_10", "interior", "proximal1_10"):
        dump("array_v11", p, st, A_POS[(p, st)])

print("\n[③] 덮음 패턴(그룹) × 상태 : 하위유닛 평균 β  (흩어진 / array v11 같은 상태)  · 칸 모양 상관 · 반쪽 재현")
for p_, ss in sorted(PG.items(), key=lambda kv: -sum(INFO[s]["n"] for s in kv[1])):
    if sum(INFO[s]["n"] for s in ss) < 2000:
        continue
    g = group_get(ss)
    for st in ("low", "high"):
        pr = prof(g, st)
        if sum(1 for v in pr if v is not None) < 5:
            continue
        h0, h1 = prof(g, st, (0,)), prof(g, st, (1,))
        sm_, sa = sub_means(pr), sub_means(A_ALL[st])
        dump("group", p_, st, pr)
        print("  %-30s %-4s | %s | 모양상관 %s · 반쪽 %s"
              % (p_, st, " ".join("%s %s/%s" % (k, "%.2f" % sm_[k] if sm_[k] is not None else "-",
                                                   "%.2f" % sa[k]) for k, _a, _z in EL if sm_[k] is not None),
                 "%.2f" % corr(pr, A_ALL[st]) if corr(pr, A_ALL[st]) is not None else "-",
                 "%.2f" % corr(h0, h1) if corr(h0, h1) is not None else "-"))
for s in KEEP:
    for st in ("low", "high"):
        pr = prof(site_get(s), st)
        if sum(1 for v in pr if v is not None) >= 5:
            dump("site", s, st, pr)
OUTP.close()

# ---------- ④ 사용자 질문 : AluY 안 기울기 · NTS1 수준 ----------
print("\n[④] AluY 안 50bp 칸 (1-300) 과 NTS1 수준")
def aluy(p):
    return [p[b] for b in range(6)]
for nm, p in (("array v11 LOW", A_ALL["low"]), ("array v11 HIGH", A_ALL["high"]),
              ("array 원위1-10 LOW", A_POS[("distal1_10", "low")]),
              ("array 내부 LOW", A_POS[("interior", "low")])):
    print("  %-22s AluY %s | NTS1 평균 %.3f" % (nm, " ".join("%.2f" % v if v is not None else " -  " for v in aluy(p)),
                                              sub_means(p)["NTS1"]))
for p_, ss in sorted(PG.items(), key=lambda kv: -sum(INFO[s]["n"] for s in kv[1])):
    if "AluY" not in p_ or sum(INFO[s]["n"] for s in ss) < 2000:
        continue
    g = group_get(ss)
    for st in ("low", "high"):
        pr = prof(g, st)
        if pr[0] is None and pr[3] is None:
            continue
        print("  %-22s AluY %s | NTS1 평균 %s" % ((p_ + " " + st)[:22],
              " ".join("%.2f" % v if v is not None else " -  " for v in aluy(pr)),
              "%.3f" % sub_means(pr)["NTS1"] if sub_means(pr)["NTS1"] is not None else "-"))
print("\nsaved DISP_SITES_TABLE.tsv · DISP_PROFILES.tsv")
