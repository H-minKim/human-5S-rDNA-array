# -*- coding: utf-8 -*-
# array 서열 자체에서 반복 주기를 직접 잰다 — 2 판
#   1 판의 결함 : 주기 하한이 40 이라 단순반복(44~86 bp)이 뽑혔다. 5S 유닛은 유전자 119 bp 를 포함하므로
#                 어떤 진짜 유닛도 120 bp 미만일 수 없다 -> 하한 130
#   방법 A : 한-핫 4채널 FFT 자기상관 -> lag 마다 글자 일치율 (참조·주석 안 씀)
#            기본 주기 = 봉우리 중 배수(2p·3p)도 봉우리인 것들 가운데 가장 작은 lag
#   방법 B : 3 회 이상 나오는 25-mer 의 연속 출현 간격 최빈값 (A 와 독립)
#   방법 C : 고른 주기 p 로 실제 반복 비율 — 50 bp 창이 p 만큼 떨어진 자리와 45/50 이상 맞는 위치의 비율
#            (진짜 탠덤이면 거의 1, 유전자만 흩어져 있으면 0 에 가깝다)
#   산출 : PERIOD_DIRECT2.tsv
# MIT License · Copyright (c) 2026 Hungmin Kim
import io
import os
import sys
import glob
import collections
import numpy as np

NL = chr(10)
TAB = chr(9)
CAN = "$RDNA_OTHERSP_DIR/data_v2/01_canonical"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lit")
WIN = 150000
MAXLAG = 20000
MINLAG = 130       # 5S 유전자 119 bp -> 진짜 유닛은 이보다 크다
W50 = 50


def load(p):
    L = io.open(p, encoding="utf8")
    hd = L.readline().strip()
    seq = "".join(x.strip() for x in L).upper()
    kv = dict(x.split("=", 1) for x in hd.split() if "=" in x)
    a, b = kv["array"].split("-")
    return seq, int(a), int(b)


def autocorr(s):
    n = len(s)
    m = 1
    while m < 2 * n:
        m *= 2
    arr = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
    tot = np.zeros(m // 2 + 1, dtype=np.complex128)
    for ch in b"ACGT":
        v = (arr == ch).astype(np.float64)
        F = np.fft.rfft(v, m)
        tot += F * np.conj(F)
    c = np.fft.irfft(tot, m)[:n]
    return c / np.maximum(n - np.arange(n), 1)


def tandem_frac(s, p):
    # p 만큼 떨어진 50 bp 창이 45 자리 이상 맞는 위치의 비율
    a = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
    if len(a) <= p + W50 + 1:
        return 0.0, 0
    eq = (a[:-p] == a[p:]).astype(np.int32)
    cs = np.concatenate(([0], np.cumsum(eq)))
    n = len(eq) - W50 + 1
    if n <= 0:
        return 0.0, 0
    win = cs[W50:W50 + n] - cs[0:n]
    return float((win >= 45).mean()), n


def kmer_gap(s, k=25):
    P = collections.defaultdict(list)
    for i in range(0, len(s) - k + 1):
        P[s[i:i + k]].append(i)
    g = collections.Counter()
    nk = 0
    for key, v in P.items():
        if len(v) < 3 or "N" in key:
            continue
        nk += 1
        for j in range(len(v) - 1):
            d = v[j + 1] - v[j]
            if MINLAG <= d <= MAXLAG:
                g[int(round(d / 10.0)) * 10] += 1
    return g, nk


ST = {}
f = io.open("$RDNA_OTHERSP_DIR/data_v2/00_catalog/5S_array_structure.tsv", encoding="utf8")
hd = f.readline().rstrip(NL).split(TAB)
ix = dict((n, i) for i, n in enumerate(hd))
for L in f:
    r = L.rstrip(NL).split(TAB)
    ST["%s__%s_%s_%s" % (r[ix["species"]], r[ix["chrom"]], r[ix["array_start"]], r[ix["array_end"]])] = r
UL = {}
for p in glob.glob(os.path.join(CAN, "*", "*.unit.fa")):
    aid = os.path.basename(p).replace(".unit.fa", "")
    UL[aid] = len("".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:]))

w = io.open(os.path.join(OUT, "PERIOD_DIRECT2.tsv"), "w", encoding="utf8", newline=NL)
w.write(TAB.join(["array_id", "array_span", "measured_bp", "my_unit_len", "cms_period_med", "copies_cms",
                  "A_period", "A_ident", "A_base", "A_signal", "A_harmonics", "B_period", "B_frac",
                  "C_tandem_frac", "C_at_my_unit", "copies_from_A", "verdict"]) + NL)
files = [p for p in sorted(glob.glob(os.path.join(CAN, "*", "*.fa"))) if not p.endswith(".unit.fa")]
print("%-50s %7s %6s %6s %6s %6s %6s %5s  %s" % ("array_id", "span", "unit", "A", "B", "A일치", "C반복", "copyA", "판정"))
res = []
for p0 in files:
    aid = os.path.basename(p0).replace(".fa", "")
    seq, a0, a1 = load(p0)
    span = a1 - a0 + 1
    fl = (len(seq) - span) // 2
    arr = seq[fl:fl + span]
    sub = arr if len(arr) <= WIN else arr[(len(arr) - WIN) // 2:(len(arr) - WIN) // 2 + WIN]
    hi = min(MAXLAG, max(len(sub) // 3, MINLAG + 5))
    idn = autocorr(sub)
    base = float(np.median(idn[MINLAG:hi])) if hi > MINLAG + 4 else 0.0
    pk = []
    seg = idn[MINLAG:hi]
    mx = float(seg.max()) if len(seg) else 0.0
    thr = max(base + 0.05, base + 0.35 * (mx - base))
    for i in range(2, len(seg) - 2):
        v = seg[i]
        if v >= thr and v >= seg[i - 1] and v >= seg[i + 1] and v > seg[i - 2] and v > seg[i + 2]:
            pk.append((MINLAG + i, float(v)))
    # 배수도 봉우리인 것 우선, 그 중 가장 작은 lag
    def harmof(q):
        h = 0
        for k in (2, 3, 4):
            d = q * k
            if d < hi and any(abs(x[0] - d) <= max(3, q // 50) for x in pk):
                h += 1
        return h
    Ap, Ai, harm = 0, 0.0, 0
    if pk:
        strong = [q for q in pk if q[1] >= base + 0.80 * (mx - base)]
        withh = [q for q in strong if harmof(q[0]) >= 2]
        pick = min(withh, key=lambda x: x[0]) if withh else min(strong, key=lambda x: x[0])
        Ap, Ai = pick[0], pick[1]
        harm = harmof(Ap)
    g, nk = kmer_gap(sub)
    Bp, Bf = 0, 0.0
    if g:
        tot = sum(g.values())
        Bp, cnt = g.most_common(1)[0]
        Bf = cnt / float(tot)
    myu = UL.get(aid, 0)
    Cf, _ = tandem_frac(sub, Ap) if Ap else (0.0, 0)
    Cu, _ = tandem_frac(sub, myu) if myu and myu < len(sub) - 60 else (0.0, 0)
    cms = ST[aid][ix["period_med"]] if aid in ST else "NA"
    cpc = ST[aid][ix["copies"]] if aid in ST else "NA"
    cpa = int(round(span / float(Ap))) if Ap else 0
    sig = Ai - base
    near = (Ap and myu and abs(Ap - myu) <= max(0.08 * myu, 20))
    AB = (Ap and Bp and abs(Ap - Bp) <= max(0.10 * Ap, 20))
    if not Ap or sig < 0.05 or Cf < 0.10:
        v = "탠덤 신호 없음"
    elif near and AB and Cf >= 0.5:
        v = "OK 직접 확인"
    elif near and AB:
        v = "OK (반복 비율 %.2f)" % Cf
    elif AB:
        v = "주기 %d != 내 유닛 %d" % (Ap, myu)
    else:
        v = "불확실 (A %d · B %d)" % (Ap, Bp)
    res.append((aid, v))
    print("%-50s %7d %6d %6d %6d %6.3f %6.2f %5d  %s" % (aid, span, myu, Ap, Bp, Ai, Cf, cpa, v))
    w.write(TAB.join([aid, str(span), str(len(sub)), str(myu), cms, cpc, str(Ap), "%.3f" % Ai, "%.3f" % base,
                      "%.3f" % sig, str(harm), str(Bp), "%.3f" % Bf, "%.3f" % Cf, "%.3f" % Cu,
                      str(cpa), v]) + NL)
    sys.stdout.flush()
w.close()
print()
print("판정 집계:", dict(collections.Counter(v.split(" (")[0].split(" !=")[0] for _, v in res)))
print("표 : lit/PERIOD_DIRECT2.tsv")
