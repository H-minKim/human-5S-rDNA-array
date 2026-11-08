# -*- coding: utf-8 -*-
"""[숙제2] 격자 전체 다중성 문턱 + [숙제1] 구역 창 통일 (2026-10-12).

숙제2 - zzn/zzm 은 문턱을 크기 칸마다 따로 세웠다.
  자리 269 개에 대해서는 family-wise 지만 칸 29 개를 가로질러서는 보정이 없다.
  귀무에서도 칸마다 5% 확률로 봉우리가 하나 나오므로 목록에 평균 1.5 개가 섞인다.
  고침 : 치환 한 번에 이동량을 한 벌 뽑아 모든 칸에 같이 쓰고(= 한 분자의 footprint 는
  어느 칸에 있든 함께 움직인다 - 더 옳은 귀무이기도 하다), 그 치환의 격자 전체 최대를
  기록해 95 분위로 문턱 하나를 세운다.

숙제1 - 구역 창이 두 벌이었다.
  Fig5E/Fig5D 는 FWHM 으로 뽑은 narrow(-28..-10 / +57..+68),
  그림1/표2-4 와 스캔은 내가 눈대중으로 넓힌 wide(-35..0 / +45..+90).
  순위는 안 바뀌고 이름표만 바뀌지만 "상위 10 중 9 ↔ 6" 으로 대표 문장이 달라진다.
  -> narrow 로 통일한다 (자료에서 수로 유도된 쪽). wide 는 표에 한 열 더 남긴다.

zzn 의 산출(scan_maps.npz · peaks_both.tsv)은 지우지 않는다. 여기는 _v2 로 낸다.
"""
import io
import os
import sys
import gzip
import glob
import collections
import numpy as np

_SELF = os.path.basename(os.path.abspath(__file__))


def _need(k):
    # No defaults are set: if a variable is missing the script stops with a readable message
    v = os.environ.get(k, "")
    if not v:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s." % (k, _SELF)
            + chr(10)
            + "The deposited copy has no default; see README.md section 2")
    return v


sys.path.insert(0, _need("RDNA_LIB_DIR"))
import a_lib as L

NL = chr(10)
TB = chr(9)
DD = _need("RDNA_DATA_ROOT")
FD = os.path.join(DD, "fhmm")
OUTD = _need("RDNA_OUT_DIR")
MLTHR = 128
PER, HALF = 2240, 1120
LO, HI = -150, 118
NP = HI - LO + 1
NB = 200
SM = 11
SEP = 15
MERGE = 12
SBIN, SMIN, SMAX = 5, 20, 164
GAP, FLANKW = 8, 17
WIDTHS = [9, 15, 21, 27, 33, 41, 51, 61, 71]
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
ZT, ZI = (-28, -10), (57, 68)          # 정본 = FWHM 으로 뽑은 narrow
ZTW, ZIW = (-35, 0), (45, 90)          # 옛 wide (표에 한 열로만 남긴다)
rel = np.arange(LO, HI + 1)
rng = np.random.default_rng(20261011)
K = np.ones(SM) / SM

# ══ FiberHMM
SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio" or x[5] == "msp":
                continue
            SEG[(x[0], x[4], x[2], x[3])].append(((int(x[6]) + int(x[7])) // 2, int(x[8])))
BINS = [(s, s + SBIN - 1) for s in range(SMIN, SMAX + 1, SBIN)]
MID = {b: ([], []) for b in BINS}
nmolF = 0
for p in sorted(glob.glob(os.path.join(FD, "FHMM_READCOV_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio" or int(x[5]) > -HALF or int(x[6]) < HALF - 1:
                continue
            for (c, sz) in SEG.get((x[0], x[4], x[2], x[3]), []):
                if not (-HALF <= c < HALF):
                    continue
                k = (sz - SMIN) // SBIN
                if 0 <= k < len(BINS) and BINS[k][0] <= sz <= BINS[k][1]:
                    MID[BINS[k]][0].append(nmolF)
                    MID[BINS[k]][1].append(c + HALF)
            nmolF += 1
i0, i1 = LO + HALF, HI + HALF + 1
BOK = [b for b in BINS if len(MID[b][1]) >= 150]
MIa = [np.array(MID[b][0], np.int64) for b in BOK]
CEa = [np.array(MID[b][1], np.int64) for b in BOK]
NBK = len(BOK)


def cv(ce):
    return np.convolve(1000.0 * np.bincount(ce, minlength=PER)[i0:i1] / nmolF, K, "same")


OBSM = np.array([cv(c) for c in CEa])
PROF = np.empty((NB, NBK, i1 - i0))
for t in range(NB):
    d = rng.integers(0, PER, nmolF)            # 한 벌을 모든 칸에 같이 쓴다
    for bi in range(NBK):
        PROF[t, bi] = cv((CEa[bi] + d[MIa[bi]]) % PER)
BASE = PROF.mean(0)
ENR = OBSM / np.maximum(BASE, 1e-9)
NR = PROF / np.maximum(BASE, 1e-9)
FTHR_BIN = np.percentile(NR.max(2), 95, axis=0)          # 옛 방식 (칸마다)
FTHRG = float(np.percentile(NR.reshape(NB, -1).max(1), 95))   # 격자 전체
print("FiberHMM 문턱 : 칸마다 %.2f~%.2f  ->  격자 전체 %.2f"
      % (FTHR_BIN.min(), FTHR_BIN.max(), FTHRG))
FH, FH_OLD = [], []
FMAP, FTHR, FSZ = [], [], []
for bi, b in enumerate(BOK):
    enr = ENR[bi]
    FMAP.append(enr); FTHR.append(FTHRG); FSZ.append([b[0], b[1]])
    for thr, box in ((FTHRG, FH), (FTHR_BIN[bi], FH_OLD)):
        pk = []
        for j in np.argsort(enr)[::-1]:
            if enr[j] <= thr:
                break
            if any(abs(rel[j] - q) < SEP for q, _ in pk):
                continue
            pk.append((rel[j], enr[j]))
        for q, e in pk:
            box.append((q, b, e))
print("FiberHMM 봉우리 : 칸마다 문턱 %d 개  ->  격자 전체 문턱 %d 개"
      % (len(FH_OLD), len(FH)))
print("FiberHMM 분자 %s · 통과 봉우리 %d" % (format(nmolF, ","), len(FH)))

# ══ Dubocanin
UNIT, _ = L.load_units()
TSSU = {k: (v[2] - v[0] + 1) for k, v in UNIT.items()}
A, M = [], []
for node in ("w3", "iB", "iA"):
    p = os.path.join(L.R, "08_methylation_s5", node, L.METH["PACBIO_FIBERSEQ"],
                     "M6A_READUNIT_%s.tsv.gz" % L.HOST[node])
    if not os.path.exists(p):
        continue
    mol = collections.defaultdict(list)
    with gzip.open(p, "rt", encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if len(x) < 14 or x[1] != "PACBIO_FIBERSEQ" or x[12] != "Y" or not x[13]:
                continue
            smp, slot, ap, rid = x[0], x[4], x[5], x[6]
            if smp in L.PB_EXCL or (smp, slot) in SHIFT:
                continue
            mol[(smp, rid, slot, ap)].append(x[13])
    for k, blocks in mol.items():
        tu = TSSU.get((k[0], k[2], k[3]))
        if tu is None:
            continue
        a = np.zeros(NP, np.int32); m = np.zeros(NP, np.int32)
        for blk in blocks:
            for t in blk.split(","):
                if not t:
                    continue
                p0, v0 = t.split(":")
                r = int(p0) - tu
                if LO <= r <= HI:
                    a[r - LO] += 1
                    if int(v0) >= MLTHR:
                        m[r - LO] += 1
        if a.sum() >= 20:
            A.append(a); M.append(m)
A = np.array(A, np.int64); M = np.array(M, np.int64)
nmolD = len(A)


def W(cs, c0, c1):
    o = np.full(NP, np.nan); i = np.arange(NP)
    x0, x1 = i + c0, i + c1 + 1
    ok = (x0 >= 0) & (x1 <= NP)
    o[ok] = cs[x1[ok]] - cs[x0[ok]]
    return o


def dprof(sa, sm, w):
    h = w // 2
    ca = np.concatenate([[0.0], np.cumsum(sa)]); cm = np.concatenate([[0.0], np.cumsum(sm)])
    cA = W(ca, -h, h); cM = W(cm, -h, h)
    fA = W(ca, -h - GAP - FLANKW, -h - GAP - 1) + W(ca, h + GAP + 1, h + GAP + FLANKW)
    fM = W(cm, -h - GAP - FLANKW, -h - GAP - 1) + W(cm, h + GAP + 1, h + GAP + FLANKW)
    s = np.full(NP, np.nan)
    ok = np.isfinite(cA) & np.isfinite(fA) & (cA >= 200) & (fA >= 200)
    s[ok] = (fM[ok] / fA[ok]) / ((cM[ok] + 0.5) / (cA[ok] + 1.0))
    return s


NULA = np.empty((NB, NP)); NULM = np.empty((NB, NP))
for t in range(NB):
    d = rng.integers(0, NP, nmolD)
    ra = np.zeros(NP); rm = np.zeros(NP)
    for o in range(NP):
        sel = d == o
        if sel.any():
            ra += np.roll(A[sel].sum(0), o)
            rm += np.roll(M[sel].sum(0), o)
    NULA[t] = ra; NULM[t] = rm
sA, sM = A.sum(0).astype(float), M.sum(0).astype(float)
DNULL = np.array([[np.nanmax(dprof(NULA[t], NULM[t], w)) for w in WIDTHS]
                  for t in range(NB)])
DTHR_W = np.nanpercentile(DNULL, 95, axis=0)                 # 옛 방식 (폭마다)
DTHRG = float(np.nanpercentile(np.nanmax(DNULL, axis=1), 95))  # 격자 전체
print("Dubocanin 문턱 : 폭마다 %.3f~%.3f  ->  격자 전체 %.3f"
      % (DTHR_W.min(), DTHR_W.max(), DTHRG))
DU, DU_OLD = [], []
DMAP, DTHR = [], []
print()
print("=" * 92)
print("## Dubocanin — core 폭마다 귀무 통과 봉우리 전부")
for wi, w in enumerate(WIDTHS):
    obs = dprof(sA, sM, w)
    DMAP.append(np.nan_to_num(obs, nan=0.0)); DTHR.append(DTHRG)
    for thr, box in ((DTHRG, DU), (DTHR_W[wi], DU_OLD)):
        pk = []
        for j in np.argsort(np.where(np.isfinite(obs), obs, -9))[::-1]:
            if not np.isfinite(obs[j]) or obs[j] <= thr:
                break
            if any(abs(rel[j] - q) < SEP for q, _ in pk):
                continue
            pk.append((rel[j], obs[j]))
        for q, v in pk:
            box.append((q, w, v))
    pkg = [(q, v) for q, w2, v in DU if w2 == w]
    print("  core %-3d : %s" % (w, "  ".join("%+d:%.2f" % t for t in pkg) or "없음"))
print("Dubocanin 봉우리 : 폭마다 문턱 %d 개  ->  격자 전체 문턱 %d 개"
      % (len(DU_OLD), len(DU)))

# ── 봉우리 목록을 남겨 둔다 (재계산 방지)
with io.open(os.path.join(DD, "peaks_both_v2.tsv"), "w", encoding="utf8", newline=NL) as g:
    g.write(TB.join(["방법", "자리", "크기또는폭_lo", "크기또는폭_hi", "값",
                     "이름표_narrow", "이름표_wide"]) + NL)

    def _lab(q, zt, zi):
        return "TFIIIB" if zt[0] <= q <= zt[1] else ("ICR" if zi[0] <= q <= zi[1] else "")

    for q, b, e in FH:
        g.write(TB.join(["FiberHMM", "%d" % q, "%d" % b[0], "%d" % b[1], "%.4f" % e,
                         _lab(q, ZT, ZI), _lab(q, ZTW, ZIW)]) + NL)
    for q, w, v in DU:
        g.write(TB.join(["Dubocanin", "%d" % q, "%d" % w, "%d" % w, "%.4f" % v,
                         _lab(q, ZT, ZI), _lab(q, ZTW, ZIW)]) + NL)

# ══ 표3 : 자리로 뭉치기 — 단일연결은 사슬처럼 다 붙는다. 밀도 봉우리에 붙인다.
allq = np.array([q for q, _, _ in FH] + [q for q, _, _ in DU])
h = np.zeros(NP)
for q in allq:
    h[q - LO] += 1
hs = np.convolve(h, np.ones(11) / 11.0, "same")
modes = []
for j in range(1, NP - 1):
    if hs[j] >= hs[j - 1] and hs[j] > hs[j + 1] and hs[j] > 0:
        if not modes or (rel[j] - modes[-1]) >= 20:
            modes.append(rel[j])
        elif hs[j] > hs[modes[-1] - LO]:
            modes[-1] = rel[j]
cent = [[] for _ in modes]
extra = []
for q in sorted(set(allq.tolist())):
    d = [abs(q - m) for m in modes]
    k = int(np.argmin(d)) if modes else -1
    if modes and d[k] <= 20:
        cent[k].append(q)
    else:
        extra.append([q])
cent = [c for c in cent if c]
# 낱개는 이미 군집 안에 든 자리면 버린다 (중복 행 방지)
for e in extra:
    q = e[0]
    if not any(min(c) - MERGE <= q <= max(c) + MERGE for c in cent):
        cent.append(e)
cent.sort(key=lambda c: min(c))

np.savez_compressed(os.path.join(DD, "scan_maps_v2.npz"),
                    fmap=np.array(FMAP), fthr=np.array(FTHR), fsz=np.array(FSZ),
                    dmap=np.array(DMAP), dthr=np.array(DTHR), dw=np.array(WIDTHS),
                    lo=LO, hi=HI)
print("지도 저장 -> scan_maps_v2.npz")
print()
print("=" * 100)
print("## 표3 — 자리별 두 방법 통과 (자리 뭉침 ±%d bp)" % MERGE)
print("  %-14s | %-34s | %-26s | %s" % ("자리", "FiberHMM (크기:배수)", "Dubocanin (폭:값)", "판정"))
ROWS = []
for g in cent:
    lo, hi = min(g), max(g)
    f = sorted([(b[0], b[1], e) for q, b, e in FH if lo - 1 <= q <= hi + 1],
               key=lambda t: -t[2])
    d = sorted([(w, v) for q, w, v in DU if lo - 1 <= q <= hi + 1], key=lambda t: -t[1])
    if not f and not d:
        continue
    jg = "둘 다" if (f and d) else ("FiberHMM 만" if f else "Dubocanin 만")
    fs = " ".join("%d-%d:%.1f" % t for t in f[:4])
    ds = " ".join("%d:%.2f" % t for t in d[:4])
    print("  %-14s | %-34s | %-26s | %s"
          % ("%+d..%+d" % (lo, hi) if lo != hi else "%+d" % lo, fs[:34], ds[:26], jg))
    ROWS.append((lo, hi, len(f), len(d), f[0][2] if f else 0, d[0][1] if d else 0, jg))

out = os.path.join(OUTD, "TF_20261011_table3_two_methods_by_position.tsv")
with io.open(out, "w", encoding="utf8", newline=NL) as g:
    g.write("# 자리 뭉침 ±%d bp · FiberHMM 배수문턱 · Dubocanin 값문턱 · 둘 다 귀무 %d 회" % (MERGE, NB) + NL)
    g.write(TB.join(["자리lo", "자리hi", "FiberHMM_통과칸수", "Dubocanin_통과폭수",
                     "FiberHMM_최대배수", "Dubocanin_최대값", "판정"]) + NL)
    for r in ROWS:
        g.write(TB.join(["%d" % r[0], "%d" % r[1], "%d" % r[2], "%d" % r[3],
                         "%.3f" % r[4], "%.3f" % r[5], r[6]]) + NL)
print()
print("적었다 ->", out)
