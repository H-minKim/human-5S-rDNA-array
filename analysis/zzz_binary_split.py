# -*- coding: utf-8 -*-
"""메틸을 양분으로 다시 — 컷 0.508 (이봉의 골) (2026-10-12).

왜 바꾸나
  지금까지 쓰던 층은 삼분위였고 경계가 0.7034 / 0.8485 였다.
  그런데 유닛 CpG β 분포는 뚜렷한 이봉이다 :
     낮은 봉우리 0.125~0.175 · 골 0.475~0.500 (가장 얕은 칸 n=148) · 높은 봉우리 0.850~0.875
  0.7034 는 높은 봉우리의 오르막 중턱을 자른다 - 같은 집단을 반으로 쪼갠 것이다.
  골에서 자르는 것이 옳다. 이 프로젝트에 기록된 컷 0.508 이 바로 그 자리다
  (경험적 최저 칸 0.475~0.500 · 0.488 과 0.508 은 나누는 비율이 25.0% ↔ 25.3% 로 거의 같다).

내는 것
  fhmm_peaks_5bp_bin.npz : FiberHMM 5 bp 칸 x 자리 · 2 층 (저 β<0.508 · 고 β>=0.508)
  occlusion_bin.npz      : flank-core m6A 빈도비용 nA / nM · 2 층
옛 산출(삼분위)은 지우지 않는다.
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
CUT = 0.508
MLTHR = 128
PER, HALF = 2240, 1120
NB = 200
SM = 11
CLS = [(30, 34, "30~34"), (35, 39, "35~39"), (40, 44, "40~44"), (45, 49, "45~49")]
NC = len(CLS)
SRC = "PACBIO_FIBERSEQ"
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
LO, HI = -200, 200
NR = HI - LO + 1

BETA = {(s, r, sl, a): v for (s, sr, r, sl, a), v in
        L.load_beta(min_cpg=10).items() if sr == SRC}
bb = np.array(list(BETA.values()))
print("β n=%s · 컷 %.3f -> 저 %d (%.1f%%) · 고 %d (%.1f%%)"
      % (format(len(bb), ","), CUT, (bb < CUT).sum(), 100 * (bb < CUT).mean(),
         (bb >= CUT).sum(), 100 * (bb >= CUT).mean()))

# ══ 1. FiberHMM 5 bp 칸 · 2 층
SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio" or x[5] == "msp":
                continue
            SEG[(x[0], x[4], x[2], x[3])].append((int(x[6]), int(x[7]), int(x[8])))

OBS = np.zeros((2, NC, PER))
NUL = np.zeros((2, NC, PER))
NMOL = np.zeros(2)
rng = np.random.default_rng(20261012)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_READCOV_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio" or int(x[5]) > -HALF or int(x[6]) < HALF - 1:
                continue
            k = (x[0], x[4], x[2], x[3])
            b = BETA.get(k)
            if b is None:
                continue
            si = 0 if b < CUT else 1
            cen, cls = [], []
            for (s0, s1, sz) in SEG.get(k, []):
                c = (s0 + s1) // 2
                if not (-HALF <= c < HALF):
                    continue
                for j, (a0, a1, _) in enumerate(CLS):
                    if a0 <= sz <= a1:
                        cen.append(c + HALF)
                        cls.append(j)
                        break
            NMOL[si] += 1
            if not cen:
                continue
            cen = np.array(cen)
            cls = np.array(cls)
            np.add.at(OBS[si], (cls, cen), 1)
            d = rng.integers(0, PER, NB)
            for t in range(NB):
                np.add.at(NUL[si], (cls, (cen + d[t]) % PER), 1)
NUL /= float(NB)
np.savez_compressed(os.path.join(DD, "fhmm_peaks_5bp_bin.npz"),
                    obs=OBS, nul=NUL, nmol=NMOL, per=PER, half=HALF, cut=CUT,
                    cls=np.array([c[2] for c in CLS]))
print("FiberHMM 분자 : 저 %s · 고 %s" % (format(int(NMOL[0]), ","),
                                        format(int(NMOL[1]), ",")))

K = np.ones(SM) / SM


def sm_(v):
    return np.convolve(np.r_[v[-SM:], v, v[:SM]], K, "same")[SM:-SM]


print()
print("## FiberHMM 봉우리 (양분) — 저 자리 기준")
print("  %-8s %8s %10s %10s %9s" % ("크기", "자리", "저 배수", "고 배수", "저/고"))
for j, (a0, a1, lab) in enumerate(CLS):
    o0, e0_ = sm_(OBS[0, j]), sm_(NUL[0, j])
    r0 = np.where(e0_ > 0, o0 / np.maximum(e0_, 1e-9), np.nan)
    i = int(np.nanargmax(r0))
    o1, e1 = sm_(OBS[1, j]), sm_(NUL[1, j])
    r1 = o1[i] / e1[i] if e1[i] > 0 else float("nan")
    print("  %-8s %+8d %10.2f %10.2f %9.2f"
          % (lab, i - HALF, r0[i], r1, r0[i] / r1 if r1 > 0 else float("nan")))

# ══ 2. occlusion (nA / nM) · 2 층
UNIT, _ = L.load_units()
TSSU = {k: (v[2] - v[0] + 1) for k, v in UNIT.items()}
nA = np.zeros((2, NR))
nM = np.zeros((2, NR))
nmolD = np.zeros(2)
for node in ("w3", "iB", "iA"):
    p = os.path.join(L.R, "08_methylation_s5", node, L.METH[SRC],
                     "M6A_READUNIT_%s.tsv.gz" % L.HOST[node])
    if not os.path.exists(p):
        continue
    MOL = collections.defaultdict(list)
    with gzip.open(p, "rt", encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if len(x) < 14 or x[1] != SRC or x[12] != "Y" or not x[13]:
                continue
            smp, slot, ap, rid = x[0], x[4], x[5], x[6]
            if smp in L.PB_EXCL or (smp, slot) in SHIFT:
                continue
            MOL[(smp, rid, slot, ap)].append(x[13])
    for k, blocks in MOL.items():
        b = BETA.get(k)
        tu = TSSU.get((k[0], k[2], k[3]))
        if b is None or tu is None:
            continue
        si = 0 if b < CUT else 1
        nmolD[si] += 1
        for blk in blocks:
            for t in blk.split(","):
                if not t:
                    continue
                a_, b_ = t.split(":")
                r = int(a_) - tu
                if LO <= r <= HI:
                    nA[si, r - LO] += 1
                    nM[si, r - LO] += (int(b_) >= MLTHR)
np.savez_compressed(os.path.join(DD, "occlusion_bin.npz"),
                    nA=nA, nM=nM, nmol=nmolD, lo=LO, hi=HI, cut=CUT)
print()
print("occlusion 분자 : 저 %s · 고 %s" % (format(int(nmolD[0]), ","),
                                         format(int(nmolD[1]), ",")))

rel = np.arange(LO, HI + 1)


def seg(a, c0, c1):
    i0, i1 = c0 - LO, c1 - LO + 1
    return np.nan if (i0 < 0 or i1 > NR) else a[i0:i1].sum()


print()
print("## flank-core m6A 빈도비 (원문 창 배치) — 양분")
print("  %-8s %6s %10s %10s %9s" % ("자리", "pos", "저", "고", "저/고"))
for nm2, pos in (("TFIIIB", -11), ("ICR", 58)):
    v = []
    for si in (0, 1):
        cA = seg(nA[si], pos - 7, pos + 1)
        cM = seg(nM[si], pos - 7, pos + 1)
        fA = np.nansum([seg(nA[si], pos - 35, pos - 27), seg(nA[si], pos + 22, pos + 27)])
        fM = np.nansum([seg(nM[si], pos - 35, pos - 27), seg(nM[si], pos + 22, pos + 27)])
        v.append((fM / fA) / ((cM + 0.5) / (cA + 1.0)))
    print("  %-8s %+6d %10.3f %10.3f %9.2f" % (nm2, pos, v[0], v[1], v[0] / v[1]))
print()
print("저장 -> fhmm_peaks_5bp_bin.npz · occlusion_bin.npz")
