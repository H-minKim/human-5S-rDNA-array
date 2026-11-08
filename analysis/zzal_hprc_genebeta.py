# -*- coding: utf-8 -*-
"""HPRC 를 세포주의 층 라벨로 갈라 보면 (2026-10-14 · 2 판).

1 판(`zzak_hprc_cellmethod.py`)에서 두 가지를 배웠다
  ① '전체 풀' 축은 HPRC 에서 안 된다.
     배치 최종이 아닌 유닛(in_final=N)은 FiberHMM 에 3,632 개 있지만
     창을 관통하는 것이 3 개뿐이다 - 짧아서 여과된 것이라 애초에 못 쓴다.
     ⇒ 이 축은 빼고, 그 수를 그림·표에 적는다.
  ② 자를 공개판에 맞춘다. `fhmm_peaks_bin10_both.npz`(표3 의 10 bp 줄)는
     분절에서 msp 만 뺀다(nuc+tf). tf 만 쓰면 49.85, 공개판 규칙이면 52.77 이다.
     ⇒ 공개판 규칙을 따른다.

그래서 이 판이 보는 것은 층 라벨 하나다
    유닛 전체 β 0.508 (본표)  ↔  유전자부(끝 119 bp · CpG>=5) β 0.5 (세포주 방식)

내는 것 : hprc_genebeta.npz
  R   (라벨 2 · 자리 2 · 층 2)  그 자리의 귀무 대비 배수
  N   (라벨 2 · 층 2)           분자 수
  CF  (2 x 2)                   두 라벨의 혼동표
  GC                            유전자부 CpG 수 분포
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

NL, TB = chr(10), chr(9)
DD = _need("RDNA_DATA_ROOT")
FD = os.path.join(DD, "fhmm")
SRC = "PACBIO_FIBERSEQ"
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
PER, HALF = 2240, 1120
NB, SM, MLTHR = 200, 11, 128
GENE_FROM, MIN_GENE_CPG = 2122, 5
CUT_UNIT, CUT_GENE = 0.508, 0.5
BINS = [("TFIIIB", 30, 39, -16), ("ICR", 40, 49, 63)]
NC = len(BINS)

# ── 유닛 전체 β : 공개판과 똑같이 load_beta
UB = {(s, r, sl, a): v for (s, sr, r, sl, a), v in
      L.load_beta(min_cpg=10).items() if sr == SRC}

# ── 유전자부 β : calls 에서 끝 119 bp 만 · 가닥 두 줄은 합친다
raw = collections.defaultdict(lambda: [0, 0])
for node in ("w3", "iB", "iA"):
    p = os.path.join(L.R, "08_methylation_s5", node, L.METH[SRC],
                     "METH_READUNIT_%s.tsv.gz" % L.HOST[node])
    if not os.path.exists(p):
        continue
    with gzip.open(p, "rt", encoding="utf8") as f:
        ix = {k: i for i, k in enumerate(f.readline().rstrip(NL).split(TB))}
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[ix["data_source"]] != SRC or x[ix["in_final"]] != "Y":
                continue
            smp, sl = x[ix["sample"]], x[ix["slot"]]
            if smp in L.PB_EXCL or (smp, sl) in SHIFT:
                continue
            k = (smp, x[ix["read_id"]], sl, x[ix["array_pos"]])
            for t in (x[ix["calls"]] or "").split(","):
                if not t:
                    continue
                q, m = t.split(":")
                if int(q) >= GENE_FROM:
                    raw[k][0] += 1
                    raw[k][1] += (int(m) >= MLTHR)
GC = collections.Counter(v[0] for v in raw.values())
GB = {k: v[1] / float(v[0]) for k, v in raw.items() if v[0] >= MIN_GENE_CPG}
print("유닛 전체 β %s · 유전자부 β %s (CpG>=%d 관문에서 %s 개 탈락)"
      % (format(len(UB), ","), format(len(GB), ","), MIN_GENE_CPG,
         format(len(raw) - len(GB), ",")))

# ── FiberHMM 분절 : 공개판 규칙 (msp 만 뺀다)
SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio" or x[5] == "msp":
                continue
            SEG[(x[0], x[4], x[2], x[3])].append(
                ((int(x[6]) + int(x[7])) // 2, int(x[8])))

OBS = np.zeros((2, 2, NC, PER))
NUL = np.zeros((2, 2, NC, PER))
NMO = np.zeros((2, 2))
CF = np.zeros((2, 2))
rng = np.random.default_rng(20261012)   # 공개판(zzac)과 같은 씨앗 - 자를 맞춘다
for p in sorted(glob.glob(os.path.join(FD, "FHMM_READCOV_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio":
                continue
            if int(x[5]) > -HALF or int(x[6]) < HALF - 1:
                continue
            k = (x[0], x[4], x[2], x[3])
            labs = []
            if k in UB:
                labs.append((0, 0 if UB[k] < CUT_UNIT else 1))
            if k in GB:
                labs.append((1, 0 if GB[k] < CUT_GENE else 1))
            if not labs:
                continue
            if k in UB and k in GB:
                CF[0 if UB[k] < CUT_UNIT else 1, 0 if GB[k] < CUT_GENE else 1] += 1
            cen, cls = [], []
            for c, sz in SEG.get(k, []):
                if not (-HALF <= c < HALF):
                    continue
                for bi, (_n, a0, a1, _p) in enumerate(BINS):
                    if a0 <= sz <= a1:
                        cen.append(c + HALF)
                        cls.append(bi)
                        break
            cen = np.array(cen, dtype=int)
            cls = np.array(cls, dtype=int)
            d = rng.integers(0, PER, NB) if len(cen) else None
            for li, si in labs:
                NMO[li, si] += 1
                if not len(cen):
                    continue
                np.add.at(OBS[li, si], (cls, cen), 1)
                for t in range(NB):
                    np.add.at(NUL[li, si], (cls, (cen + d[t]) % PER), 1)
NUL /= float(NB)

K = np.ones(SM) / SM


def sm(v):
    return np.convolve(np.r_[v[-SM:], v, v[:SM]], K, "same")[SM:-SM]


R = np.zeros((2, NC, 2))
for li in (0, 1):
    for bi in range(NC):
        j = BINS[bi][3] + HALF
        for si in (0, 1):
            o, e = sm(OBS[li, si, bi]), sm(NUL[li, si, bi])
            R[li, bi, si] = o[j] / e[j] if e[j] > 0 else np.nan

np.savez_compressed(os.path.join(DD, "hprc_genebeta.npz"),
                    R=R, N=NMO, CF=CF, obs=OBS, nul=NUL,
                    gc=np.array(sorted(GC.items())),
                    sites=np.array([b[0] for b in BINS]),
                    pos=np.array([b[3] for b in BINS]))

LN = ["유닛 전체 β 0.508 (본표)", "유전자부 β 0.5 (세포주 방식)"]
print()
print("## 그 자리 귀무 대비 배수 — 자리 고정 (TFIIIB -16 · ICR +63 · 10 bp 칸)")
print("%-26s %-8s %9s %9s %9s %14s" % ("층 라벨", "자리", "저", "고", "저÷고", "분자 저/고"))
for li in (0, 1):
    for bi in range(NC):
        v = R[li, bi]
        print("%-26s %-8s %9.2f %9.2f %9.2f %14s"
              % (LN[li], BINS[bi][0], v[0], v[1],
                 v[0] / v[1] if v[1] > 0 else float("nan"),
                 "%s / %s" % (format(int(NMO[li, 0]), ","), format(int(NMO[li, 1]), ","))))
print()
print("## 두 라벨 혼동표 (둘 다 있는 유닛)")
print("                      유전자부 저    유전자부 고")
for i, nm in enumerate(("유닛 저", "유닛 고")):
    print("   %-16s %11s %13s" % (nm, format(int(CF[i, 0]), ","), format(int(CF[i, 1]), ",")))
tot = CF.sum()
print("   일치 %.1f%% (n %s)" % (100 * (CF[0, 0] + CF[1, 1]) / tot, format(int(tot), ",")))
print()
print("## 유전자부 CpG 수 분포")
for k in sorted(GC):
    if GC[k] > 50:
        print("   CpG %2d 개 : 유닛 %s" % (k, format(GC[k], ",")))
print()
print("저장 -> hprc_genebeta.npz")
