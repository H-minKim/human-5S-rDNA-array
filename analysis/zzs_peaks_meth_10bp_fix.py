# -*- coding: utf-8 -*-
"""[고침] zz7 과 같은 계산 · 귀무 배율 버그만 수정 (2026-10-11).

버그 : NSUB 는 footprint 가 있는 분자만 센다 (`if not cen: continue` 가 그 위에 있다).
        NMOL 은 전부 센다. 코드는 귀무를 NMOL/NSUB 로 늘린다.
        10 bp 8 칸(10~200 bp)에서는 거의 모든 분자가 footprint 를 하나는 가져
        배율이 3.28 / 3.28 / 3.44 (저/중/고) 로 SUB=3 과 거의 같았다 - 거의.
        층마다 달라서 저÷고 비가 약 5% 틀린다.
        (같은 버그가 좁은 칸에서는 5.2 / 7.2 / 8.6 배로 터졌다 - zzr 참고)
고침 : 귀무를 솎지 않는다 (SUB=1) -> 배율이 정확히 1 이므로 곱하지 않는다.
zz7 의 산출(fhmm_peaks.npz)은 그대로 둔다. 여기는 fhmm_peaks_v2.npz 로 낸다.
"""
import io
import os
import sys
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

DD = _need("RDNA_DATA_ROOT")
FD = os.path.join(DD, "fhmm")
NL = chr(10)
TB = chr(9)
PER = 2240                      # 유닛 주기
HALF = PER // 2                 # 1120
NB = 200                        # 귀무 반복
SUB = 1                         # 고침 : 솎지 않는다 (배율 = 1)
SM = 11                         # 자리 평활 (bp)
CLS = [(10, 19, "10~19"), (20, 29, "20~29"), (30, 39, "30~39"), (40, 49, "40~49"),
       (50, 59, "50~59"), (60, 89, "60~89"), (90, 140, "90~140"), (141, 200, "141~200")]
NC = len(CLS)
PLAT = {"pacbio": 1, "nanopore": 0}
SRCN = {"pacbio": "PACBIO_FIBERSEQ", "nanopore": "UCSC_ONT_FIBERSEQ"}

old = np.load(os.path.join(DD, "perbase_strat.npz"), allow_pickle=True)
CUTB = old["cuts"][:, 0, :]
BETA = {}
for plat, src in SRCN.items():
    for (s, sr, r, sl, a), v in L.load_beta(min_cpg=10).items():
        if sr == src:
            BETA[(plat, s, r, sl, a)] = v

SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[5] == "msp":
                continue
            SEG[(x[1], x[0], x[4], x[2], x[3])].append((int(x[6]), int(x[7]), int(x[8])))
print("구간 보유 read x 유닛 %s" % format(len(SEG), ","))

OBS = np.zeros((2, 3, NC, PER))         # 플랫폼 x 층 x 크기무리 x 자리
NUL = np.zeros((2, 3, NC, PER))
NMOL = np.zeros((2, 3))
NSUB = np.zeros((2, 3))
rng = np.random.default_rng(20261011)
nskip = 0

for p in sorted(glob.glob(os.path.join(FD, "FHMM_READCOV_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            plat, smp, sl, ap, rid = x[1], x[0], x[2], x[3], x[4]
            lo, hi = int(x[5]), int(x[6])
            if lo > -HALF or hi < HALF - 1:      # 창을 온전히 덮는 분자만
                nskip += 1
                continue
            b = BETA.get((plat, smp, rid, sl, ap))
            if b is None:
                continue
            pi = PLAT[plat]
            si = 0 if b < CUTB[pi, 0] else (1 if b < CUTB[pi, 1] else 2)
            fp = SEG.get((plat, smp, rid, sl, ap), [])
            cen, cls = [], []
            for (s0, s1, sz) in fp:
                c = (s0 + s1) // 2
                if not (-HALF <= c < HALF):
                    continue
                for k, (a0, a1, _) in enumerate(CLS):
                    if a0 <= sz <= a1:
                        cen.append(c + HALF)
                        cls.append(k)
                        break
            NMOL[pi, si] += 1
            if not cen:
                continue
            cen = np.array(cen)
            cls = np.array(cls)
            np.add.at(OBS[pi, si], (cls, cen), 1)
            if int(NMOL[pi, si]) % SUB:
                continue
            NSUB[pi, si] += 1
            d = rng.integers(0, PER, NB)
            for t in range(NB):
                np.add.at(NUL[pi, si], (cls, (cen + d[t]) % PER), 1)

NUL /= float(NB)
np.savez_compressed(os.path.join(DD, "fhmm_peaks_v2.npz"),
                    obs=OBS, nul=NUL, nmol=NMOL, nsub=NSUB,
                    per=PER, half=HALF, cls=np.array([c[2] for c in CLS]))
print("창을 다 안 덮어 버린 분자 %s" % format(nskip, ","))
print("분자 수  저 / 중 / 고")
for nm, pi in (("PacBio", 1), ("ONT", 0)):
    print("  %-8s %s   (귀무 표본 %s)"
          % (nm, " / ".join("%d" % NMOL[pi, s] for s in range(3)),
             " / ".join("%d" % NSUB[pi, s] for s in range(3))))

# ── 봉우리 표 : 관측 ÷ 귀무 (평활 후)
K = np.ones(SM) / SM


def sm(v):
    return np.convolve(np.r_[v[-SM:], v, v[:SM]], K, "same")[SM:-SM]


print()
print("=" * 86)
print("## 봉우리 — 관측 ÷ 자리 귀무 (평활 %d bp · 분자 규모로 맞춤)" % SM)
for nm, pi in (("PacBio", 1), ("ONT", 0)):
    print()
    print("  %s" % nm)
    print("    %-9s %10s %10s %10s %10s %10s"
          % ("크기", "저 자리", "저 배수", "중 배수", "고 배수", "저 관측/1000"))
    for k, (a0, a1, lab) in enumerate(CLS):
        o0 = sm(OBS[1 if nm == "PacBio" else 0, 0, k])
        n0 = sm(NUL[1 if nm == "PacBio" else 0, 0, k])     # 배율 1
        r0 = np.where(n0 > 0, o0 / np.maximum(n0, 1e-9), np.nan)
        j = int(np.nanargmax(np.where(np.arange(PER) < PER, r0, np.nan)))
        pos = j - HALF
        row = [r0[j]]
        for si in (1, 2):
            oo = sm(OBS[pi, si, k])
            nn = sm(NUL[pi, si, k])                            # 배율 1
            row.append(oo[j] / nn[j] if nn[j] > 0 else np.nan)
        print("    %-9s %+10d %10.2f %10.2f %10.2f %10.2f"
              % (lab, pos, row[0], row[1], row[2], 1000.0 * o0[j] / max(NMOL[pi, 0], 1)))
print()
print("저장 -> fhmm_peaks_v2.npz")
