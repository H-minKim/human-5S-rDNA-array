# -*- coding: utf-8 -*-
"""뉴클레오솜 제거 판 + 10 bp 양분 (2026-10-12).

둘을 한 패스로 낸다.
 ① fhmm_peaks_bin10.npz  : 10 bp 크기 무리 x 자리 · 양분(β 0.508) - Fig5Fa/Fb 용
 ② nucfree_bin.npz       : 자리마다 뉴클레오솜이 없는 분자만 추려 TF 점유를 다시 센다
      NCOV   그 자리를 덮는 분자 수
      NNONUC 그 자리에 뉴클레오솜이 없는 분자 수
      NTF    TF 콜이 그 자리를 덮는 분자 수
      NTFNN  뉴클레오솜이 없으면서 TF 콜이 덮는 분자 수
    -> 생판 점유 = NTF/NCOV · 뉴클 제거 점유 = NTFNN/NNONUC
 ③ 표S7 : 두 자리의 생판 ↔ 뉴클 제거 (zzb 정의 : 중심이 창 안 + 크기 범위)

왜 : 저메틸 분자가 더 '열려' 있으면 TF 콜 기회 자체가 많다. 뉴클레오솜이 없는 분자만
  놓고 다시 세면 그 기회 차이를 걷어낸 값이 된다.
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

NL = chr(10)
TB = chr(9)
DD = _need("RDNA_DATA_ROOT")
FD = os.path.join(DD, "fhmm")
OUTD = _need("RDNA_CAPTION_DIR")
CUT = 0.508
PER, HALF = 2240, 1120
NB = 200
LO, HI = -200, 160
NP = HI - LO + 1
CLS10 = [(10, 19, "10~19"), (20, 29, "20~29"), (30, 39, "30~39"), (40, 49, "40~49"),
         (50, 59, "50~59"), (60, 89, "60~89"), (90, 140, "90~140"), (141, 200, "141~200")]
SITE = {"TFIIIB": ((-28, -10), (30, 59)), "ICR": ((57, 68), (20, 49))}

SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio":
                continue
            SEG[(x[0], x[4], x[2], x[3])].append((x[5], int(x[6]), int(x[7]), int(x[8])))
BETA = {(s, r, sl, a): v for (s, sr, r, sl, a), v in
        L.load_beta(min_cpg=10).items() if sr == "PACBIO_FIBERSEQ"}

NC = len(CLS10)
OBS = np.zeros((2, NC, PER))
NUL = np.zeros((2, NC, PER))
NMOL = np.zeros(2)
NCOV = np.zeros((2, NP))
NNONUC = np.zeros((2, NP))
NTF = np.zeros((2, NP))
NTFNN = np.zeros((2, NP))
CNT = {nm: {si: [0, 0, 0, 0] for si in (0, 1)} for nm in SITE}
rng = np.random.default_rng(20261012)
idx = np.arange(LO, HI + 1)

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
            sg = SEG.get(k, [])
            NMOL[si] += 1
            # ① 10 bp 봉우리
            cen, cls = [], []
            for (c, s0, e0, sz) in sg:
                if c == "msp":
                    continue
                cc = (s0 + e0) // 2
                if not (-HALF <= cc < HALF):
                    continue
                for j, (a0, a1, _) in enumerate(CLS10):
                    if a0 <= sz <= a1:
                        cen.append(cc + HALF)
                        cls.append(j)
                        break
            if cen:
                cen = np.array(cen)
                cls = np.array(cls)
                np.add.at(OBS[si], (cls, cen), 1)
                d = rng.integers(0, PER, NB)
                for t in range(NB):
                    np.add.at(NUL[si], (cls, (cen + d[t]) % PER), 1)
            # ② 자리마다 뉴클레오솜 유무 x TF 유무
            nuc = np.zeros(NP, bool)
            tf = np.zeros(NP, bool)
            for (c, s0, e0, sz) in sg:
                a = max(LO, s0) - LO
                bq = min(HI, e0) - LO
                if bq < a:
                    continue
                if c == "nuc":
                    nuc[a:bq + 1] = True
                elif c == "tf" and 20 <= sz <= 164:
                    tf[a:bq + 1] = True
            NCOV[si] += 1
            NNONUC[si] += ~nuc
            NTF[si] += tf
            NTFNN[si] += (tf & ~nuc)
            # ③ 자리별 (zzb 정의)
            for nm, ((a0, a1), (z0, z1)) in SITE.items():
                hit = any(c == "tf" and z0 <= sz <= z1 and a0 <= (s0 + e0) // 2 <= a1
                          for c, s0, e0, sz in sg)
                nu = any(c == "nuc" and e0 >= a0 and s0 <= a1 for c, s0, e0, sz in sg)
                d2 = CNT[nm][si]
                d2[0] += 1
                d2[1] += hit
                if not nu:
                    d2[2] += 1
                    d2[3] += hit

NUL /= float(NB)
np.savez_compressed(os.path.join(DD, "fhmm_peaks_bin10.npz"),
                    obs=OBS, nul=NUL, nmol=NMOL, per=PER, half=HALF, cut=CUT,
                    cls=np.array([c[2] for c in CLS10]))
np.savez_compressed(os.path.join(DD, "nucfree_bin.npz"),
                    ncov=NCOV, nnonuc=NNONUC, ntf=NTF, ntfnn=NTFNN,
                    lo=LO, hi=HI, cut=CUT, nmol=NMOL)
print("분자 : 저 %s · 고 %s (컷 %.3f)"
      % (format(int(NMOL[0]), ","), format(int(NMOL[1]), ","), CUT))

print()
print("## 자리별 — 생판 ↔ 뉴클레오솜 제거 (zzb 정의)")
print("  %-8s %-4s %8s %9s %10s %12s" % ("site", "stratum", "n", "점유%", "n_nucleosome_free", "뉴클제거 점유%"))
ROWS = []
for nm in ("TFIIIB", "ICR"):
    r = []
    for si in (0, 1):
        a, t, ao, to = CNT[nm][si]
        print("  %-8s %-4s %8d %9.3f %10d %12.3f"
              % (nm, "low" if si == 0 else "high", a, 100.0 * t / a, ao,
                 100.0 * to / max(ao, 1)))
        r.append((a, t, ao, to))
    raw = (r[0][1] / r[0][0]) / max(r[1][1] / r[1][0], 1e-12)
    nf = (r[0][3] / max(r[0][2], 1)) / max(r[1][3] / max(r[1][2], 1), 1e-12)
    op = (r[0][2] / r[0][0]) / max(r[1][2] / r[1][0], 1e-12)
    print("    -> 저÷고 생판 %.2f · 뉴클 제거 %.2f · 열린비율 차 %.2f 배" % (raw, nf, op))
    ROWS.append((nm, r, raw, nf, op))

os.makedirs(OUTD, exist_ok=True)
out = os.path.join(OUTD, "TableS7_nucleosome_removed_occupancy.tsv")
with io.open(out, "w", encoding="utf8", newline=NL) as g:
    g.write("# 표S7 · 뉴클레오솜 가림을 걷어낸 점유 · 메틸 양분 (β %.3f)" % CUT + NL)
    g.write("# 자리 정의 : TF 콜의 중심이 창 안 + 크기 범위 안" + NL)
    g.write("#   TFIIIB 중심 -28..-10 · 크기 30~59   /   ICR 중심 +57..+68 · 크기 20~49" + NL)
    g.write("# 뉴클레오솜 없음 : 그 창에 겹치는 nuc 구간이 없는 분자" + NL)
    g.write("# 옛 기록 19.1 은 재현되지 않았다 - '뉴클레오솜 없음' 정의가 달랐던 것으로 보인다." + NL)
    g.write(TB.join(["site", "stratum", "n_molecules", "n_tf_calls", "occupancy", "n_nucleosome_free", "n_tf_calls_nucleosome_free",
                     "occupancy_nucleosome_free"]) + NL)
    for nm, r, raw, nf, op in ROWS:
        for si, lab in ((0, "low"), (1, "high")):
            a, t, ao, to = r[si]
            g.write(TB.join([nm, lab, "%d" % a, "%d" % t, "%.5f" % (t / a),
                             "%d" % ao, "%d" % to, "%.5f" % (to / max(ao, 1))]) + NL)
    g.write("#" + NL + TB.join(["site", "raw_low_over_high", "nucleosome_free_low_over_high", "open_fraction_difference"]) + NL)
    for nm, r, raw, nf, op in ROWS:
        g.write(TB.join([nm, "%.4f" % raw, "%.4f" % nf, "%.4f" % op]) + NL)
print()
print("적었다 ->", out)
print("저장 -> fhmm_peaks_bin10.npz · nucfree_bin.npz")
