# -*- coding: utf-8 -*-
"""뉴클레오솜 없는 분자만으로 flank/core 를 다시 (2026-10-13 · 사용자 질문).

왜
  본표의 flank/core (MAIN f) 는 모든 분자를 풀링한다.
  그런데 한 분자에서 그 자리에 뉴클레오솜이 있으면 TF 는 못 붙는다 - 둘은 배타다.
  그러니 "ICR 비 2.08 이 뉴클레오솜 때문인가" 는 뉴클레오솜 없는 분자만 놓고
  다시 재면 바로 갈린다. (SUPP (h)(i) 의 호출기 비사용 판이다)

내는 것 : nucfree_flankcore_span.npz + 화면 표
  nA/nM  (자리 2 · 층 2 · 뉴클없음 2 · rel 401)
  nmol   (자리 2 · 층 2 · 뉴클없음 2)
  ntf    (자리 2 · 층 2 · 뉴클없음 2)   ← 그 무리 안의 점유 (천장 계산용)

정의 (zzab 과 같게 맞춘다)
  '뉴클레오솜 없음' = core 창에 겹치는 `nuc` 분절이 하나도 없음
  '점유'            = 중심이 구역 창 안 + 크기 범위인 `tf` 호출이 있음
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
SRC = "PACBIO_FIBERSEQ"
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
LO, HI = -200, 200
NR = HI - LO + 1

# 자리 : (중심, 구역창, TF 크기범위)  ← 본표와 같은 값
SITE = [("TFIIIB", -11, (-28, -10), (30, 59)),
        ("ICR", 58, (57, 68), (20, 49))]
NS = len(SITE)

BETA = {(s, r, sl, a): v for (s, sr, r, sl, a), v in
        L.load_beta(min_cpg=10).items() if sr == SRC}
print("β n=%s" % format(len(BETA), ","))

NUC = collections.defaultdict(list)
TFS = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FD, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip(NL).split(TB)
            if x[1] != "pacbio":
                continue
            k = (x[0], x[4], x[2], x[3])
            if x[5] == "nuc":
                NUC[k].append((int(x[6]), int(x[7])))
            elif x[5] == "tf":
                TFS[k].append((int(x[6]), int(x[7]), int(x[8])))
print("뉴클 가진 분자 %s · TF 가진 분자 %s"
      % (format(len(NUC), ","), format(len(TFS), ",")))

UNIT, _ = L.load_units()
TSSU = {k: (v[2] - v[0] + 1) for k, v in UNIT.items()}

nA = np.zeros((NS, 2, 2, NR))
nM = np.zeros((NS, 2, 2, NR))
nmol = np.zeros((NS, 2, 2))
ntf = np.zeros((NS, 2, 2))

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
    print("  %s : 분자 %s" % (node, format(len(MOL), ",")))
    for k, blocks in MOL.items():
        b = BETA.get(k)
        tu = TSSU.get((k[0], k[2], k[3]))
        if b is None or tu is None:
            continue
        si = 0 if b < CUT else 1
        nu = NUC.get(k, ())
        tf = TFS.get(k, ())
        # 자리마다 뉴클 유무 / TF 유무
        flag = []
        for _, c0, zw, szr in SITE:
            # 대칭 조건 : core 만 걸면 flank 가 안 걸려 비가 아래로 기울어진다.
            #   core 와 flank 둘 다 덮는 구간(c0-35 .. c0+27) 전체를 걷는다.
            a0, a1 = c0 - 35, c0 + 27
            cov = any(e >= a0 and s <= a1 for s, e in nu)
            hit = any(szr[0] <= sz <= szr[1] and zw[0] <= (s + e) // 2 <= zw[1]
                      for s, e, sz in tf)
            flag.append((0 if cov else 1, hit))     # 0=구간에 뉴클있음 · 1=구간 전체 깨끗
        for j, (fr, hit) in enumerate(flag):
            nmol[j, si, fr] += 1
            ntf[j, si, fr] += hit
        # m6A 쌓기
        for blk in blocks:
            for t in blk.split(","):
                if not t:
                    continue
                a_, b_ = t.split(":")
                r = int(a_) - tu
                if LO <= r <= HI:
                    m = int(b_) >= MLTHR
                    i = r - LO
                    for j, (fr, _h) in enumerate(flag):
                        nA[j, si, fr, i] += 1
                        nM[j, si, fr, i] += m

np.savez_compressed(os.path.join(DD, "nucfree_flankcore_span.npz"),
                    nA=nA, nM=nM, nmol=nmol, ntf=ntf, lo=LO, hi=HI, cut=CUT,
                    sites=np.array([s[0] for s in SITE]))

rel = np.arange(LO, HI + 1)


def ratio(j, si, fr):
    _, c0, _, _ = SITE[j]
    mc = (rel >= c0 - 7) & (rel <= c0 + 1)
    mf = ((rel >= c0 - 35) & (rel <= c0 - 27)) | ((rel >= c0 + 22) & (rel <= c0 + 27))
    cA, cM = nA[j, si, fr][mc].sum(), nM[j, si, fr][mc].sum()
    fA, fM = nA[j, si, fr][mf].sum(), nM[j, si, fr][mf].sum()
    cr = (cM + 0.5) / (cA + 0.5)
    fr_ = (fM + 0.5) / (fA + 0.5)
    return fr_ / cr, int(cA), int(fA)


print()
print("## flank/core — 뉴클레오솜 있는 분자 ↔ 없는 분자")
print("%-8s %-4s %-10s %9s %9s %9s %9s %9s"
      % ("자리", "층", "무리", "분자", "점유", "천장", "관측 비", "관측÷천장"))
for j, (nm, c0, _, _) in enumerate(SITE):
    for si, sl in ((0, "저"), (1, "고")):
        for fr, fn in ((0, "구간에 뉴클"), (1, "구간 깨끗")):
            n = nmol[j, si, fr]
            f = ntf[j, si, fr] / max(n, 1)
            ceil = 1.0 / max(1.0 - f, 1e-9)
            r, cA, fA = ratio(j, si, fr)
            print("%-8s %-4s %-10s %9s %8.4f %9.3f %9.3f %9.3f"
                  % (nm, sl, fn, format(int(n), ","), f, ceil, r, r / ceil))
print()
print("저장 -> nucfree_flankcore_span.npz")
