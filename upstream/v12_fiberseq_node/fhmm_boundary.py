# -*- coding: utf-8 -*-
"""레퍼런스 창을 무엇이 채우는가 — 경계 분포로 직접 (2026-10-11).

레퍼런스 창 = 크기 20~180 bp footprint 의 경계(시작 또는 끝)가 tx -54..-46.
(Data S4 cell 47 : find_fp_starts(track, 20, 180, -55+250, -45+250) · 부등호 strict)

내는 것
  ① 자리별 경계 밀도 (시작 / 끝 따로) — 크기 무리별
  ② 레퍼런스 창을 채운 footprint 의 크기·중심 분포
  ③ 창을 만족시킨 것이 시작인가 끝인가
"""
import io
import os
import collections
import numpy as np

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
W = H("~/fhmm")
TLO, THI = -250, 249
NT = THI - TLO + 1
NOCOV = 2
R_SMIN, R_SMAX = 20, 180
WLO, WHI = -54, -46                     # diff 인덱스 196..204 에 해당
CLS = [(20, 60, "20~60"), (61, 90, "61~90"), (91, 140, "91~140"), (141, 180, "141~180")]

SEG = collections.defaultdict(list)
f = io.open(os.path.join(W, "FHMM_UNITS.tsv"), encoding="utf8")
f.readline()
for L in f:
    x = L.rstrip(NL).split(TB)
    if x[4] in ("tf", "nuc"):
        SEG[(x[0], x[3], x[1], x[2])].append((int(x[5]), int(x[6]), int(x[7]), x[4]))

COV = {}
f = io.open(os.path.join(W, "FHMM_READCOV.tsv"), encoding="utf8")
f.readline()
for L in f:
    x = L.rstrip(NL).split(TB)
    COV[(x[0], x[3], x[1], x[2])] = (int(x[4]), int(x[5]))

NCO = np.zeros(NT)
ST = {c[2]: np.zeros(NT) for c in CLS}
EN = {c[2]: np.zeros(NT) for c in CLS}
hit_sz, hit_cen, hit_side, hit_cls = [], [], collections.Counter(), collections.Counter()
nref = ntot = 0

for k, (lo, hi) in COV.items():
    a0, a1 = max(TLO, lo), min(THI, hi)
    if a1 < a0:
        continue
    ntot += 1
    NCO[a0 - TLO:a1 - TLO + 1] += 1
    ok = False
    for (s0, s1, sz, cls) in SEG.get(k, []):
        if not (R_SMIN <= sz <= R_SMAX):
            continue
        lab = None
        for (c0, c1, nm) in CLS:
            if c0 <= sz <= c1:
                lab = nm
                break
        if lab:
            if TLO <= s0 <= THI:
                ST[lab][s0 - TLO] += 1
            if TLO <= s1 <= THI:
                EN[lab][s1 - TLO] += 1
        bs = WLO <= s0 <= WHI
        be = WLO <= s1 <= WHI
        if bs or be:
            ok = True
            hit_sz.append(sz)
            hit_cen.append((s0 + s1) // 2)
            hit_side["시작" if bs else ""] += 1 if bs else 0
            hit_side["끝"] += 1 if be else 0
            hit_cls[lab or "기타"] += 1
    if ok:
        nref += 1

print("read x 유닛 %d · 레퍼런스 창 양성 %d (%.2f%%)" % (ntot, nref, 100.0 * nref / max(ntot, 1)))
hs = np.array(hit_sz)
hc = np.array(hit_cen)
if len(hs):
    print()
    print("== 창을 채운 footprint 가 무엇인가 (n=%d) ==" % len(hs))
    print("  크기  중앙 %.0f  (25~75%% %.0f~%.0f)  최소 %d  최대 %d"
          % (np.median(hs), *np.percentile(hs, [25, 75]), hs.min(), hs.max()))
    print("  중심  중앙 %+.0f (25~75%% %+.0f~%+.0f)" % (np.median(hc), *np.percentile(hc, [25, 75])))
    print("  경계  시작 %d · 끝 %d" % (hit_side["시작"], hit_side["끝"]))
    print("  크기무리 " + " · ".join("%s=%d" % (k, v) for k, v in sorted(hit_cls.items())))

print()
print("== 자리별 경계 밀도 (20 bp 칸 · 섬유당) — 시작 / 끝 ==")
print("  %-14s" % "자리" + "".join("%18s" % c[2] for c in CLS))
for a in range(-160, 160, 20):
    i0, i1 = a - TLO, a + 20 - TLO
    c = NCO[i0:i1].mean()
    row = "  %+6d..%+5d" % (a, a + 19)
    for (_, _, nm) in CLS:
        row += "%9.4f/%-8.4f" % (ST[nm][i0:i1].sum() / max(c, 1),
                                 EN[nm][i0:i1].sum() / max(c, 1))
    mk = "   <= 레퍼런스 창" if a == -60 else ""
    print(row + mk)
