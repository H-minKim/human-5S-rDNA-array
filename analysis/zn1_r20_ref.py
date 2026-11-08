# -*- coding: utf-8 -*-
"""r20 을 레퍼런스 문턱으로 다시 만든다 — F2 (a)(b) 가 이 표를 읽는다 (2026-10-11).

바꾸는 것 하나뿐 : `actuation` 의 정의
   옛 : MSP 가 덮은 염기 / 호출된 염기      ← 문턱 없음 (레퍼런스에 없는 지표)
   새 : large MSP(>150 bp) 가 덮은 염기 / 호출된 염기
        *"The cutoff for identifying di-nucleosomes and large MSPs was 210 and 150"*

나머지는 r20 과 똑같이 둔다 — 그래야 전/후를 같은 자리에서 맞댈 수 있다
   · 구역 : far −1500..−1000 · mid −900..−400 · up −300..−1 · gene 0..+118  (TSS 기준)
   · 자리별 윤곽 : 유닛 안 1..2241 에 유닛마다 자기 TSS 를 OFF=2123 에 맞춰 쌓는다
   · 검정 : 개체 단위 Wilcoxon 부호순위 (저메틸 ↔ 고메틸)
   · 제외 : PacBio fibertools 2 개체 · PacBio 자리이동 슬롯 2 개
배치는 본표 a23 과 같다(정본 좌표).
"""
import io
import os
import sys
import collections
import numpy as np
from scipy.stats import wilcoxon

# The data and output roots are environment variables with no default (CONVENTIONS, "the deposited copy wins").
#  저장소 안 파일로 가는 길은 `__file__` 에서 상대로 푼다 — `a_lib`·`bamlite`·`p02_lib` 가 이 폴더에 같이 예치돼 있다.
#  변수 표는 `README.md` §2 에 있다.
_HERE = os.path.dirname(os.path.abspath(__file__))


def _need(k):
    try:
        return os.environ[k]
    except KeyError:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s. "
            "The deposited copy has no default; see README.md section 2."
            % (k, os.path.basename(__file__)))

#  쓰는 폴더다 — 발행 입력 폴더를 가리키면 그 표를 덮는다.
W = _need("RDNA_RESULTMAP_DIR")
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

ULEN, OFF, LARGE = 2241, 2123, 150
REG = [("far", -1500, -1000), ("mid", -900, -400), ("up", -300, -1), ("gene", 0, 118)]
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
UNIT, _ = L.load_units()
GMM = L.load_state_gmm()

# (plat, smp, state, region) -> [large bp, 호출 bp]
RG = collections.defaultdict(lambda: np.zeros(2))
# (plat, state, unit_pos) -> [large bp, 호출 bp, 분자]
PR = collections.defaultdict(lambda: np.zeros(3))
nread = 0
for smp, src, rid, slot, arr, lo, us in L.iter_reads(UNIT, 0, 0):
    plat = "ONT" if src.startswith("UCSC") else "PacBio"   # 그림이 찾는 이름 (2026-10-12 고침)
    if src == "PACBIO_FIBERSEQ" and (smp in L.PB_EXCL or (smp, slot) in SHIFT):
        continue
    nread += 1
    cod = np.zeros(len(arr), dtype=np.int8)
    cod[arr == ord("N")] = 1
    cod[arr == ord("M")] = 2
    ch = np.flatnonzero(np.diff(cod)) + 1
    bb = np.concatenate(([0], ch, [len(cod)]))
    segs = [(bb[i], bb[i + 1], cod[bb[i]]) for i in range(len(bb) - 1)]
    big = np.zeros(len(arr), dtype=bool)
    for i in range(1, len(segs) - 1):
        a0, e0, c = segs[i]
        if c == 2 and segs[i - 1][2] != 0 and segs[i + 1][2] != 0 and (e0 - a0) > LARGE:
            big[a0:e0] = True
    for ap, _t in us:
        s_, e_, tss = UNIT[(smp, slot, str(ap))]
        st = GMM.get((smp, src, rid, slot, str(ap)), "NA")
        if st not in ("L", "H"):
            continue
        a0, b0 = s_ - 1 - lo, e_ - lo
        sc, sb = cod[a0:b0], big[a0:b0]
        if len(sc) == 0:
            continue
        tss_u = tss - s_ + 1
        up = np.arange(1, len(sc) + 1)
        rel = up - tss_u
        ok = sc > 0
        for nm, r0, r1 in REG:
            m = ok & (rel >= r0) & (rel <= r1)
            k = int(m.sum())
            if k:
                v = RG[(plat, smp, st, nm)]
                v[0] += int(sb[m].sum()); v[1] += k
        # 자리별 윤곽 : TSS 를 OFF 에 맞춘다
        tgt = OFF + rel
        inr = ok & (tgt >= 1) & (tgt <= ULEN)
        if inr.any():
            tt = tgt[inr]; bbv = sb[inr]
            for t_, b_ in zip(tt, bbv):
                p = PR[(plat, st, int(t_))]
                p[0] += int(b_); p[1] += 1
            for t_ in np.unique(tt):
                PR[(plat, st, int(t_))][2] += 1
print("read %d · 구역 칸 %d · 자리 칸 %d" % (nread, len(RG), len(PR)))

# ───────── unified (구역 · 개체 Wilcoxon)
OUT1 = os.path.join(W, "r20_actuation_unified_REF.tsv")
rows = []
for plat in ("ONT", "PacBio"):
    for nm, _a, _b in REG:
        aL, aH = [], []
        smps = {s for (p_, s, t_, r_) in RG if p_ == plat and r_ == nm}
        for s in sorted(smps):
            vL, vH = RG.get((plat, s, "L", nm)), RG.get((plat, s, "H", nm))
            if vL is None or vH is None or vL[1] < 500 or vH[1] < 500:
                continue
            aL.append(vL[0] / vL[1]); aH.append(vH[0] / vH[1])
        if len(aL) < 5:
            continue
        aL, aH = np.array(aL), np.array(aH)
        _s, p = wilcoxon(aL, aH)
        rows.append(("region", plat, nm, len(aL), aL.mean(), aH.mean(),
                     aL.mean() - aH.mean(), p))
with io.open(OUT1, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# r20 을 레퍼런스 문턱으로 : actuation = large MSP(>150bp) 덮은 염기 / 호출 염기"
            + chr(10))
    g.write("# 구역·검정·제외는 r20 과 같다 (개체 단위 Wilcoxon)" + chr(10))
    g.write(chr(9).join(["kind", "platform", "key", "n", "v_L", "v_H", "diff", "p"]) + chr(10))
    for r in rows:
        g.write(chr(9).join([r[0], r[1], r[2], str(r[3]), "%.5f" % r[4], "%.5f" % r[5],
                             "%.5f" % r[6], "%.3g" % r[7]]) + chr(10))

print()
print("=== 구역별 (large MSP >150 bp)")
print("   %-5s %-6s %4s %9s %9s %9s %10s" % ("플랫", "구역", "n", "저메틸", "고메틸", "차", "p"))
for r in rows:
    print("   %-5s %-6s %4d %9.5f %9.5f %9.5f %10.2e%s"
          % (r[1], r[2], r[3], r[4], r[5], r[6], r[7], "" if r[7] < 0.01 else "  n.s."))

# ───────── profile (자리별)
OUT2 = os.path.join(W, "r20_actuation_profile_REF.tsv")
with io.open(OUT2, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# 자리별 actuation = large MSP(>150bp) 덮은 염기 / 호출 염기 "
            "(유닛마다 자기 TSS 를 %d 에 맞춤)" % OFF + chr(10))
    g.write(chr(9).join(["platform", "state", "unit_pos", "n_mol", "actuation"]) + chr(10))
    n2 = 0
    for (plat, st, upx), v in sorted(PR.items()):
        if v[1] < 200:
            continue
        g.write(chr(9).join([plat, st, str(upx), str(int(v[2])),
                             "%.5f" % (v[0] / v[1])]) + chr(10))
        n2 += 1
print()
print("자리 줄 %d -> %s" % (n2, os.path.basename(OUT2)))
print("구역 줄 %d -> %s" % (len(rows), os.path.basename(OUT1)))
