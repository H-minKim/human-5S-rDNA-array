# -*- coding: utf-8 -*-
"""자리별 집계 한 번 + 제자리 음성대조(sliding) + 두 플랫폼 합쳐 판정 (2026-10-09).

`zf6` 이 틀린 곳 (사용자 지적으로 적발)
  ① 4 조합을 독립 반복처럼 다뤘다. 그런데 ONT·유닛β 와 ONT·유전자부β 는
     같은 read(91,819 ↔ 91,774) 다 ⇒ 독립은 ONT · PacBio 둘뿐이다.
  ② ○× 로만 봤다(5% 밖인가). 두 플랫폼이 모두 하위 15% 라는 사실을 버렸다.
  ③ 대조를 NTS(유전자에서 400~2,000bp 떨어진 곳) 에 깔았다.
     그런데 A-box·IE·C-box 는 유전자 안에 있다. 유전자는 전사되는 곳이라 NTS 와 애초에 다르다
     ⇒ 요소가 특별한지 보려면 대조도 유전자 안이어야 한다.

이 글이 하는 것
  · 한 패스에 자리별(rel) nA·nM 을 층마다 쌓는다 (이 배열 하나로 뒤의 모든 계산을 한다)
  · 대조를 세 구역에서 각각 sliding 으로 만든다
      유전자안(요소 뺀) / 상류90 / 먼 NTS
  · 두 플랫폼을 합쳐 판정한다 (한쪽 백분위 곱이 아니라 양쪽 모두 같은 꼬리인지)
"""
import io
import os
import sys
import gzip
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

OUTD = _need("RDNA_OUT_DIR")
DD = _need("RDNA_DATA_ROOT")
NPZ = os.path.join(DD, "perbase_strat.npz")
NL = chr(10)
MLTHR = 128
LO, HI = -2100, 170
BGW = (-2000, -100)
MINCPG_GENE = 5
SHIFT = {("HG02559", "2"), ("HG03672", "1")}
ST = ["저", "중", "고"]
AX = ["유닛 β", "유전자부 β"]
NW = HI - LO + 1
SRCS = [("UCSC_ONT_FIBERSEQ", "ONT"), ("PACBIO_FIBERSEQ", "PacBio")]
ELEM = [("TFIIIB", -30, -10), ("A-box", 49, 63), ("IE", 66, 71),
        ("C-box", 79, 96), ("ICR", 49, 96), ("유전자", 0, 118)]
# 대조 구역 : (이름, 시작, 끝, 요소 겹침 허용 안 하는 창)
CTLZ = [("유전자안", 0, 118), ("상류90", -90, -1), ("먼 NTS", -2000, -400)]

if os.path.exists(NPZ) and os.environ.get("REBUILD", "") != "1":
    z = np.load(NPZ, allow_pickle=True)
    NAr, NMr, CUT, NU = z["nA"], z["nM"], z["cuts"], z["nunit"]
    print("자리별 집계를 다시 쓴다 ->", NPZ)
else:
    UNIT, _ = L.load_units()
    TSSU = {k: (v[2] - v[0] + 1) for k, v in UNIT.items()}
    NAr = np.zeros((2, 2, 3, NW))          # 플랫폼 x 층화축 x 층 x 자리
    NMr = np.zeros((2, 2, 3, NW))
    NU = np.zeros((2, 2, 3))
    CUT = np.zeros((2, 2, 2))
    for pi, (SRC, nm_src) in enumerate(SRCS):
        UB = {(s, r, sl, a): v for (s, sr, r, sl, a), v in L.load_beta(min_cpg=10).items() if sr == SRC}
        GB = {}
        for node in ("w3", "iB", "iA"):
            p = os.path.join(L.R, "08_methylation_s5", node, L.METH[SRC],
                             "METH_READUNIT_%s.tsv.gz" % L.HOST[node])
            if not os.path.exists(p):
                continue
            with gzip.open(p, "rt", encoding="utf8") as f:
                f.readline()
                for l in f:
                    x = l.rstrip(NL).split(chr(9))
                    if len(x) < 13 or x[1] != SRC or x[11] != "Y" or not x[12]:
                        continue
                    tu = TSSU.get((x[0], x[4], x[5]))
                    if tu is None:
                        continue
                    v = np.array(x[12].replace(":", ",").split(","), dtype=np.int64).reshape(-1, 2)
                    rel = v[:, 0] - tu
                    m = (rel >= 0) & (rel <= 118)
                    n = int(m.sum())
                    if n >= MINCPG_GENE:
                        GB[(x[0], x[6], x[4], x[5])] = (v[m, 1] >= MLTHR).sum() / float(n)
        for ax, d in enumerate((UB, GB)):
            CUT[pi, ax] = np.percentile(np.array(list(d.values())), [33.333, 66.667])
        print("%-7s 유닛β n %d 절단 %.4f/%.4f · 유전자부β n %d 절단 %.4f/%.4f"
              % (nm_src, len(UB), CUT[pi, 0, 0], CUT[pi, 0, 1], len(GB), CUT[pi, 1, 0], CUT[pi, 1, 1]))
        used = 0
        for node in ("w3", "iB", "iA"):
            p = os.path.join(L.R, "08_methylation_s5", node, L.METH[SRC],
                             "M6A_READUNIT_%s.tsv.gz" % L.HOST[node])
            if not os.path.exists(p):
                continue
            with gzip.open(p, "rt", encoding="utf8") as f:
                f.readline()
                for l in f:
                    x = l.rstrip(NL).split(chr(9))
                    if len(x) < 14 or x[1] != SRC or x[12] != "Y" or not x[13]:
                        continue
                    smp, slot, ap, rid = x[0], x[4], x[5], x[6]
                    if SRC == "PACBIO_FIBERSEQ" and (smp in L.PB_EXCL or (smp, slot) in SHIFT):
                        continue
                    tu = TSSU.get((smp, slot, ap))
                    if tu is None:
                        continue
                    bs = (UB.get((smp, rid, slot, ap)), GB.get((smp, rid, slot, ap)))
                    if bs[0] is None and bs[1] is None:
                        continue
                    v = np.array(x[13].replace(":", ",").split(","), dtype=np.int64).reshape(-1, 2)
                    rel = v[:, 0] - tu
                    m = (rel >= LO) & (rel <= HI)
                    if not m.any():
                        continue
                    idx = rel[m] - LO
                    hit = (v[m, 1] >= MLTHR)
                    for ax in (0, 1):
                        b = bs[ax]
                        if b is None:
                            continue
                        si = 0 if b < CUT[pi, ax, 0] else (1 if b < CUT[pi, ax, 1] else 2)
                        np.add.at(NAr[pi, ax, si], idx, 1)
                        if hit.any():
                            np.add.at(NMr[pi, ax, si], idx[hit], 1)
                        NU[pi, ax, si] += 1
                    used += 1
            print("   읽음 %s · 누적 %d" % (node, used))
    np.savez_compressed(NPZ, nA=NAr, nM=NMr, cuts=CUT, nunit=NU, lo=LO, hi=HI)
    print("저장 ->", NPZ)


def inter(pi, ax, a, b):
    """창 [a,b] 의 교호작용 = [rate(고)/rate(저)] / [바탕 rate(고)/rate(저)] · A 관측 최소값"""
    i0, i1 = a - LO, b - LO + 1
    A = NAr[pi, ax, :, i0:i1].sum(1)
    M = NMr[pi, ax, :, i0:i1].sum(1)
    j0, j1 = BGW[0] - LO, BGW[1] - LO + 1
    bA = NAr[pi, ax, :, j0:j1].sum(1)
    bM = NMr[pi, ax, :, j0:j1].sum(1)
    if A[0] < 300 or A[2] < 300:
        return float("nan"), min(A[0], A[2])
    r = (M / A)
    br = (bM / bA)
    return (r[2] / r[0]) / (br[2] / br[0]), min(A[0], A[2])


def nulls(pi, ax, size, zname, z0, z1, excl):
    """구역 [z0,z1] 을 size 창으로 쓸어가며 교호작용 분포 · excl 과 겹치는 창은 뺀다"""
    out = []
    for s in range(z0, z1 - size + 2):
        e = s + size - 1
        if any(not (e < ea or s > eb) for ea, eb in excl):
            continue
        v, n = inter(pi, ax, s, e)
        if v == v:
            out.append(v)
    return np.array(out)


EXCL = [(-30, -10), (49, 96)]        # 요소가 든 자리는 대조에서 뺀다
print()
print("=== 요소 교호작용과 제자리 음성대조 (sliding · 두 플랫폼 · 층화축 2)")
rows = []
for ax in (0, 1):
    print()
    print("---- 층화축 : %s" % AX[ax])
    print("  %-8s %5s | %-26s | %-26s" % ("요소", "크기", "ONT  값 / 대조중앙 / 백분위", "PacBio  값 / 대조중앙 / 백분위"))
    for nm, a, b in ELEM:
        size = b - a + 1
        zname, z0, z1 = ("상류90", -90, -1) if nm == "TFIIIB" else ("유전자안", 0, 118)
        if size > (z1 - z0 + 1) - 20:
            zname, z0, z1 = "먼 NTS", -2000, -400
        cells = []
        for pi, (_, src) in enumerate(SRCS):
            v, n = inter(pi, ax, a, b)
            cn = nulls(pi, ax, size, zname, z0, z1, EXCL)
            if v != v or len(cn) < 8:
                cells.append("  값없음 (대조 %d)" % len(cn)); continue
            pc = 100.0 * (cn < v).mean()
            cells.append("%6.3f / %6.3f / %5.1f%% (n%d)" % (v, np.median(cn), pc, len(cn)))
            rows.append((AX[ax], src, nm, size, zname, v, float(np.median(cn)),
                         float(np.percentile(cn, 5)), float(np.percentile(cn, 95)), pc, len(cn), int(n)))
        print("  %-8s %5d | %-26s | %-26s" % (nm, size, cells[0], cells[1]))
        print("           %5s   대조 구역 = %s" % ("", zname))

print()
print("=== 두 플랫폼 합쳐 판정 (독립 반복은 ONT · PacBio 둘)")
print("  %-8s %-10s %10s %10s %10s" % ("요소", "층화축", "ONT 백분위", "PB 백분위", "판정"))
for ax in (0, 1):
    for nm, a, b in ELEM:
        pcs = [r[9] for r in rows if r[0] == AX[ax] and r[2] == nm]
        if len(pcs) < 2:
            continue
        lo_both = all(p <= 15 for p in pcs)
        hi_both = all(p >= 85 for p in pcs)
        vd = "양쪽 하위 15%" if lo_both else ("양쪽 상위 15%" if hi_both else "어긋남/안")
        print("  %-8s %-10s %9.1f%% %9.1f%% %10s" % (nm, AX[ax], pcs[0], pcs[1], vd))

out = os.path.join(OUTD, "TF_20261010_in_place_negative_control_sliding.tsv")
with io.open(out, "w", encoding="utf8", newline=NL) as g:
    g.write(chr(9).join(["층화축", "플랫폼", "요소", "창크기", "대조구역", "교호작용",
                         "대조중앙", "대조5%", "대조95%", "백분위", "대조창수", "A관측최소"]) + NL)
    for r in rows:
        g.write(chr(9).join([r[0], r[1], r[2], "%d" % r[3], r[4]]
                            + ["%.4f" % r[i] for i in range(5, 10)]
                            + ["%d" % r[10], "%d" % r[11]]) + NL)
print()
print("적었다 ->", out)
