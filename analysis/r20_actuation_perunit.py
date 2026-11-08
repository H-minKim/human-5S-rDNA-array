# -*- coding: utf-8 -*-
"""R20 — R18 과 같되 구역을 유닛마다 자기 TSS 에 맞춰 자른다.  [for_paper · 2026-10-06]

왜 : 트랙은 유닛 안 위치로 색인되는데 r18 은 구역을 고정 OFF=2123 에서 잘랐다.
      실측 TSS 는 2113~2126 으로 흩어져(z03) gene 구역이 ~2% 과소추정됐다(z04).
고친 곳 : 트랙을 쌓는 자리 한 군데 — 자리 i 를 OFF+(i−tss_u) 에 쌓는다.
           ⇒ OFF 가 그 유닛의 TSS 를 가리키므로 아래 계산은 r18 과 같은 코드다.
r18 은 그대로 둔다(옛 판 보존). 출력만 r20_* 로 나간다.


레퍼런스(FIRE/fibertools) 정량 단위
  percent actuation = 자리마다 : (그 자리에 MSP 를 가진 분자) / (그 자리를 덮는 분자)
  ※ 우리가 쓰던 nucfrac(=N bp/(N+M) bp)은 bp 비율이라 레퍼런스 단위가 아니다.

주 플랫폼 = ONT · PacBio(HiFi)는 병기만. PacBio 전용 검정(mono/kb · nucA)은 본 결과에서 뺀다.

한 번에 내는 것 (전부 같은 단위로)
  ① 구역별 actuation (개체 단위 · Wilcoxon)
  ② 창 의존성 : 구역 창 ↔ 유닛 전체
  ③ 층 의존성 : 자리당 분자 수 층
  ④ 품질 교란 : 분자를 QUAL 1 점 칸으로 H 분포에 맞춰 재가중
  ⑤ 미처리 배경 : 미처리 대조에서 같은 계산 (원값 · 품질 맞춤 후)
  ⑥ 자리별 윤곽 (그림용)
출력 r20_actuation_unified.tsv · r20_actuation_profile.tsv
"""
import io, os, sys, gzip, glob, struct, bisect, collections
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

sys.path.insert(0, _HERE)
from p02_lib import seq2ref

#  쓰는 폴더다 — 발행 입력 폴더를 가리키면 그 표를 덮는다.
W = _need("RDNA_RESULTMAP_DIR")
R = _need("RDNA_V12_CANON_DIR")
S6 = _need("RDNA_S6_DIR")
#  파일 이름 라벨 (RFC1918 사설 대역 · 주소가 아니다) · `NODE_LABELS` 로 덮어쓸 수 있다
HOST = {"w3": "ip-172-31-58-58", "iB": "ip-172-31-32-12", "iA": "ip-172-31-57-143"}
if os.environ.get("NODE_LABELS"):
    HOST = dict(kv.split("=", 1) for kv in os.environ["NODE_LABELS"].split(",") if "=" in kv)
SUB = {"UCSC_ONT_FIBERSEQ": "ont_meth", "PACBIO_FIBERSEQ": "pb_meth"}
EXCL = {"PACBIO_FIBERSEQ": {"HG02630", "hg002", "HG02624", "HG03672"}, "UCSC_ONT_FIBERSEQ": set()}
OFF, ULEN = 2123, 2241
REG = [("far", -1500, -1000), ("mid", -900, -400), ("up", -300, -1), ("gene", 0, 118)]
MINCOV = 20
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
_CIG = "MIDNSHP=X"


def read_bam_q(path):
    with gzip.open(path, "rb") as f:
        data = f.read()
    assert data[:4] == b"BAM\x01"
    i = 4
    lt = struct.unpack_from("<i", data, i)[0]; i += 4 + lt
    nr = struct.unpack_from("<i", data, i)[0]; i += 4
    refs = []
    for _ in range(nr):
        ln = struct.unpack_from("<i", data, i)[0]; i += 4
        refs.append(data[i:i + ln - 1].decode()); i += ln + 4
    while i < len(data):
        bs = struct.unpack_from("<i", data, i)[0]
        rec = data[i + 4:i + 4 + bs]; i += 4 + bs
        refid, pos, lrn, mapq, _b, ncig, flag, lseq = struct.unpack_from("<iiBBHHHi", rec, 0)
        o = 32
        name = rec[o:o + lrn - 1].decode(); o += lrn
        cig = []
        for k in range(ncig):
            v = struct.unpack_from("<I", rec, o + 4 * k)[0]
            cig.append((_CIG[v & 15], v >> 4))
        o += 4 * ncig + (lseq + 1) // 2
        qual = np.frombuffer(rec[o:o + lseq], dtype=np.uint8).astype(np.int16)
        yield dict(name=name, flag=flag, ref=refs[refid] if refid >= 0 else None,
                   pos=pos, cigar=cig, qual=qual)


UNITS = collections.defaultdict(list)
with io.open(os.path.join(R, "02_canonical", "MAP_V11.tsv"), encoding="utf8") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for l in f:
        x = l.rstrip("\n").split("\t")
        if x[ix["seg_type"]] == "UNIT":
            UNITS[(x[ix["sample"]], x[ix["slot"]])].append(
                (int(x[ix["unit_start"]]), int(x[ix["unit_end"]]), x[ix["array_pos"]]))

ST = {}
with gzip.open(os.path.join(S6, "a13_unit_state_gmm.tsv.gz"), "rt", encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip("\n").split("\t")
        ST[(x[1], x[0], x[2], x[3], x[4])] = x[10]

# ---------- 1 차 : 유닛별 평균 QUAL
import sys as _sys
_sys.path.insert(0, _HERE)
import a_lib as _alib
print("0 차 : 유닛별 TSS 적재", flush=True)
_UN, _ = _alib.load_units()
TSSU = {k: (v[2] - v[0] + 1) for k, v in _UN.items()}
del _UN
print("  유닛 %d · TSS 중앙 %d" % (len(TSSU), int(np.median(list(TSSU.values())))), flush=True)

print("1 차 : BAM 에서 유닛별 QUAL", flush=True)
QU = {}
BAMS = []
for src, sub in SUB.items():
    for node in HOST:
        BAMS += [(src, q) for q in sorted(glob.glob(os.path.join(R, "08_methylation_s6", node, "s6", sub, "bam", "*.G8v3.bam")))]
for src, bam in BAMS:
    smp = os.path.basename(bam).split(".")[0]
    if smp in EXCL.get(src, set()):
        continue
    for r in read_bam_q(bam):
        if r["flag"] & 0x900 or r["ref"] is None:
            continue
        try:
            sm, slot = r["ref"].rsplit("_h", 1)
        except ValueError:
            continue
        us = sorted(UNITS.get((sm, slot), []))
        if not us:
            continue
        ustart = [u[0] for u in us]
        acc = collections.defaultdict(lambda: [0, 0])
        for sp, g in seq2ref(r).items():
            j = bisect.bisect_right(ustart, g) - 1
            if j < 0:
                continue
            a, b, ap = us[j]
            if a <= g <= b and sp < len(r["qual"]):
                v = acc[ap]; v[0] += int(r["qual"][sp]); v[1] += 1
        for ap, v in acc.items():
            if v[1] >= 500:
                QU[(smp, src, r["name"], slot, ap)] = v[0] / v[1]
print("  QUAL 유닛 %d" % len(QU), flush=True)

# ---------- 2 차 : 트랙을 흘리며 자리별 누적
NB = ULEN + 2
P_SMP = collections.defaultdict(lambda: np.zeros((NB, 2)))      # (src,smp,state)
P_Q = collections.defaultdict(lambda: np.zeros((NB, 2)))        # (src,state,qbin)
QCNT = collections.Counter()                                     # (src,state,qbin) -> 분자 수
NOTSS = collections.Counter() # units with no TSS found
print("2 차 : 트랙 누적", flush=True)
for p in sorted(glob.glob(os.path.join(S6, "a01_tracks", "*.tsv.gz"))):
    with gzip.open(p, "rt", encoding="utf8") as f:
        ix = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
        for l in f:
            x = l.rstrip("\n").split("\t")
            smp, src = x[ix["sample"]], x[ix["src"]]
            if smp in EXCL.get(src, set()):
                continue
            key = (smp, src, x[ix["read_id"]], x[ix["slot"]], x[ix["array_pos"]])
            st = ST.get(key, "U")
            if st not in ("L", "H"):
                continue
            tss_u = TSSU.get((smp, x[ix["slot"]], x[ix["array_pos"]]))
            if tss_u is None:
                NOTSS[src] += 1
                continue
            a = np.frombuffer(x[ix["track"]].encode(), dtype=np.uint8)
            cov = (a != ord(".")).astype(np.float64)
            msp = (a == ord("M")).astype(np.float64)
            # 자리 i(1기반) -> OFF + (i - tss_u) · 배열 밖은 잘라낸다
            sh = OFF - tss_u
            j0 = 1 + sh                      # 자리 1 이 놓이는 칸
            k0 = max(0, 1 - j0)   # 하한 1 — 칸 0 에 쓰지 않는다(r18 과 같은 경계)                 # 트랙에서 버릴 앞쪽
            k1 = min(len(a), ULEN + 1 - j0)   # 상한 ULEN(2241) — r18 과 같은 경계        # 트랙에서 버릴 뒤쪽
            if k1 <= k0:
                continue
            d0, d1 = j0 + k0, j0 + k1
            A = P_SMP[(src, smp, st)]
            A[d0:d1, 0] += cov[k0:k1]; A[d0:d1, 1] += msp[k0:k1]
            q = QU.get(key)
            if q is not None:
                qb = int(round(q))
                B = P_Q[(src, st, qb)]
                B[d0:d1, 0] += cov[k0:k1]; B[d0:d1, 1] += msp[k0:k1]
                QCNT[(src, st, qb)] += 1
print("  개체×상태 %d · QUAL 칸 %d · TSS 없음 %s" % (len(P_SMP), len(QCNT), dict(NOTSS)), flush=True)


def act(A, a0, b0, mincov=MINCOV):
    m = A[a0:b0, 0] >= mincov
    if m.sum() < 30:
        return float("nan")
    return float((A[a0:b0, 1][m] / A[a0:b0, 0][m]).mean())


rows = []
SRCS = (("UCSC_ONT_FIBERSEQ", "ONT"), ("PACBIO_FIBERSEQ", "PacBio"))

print("\n" + "=" * 100)
print("① 구역별 percent actuation (개체 단위)")
print("%-7s %-5s %6s %9s %9s %11s %11s" % ("src", "구역", "개체", "act_L", "act_H", "L-H", "p"))
for src, s6 in SRCS:
    smps = sorted({k[1] for k in P_SMP if k[0] == src})
    for rn, ra, rb in REG:
        a0, b0 = OFF + ra, OFF + rb + 1
        dL, dH = [], []
        for smp in smps:
            AL = P_SMP.get((src, smp, "L")); AH = P_SMP.get((src, smp, "H"))
            if AL is None or AH is None:
                continue
            vl, vh = act(AL, a0, b0), act(AH, a0, b0)
            if np.isfinite(vl) and np.isfinite(vh):
                dL.append(vl); dH.append(vh)
        if len(dL) < 8:
            continue
        dL, dH = np.array(dL), np.array(dH)
        _, pv = wilcoxon(dL, dH)
        print("%-7s %-5s %6d %9.4f %9.4f %+11.4f %11.2e" %
              (s6, rn, len(dL), dL.mean(), dH.mean(), (dL - dH).mean(), pv))
        rows.append(("region", s6, rn, len(dL), round(float(dL.mean()), 5),
                     round(float(dH.mean()), 5), round(float((dL - dH).mean()), 5), "%.3g" % pv))

print("\n" + "=" * 100)
print("② 창 의존성 — 구역 창 ↔ 유닛 전체")
print("%-7s %-12s %6s %11s %11s" % ("src", "창", "개체", "L-H", "p"))
for src, s6 in SRCS:
    smps = sorted({k[1] for k in P_SMP if k[0] == src})
    for wn, a0, b0 in (("구역 창", OFF - 1500, OFF + 119), ("유닛 전체", 1, ULEN + 1)):
        dL, dH = [], []
        for smp in smps:
            AL = P_SMP.get((src, smp, "L")); AH = P_SMP.get((src, smp, "H"))
            if AL is None or AH is None:
                continue
            vl, vh = act(AL, a0, b0), act(AH, a0, b0)
            if np.isfinite(vl) and np.isfinite(vh):
                dL.append(vl); dH.append(vh)
        if len(dL) < 8:
            continue
        dL, dH = np.array(dL), np.array(dH)
        _, pv = wilcoxon(dL, dH)
        print("%-7s %-12s %6d %+11.4f %11.2e" % (s6, wn, len(dL), (dL - dH).mean(), pv))
        rows.append(("window", s6, wn, len(dL), round(float(dL.mean()), 5),
                     round(float(dH.mean()), 5), round(float((dL - dH).mean()), 5), "%.3g" % pv))

print("\n" + "=" * 100)
print("③ 층 의존성 — 자리당 분자 수 층 (유닛 전체 창)")
print("%-7s %-12s %10s %10s %11s" % ("src", "분자 수 층", "act_L", "act_H", "L-H"))
CB = [(20, 50), (50, 100), (100, 200), (200, 500), (500, 10 ** 9)]
for src, s6 in SRCS:
    TL = np.zeros((NB, 2)); TH = np.zeros((NB, 2))
    for k, A in P_SMP.items():
        if k[0] != src:
            continue
        (TL if k[2] == "L" else TH)[:] += A
    for lo, hi in CB:
        m = (TL[:, 0] >= lo) & (TL[:, 0] < hi) & (TH[:, 0] >= lo) & (TH[:, 0] < hi)
        if m.sum() < 30:
            continue
        vl = float((TL[m, 1] / TL[m, 0]).mean()); vh = float((TH[m, 1] / TH[m, 0]).mean())
        print("%-7s %-12s %10.4f %10.4f %+11.4f   (자리 %d)" %
              (s6, "%d-%s" % (lo, "" if hi > 10 ** 8 else str(hi)), vl, vh, vl - vh, int(m.sum())))
        rows.append(("by_cov", s6, "%d-%s" % (lo, "" if hi > 10 ** 8 else str(hi)), int(m.sum()),
                     round(vl, 5), round(vh, 5), round(vl - vh, 5), ""))

print("\n" + "=" * 100)
print("④ 품질 교란 — QUAL 1 점 칸으로 H 분포에 맞춰 L 재가중 (유닛 전체 창)")
print("%-7s %11s %11s %9s" % ("src", "맞춤 전", "맞춤 후", "유지율"))
for src, s6 in SRCS:
    qbs = sorted({k[2] for k in QCNT if k[0] == src})
    TL = np.zeros((NB, 2)); TH = np.zeros((NB, 2))
    WL = np.zeros((NB, 2)); WH = np.zeros((NB, 2))
    for qb in qbs:
        nl = QCNT.get((src, "L", qb), 0); nh = QCNT.get((src, "H", qb), 0)
        AL = P_Q.get((src, "L", qb)); AH = P_Q.get((src, "H", qb))
        if AL is None or AH is None:
            continue
        TL += AL; TH += AH
        if nl >= 10 and nh >= 10:
            WL += AL * (nh / nl)           # L 을 H 의 QUAL 분포로
            WH += AH
    def _a(Z):
        m = Z[:, 0] >= MINCOV
        return float((Z[m, 1] / Z[m, 0]).mean()) if m.sum() > 30 else float("nan")
    raw = _a(TL) - _a(TH); adj = _a(WL) - _a(WH)
    print("%-7s %11.4f %11.4f %8.0f%%" % (s6, raw, adj, 100 * adj / raw))
    rows.append(("qual", s6, "matched", 0, round(raw, 5), round(adj, 5),
                 round(100 * adj / raw, 1), ""))

with io.open(os.path.join(W, "r20_actuation_unified.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("kind\tplatform\tkey\tn\tv_L\tv_H\tdiff\tp\n")
    for x in rows:
        g.write("\t".join(str(q) for q in x) + "\n")

# ---------- 자리별 윤곽
prof = []
for src, s6 in SRCS:
    for st in ("L", "H"):
        T = np.zeros((NB, 2))
        for k, A in P_SMP.items():
            if k[0] == src and k[2] == st:
                T += A
        for i in range(1, ULEN + 1):
            if T[i, 0] >= 200:
                prof.append((s6, st, i, int(T[i, 0]), T[i, 1] / T[i, 0]))
with io.open(os.path.join(W, "r20_actuation_profile.tsv"), "w", encoding="utf8", newline="\n") as g:
    g.write("platform\tstate\tunit_pos\tn_mol\tactuation\n")
    for r in prof:
        g.write("%s\t%s\t%d\t%d\t%.5f\n" % r)
print("\n저장 : r20_actuation_unified.tsv · r20_actuation_profile.tsv")
