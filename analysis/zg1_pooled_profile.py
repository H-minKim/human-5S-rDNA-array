# -*- coding: utf-8 -*-
"""모든 read 를 합친 자리별 프로파일 (2026-10-10).

 : "모든 read 합쳐서 봐보기 — 근위·원위는 절대, 가운데는 상대"
왜 : 슬롯 하나만 그리면 자리마다 분자가 5~20 개뿐이라 di-nuc·patch 트랙이 너무 뾰족하다.

자리 틀 (이 프로젝트의 기존 규칙 · unit-filter-and-bin-conventions)
   array_pos 1 = genomic distal (DJ), array_pos N = proximal (centromeric)
   · 원위 1..E      -> 절대 번호 d = ap            (끝에서 센다)
   · 근위 1..E      -> 절대 번호 p = N - ap + 1    (끝에서 센다)
   · 가운데         -> 상대 r = (ap - E) / (N - 2E)  를 NB 칸으로
   ⇒ 가로축 = [원위 E칸] [가운데 NB칸] [근위 E칸]

지표 (레퍼런스 Dubocanin 2025 식)
   · di-nucleosome 밀도 = (>210 bp 뉴클레오솜이 덮은 염기) / (호출된 염기)
   · 접근 가능 패치 밀도 = (>150 bp large MSP 가 덮은 염기) / (호출된 염기)
   · 개수 밀도(per kb) 도 같이 낸다
합치는 법 : 자리칸마다 분자를 다 모아 분자(염기) 합 / 분모(염기) 합  = pooled
   ⇒ 뾰족함이 사라진다. 분자 수(n) 를 같이 적는다.

배치 : 본표 a23 과 같다(정본 좌표 · 틈은 '.')
두 플랫폼 따로 낸다 (절대 수준은 플랫폼마다 다르다 — 맞대지 않는다)
"""
import io
import os
import sys
import gzip
import collections
import numpy as np

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

A = _need("RDNA_ACCESS_DIR")
#  쓰는 폴더다 — 발행 입력 폴더를 가리키면 그 표를 덮는다.
HF = _need("RDNA_FIBERSEQ_DIR")
OUT = os.path.join(HF, "01_latest_results", "REF2_20261010c_merged_profile.tsv")
NPZ = os.path.join(_need("RDNA_FHMM_DIR"), "zg1_profile.npz")   # 여기에 쓴다
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

UNIT, _ = L.load_units()
DINUC, LARGE = 210, 150
E = 20          # 양 끝에서 절대 번호로 볼 유닛 수 ( : 15 -> 20)
NB = 10         # 가운데 상대 칸 수

NMAX = collections.Counter()
for (smp, slot, ap) in UNIT:
    NMAX[(smp, slot)] = max(NMAX[(smp, slot)], int(ap))

# ── beta (자리칸마다)
BETA = collections.defaultdict(list)


def slot_cell(ap, n):
    """자리 -> 칸 번호 (0..E-1 원위 · E..E+NB-1 가운데 · E+NB..2E+NB-1 근위) · 없으면 None"""
    if n < 2 * E + 4:
        return None
    if ap <= E:
        return ap - 1
    if ap > n - E:
        return E + NB + (E - (n - ap) - 1)
    r = (ap - E - 1) / float(n - 2 * E)
    b = min(NB - 1, int(r * NB))
    return E + b


with gzip.open(os.path.join(A, "a13_unit_state_gmm.tsv.gz"), "rt", encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(chr(10)).split(chr(9))
        try:
            ap = int(x[4]); b = float(x[7])
        except Exception:
            continue
        n = NMAX.get((x[1], x[3]), 0)
        c = slot_cell(ap, n)
        if c is not None:
            BETA[(x[0], c)].append(b)
print("beta 칸 %d" % len(BETA))

NCELL = 2 * E + NB
din_bp = collections.Counter(); lmsp_bp = collections.Counter(); cov_bp = collections.Counter()
n_din = collections.Counter(); n_patch = collections.Counter(); n_mol = collections.Counter()

for smp, src, rid, slot, arr, lo, us in L.iter_reads(UNIT, 0, 0):
    if src == "PACBIO_FIBERSEQ" and smp in L.PB_EXCL:
        continue
    n = NMAX.get((smp, slot), 0)
    if n < 2 * E + 4:
        continue
    # 염기마다 어느 유닛인지
    own = np.full(len(arr), -1, dtype=np.int32)
    for ap, _t in us:
        s, e, _ = UNIT[(smp, slot, ap)]
        own[s - 1 - lo:e - lo] = int(ap)
    cod = np.zeros(len(arr), dtype=np.int8)
    cod[arr == ord("N")] = 1
    cod[arr == ord("M")] = 2
    ch = np.flatnonzero(np.diff(cod)) + 1
    b = np.concatenate(([0], ch, [len(cod)]))
    segs = [(b[i], b[i + 1], cod[b[i]]) for i in range(len(b) - 1)]
    seen = set()
    for i in range(1, len(segs) - 1):
        a0, e0, c = segs[i]
        if c == 0 or segs[i - 1][2] == 0 or segs[i + 1][2] == 0:
            continue
        ln = e0 - a0
        big = (c == 1 and ln > DINUC) or (c == 2 and ln > LARGE)
        if not big:
            continue
        for p in range(a0, e0):
            ap = own[p]
            if ap < 0:
                continue
            cc = slot_cell(int(ap), n)
            if cc is None:
                continue
            if c == 1:
                din_bp[(src, cc)] += 1
            else:
                lmsp_bp[(src, cc)] += 1
        cm = slot_cell(int(own[(a0 + e0) // 2]), n) if own[(a0 + e0) // 2] >= 0 else None
        if cm is not None:
            (n_din if c == 1 else n_patch)[(src, cm)] += 1
    ok = (cod > 0) & (own >= 0)
    for p in np.flatnonzero(ok):
        cc = slot_cell(int(own[p]), n)
        if cc is not None:
            cov_bp[(src, cc)] += 1
            seen.add(cc)
    for cc in seen:
        n_mol[(src, cc)] += 1

print("칸 %d (원위 %d · 가운데 %d · 근위 %d)" % (NCELL, E, NB, E))


def lab(c):
    if c < E:
        return "원위%d" % (c + 1)
    if c < E + NB:
        return "가운데 %d/%d" % (c - E + 1, NB)
    return "근위%d" % (2 * E + NB - c)


rows = []
for src in ("UCSC_ONT_FIBERSEQ", "PACBIO_FIBERSEQ"):
    plat = "ONT" if src.startswith("UCSC") else "PacBio"
    for c in range(NCELL):
        cv = cov_bp[(src, c)]
        if cv < 1000:
            continue
        bs = BETA.get((src, c), [])
        rows.append((plat, lab(c), c,
                     100.0 * din_bp[(src, c)] / cv,
                     100.0 * lmsp_bp[(src, c)] / cv,
                     1000.0 * n_din[(src, c)] / cv,
                     1000.0 * n_patch[(src, c)] / cv,
                     np.mean(bs) if len(bs) >= 20 else float("nan"),
                     n_mol[(src, c)], cv, len(bs)))

with io.open(OUT, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# 모든 read 를 합친 자리별 프로파일 (pooled) — 분자 합 / 분모 합" + chr(10))
    g.write("# 자리 : 원위 1..%d 절대 · 가운데 %d 칸 상대 · 근위 1..%d 절대 "
            "(array_pos 1 = 원위/DJ · N = 근위/동원체)" % (E, NB, E) + chr(10))
    g.write("# di-nuc = >210bp 뉴클레오솜이 덮은 염기 비율 · patch = >150bp large MSP 가 덮은 염기 비율"
            + chr(10))
    g.write(chr(9).join(["플랫폼", "자리", "칸번호", "dinuc_덮은%", "patch_덮은%",
                         "dinuc_개수_per_kb", "patch_개수_per_kb", "beta",
                         "분자수", "호출염기", "유닛수"]) + chr(10))
    for r in rows:
        g.write(chr(9).join([r[0], r[1], str(r[2]), "%.3f" % r[3], "%.3f" % r[4],
                             "%.4f" % r[5], "%.4f" % r[6],
                             ("%.4f" % r[7]) if r[7] == r[7] else "NA",
                             str(r[8]), str(r[9]), str(r[10])]) + chr(10))

print()
print("=== 합친 프로파일 (덮은 염기 %)")
for plat in ("ONT", "PacBio"):
    sel = [r for r in rows if r[0] == plat]
    if not sel:
        continue
    print("   [%s]  %-14s %8s %8s %8s %8s" % (plat, "자리", "di-nuc%", "patch%", "beta", "분자"))
    for r in sel:
        if r[2] < 5 or (E <= r[2] < E + NB and (r[2] - E) % 3 == 0) or r[2] >= 2 * E + NB - 5:
            print("          %-14s %8.3f %8.3f %8s %8d"
                  % (r[1], r[3], r[4], ("%.3f" % r[7]) if r[7] == r[7] else "NA", r[8]))
    print()

np.savez_compressed(NPZ, E=E, NB=NB,
                    **{"%s|%s" % (r[0], k): np.array([x[i] for x in rows if x[0] == r[0]])
                       for r in rows[:1] for k, i in (("cell", 2), ("din", 3), ("patch", 4),
                                                      ("dinkb", 5), ("patchkb", 6), ("beta", 7),
                                                      ("nmol", 8))})
print("적었다 ->", OUT)
