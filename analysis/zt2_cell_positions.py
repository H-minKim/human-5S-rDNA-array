# -*- coding: utf-8 -*-
"""칸마다 유닛 자리 수를 센다 — '자리 하나당 평균' 을 내려고.

사용자 지적 : *"(a) 가운데 부분에서 갯수 접힌 거 아니야? 유닛당 평균으로 해 줘."*

맞는 지적이다. 지금 (a) 의 `분자수` 는 칸이 넓을수록 커진다 —
   가운데 칸은 유닛 자리를 4.9 개 묶으므로, 그 중 하나만 걸친 분자도 1 로 센다.
   ⇒ 끝 칸(자리 1 개)과 같은 자가 아니다.

고치는 법 : 자리 하나당 평균 깊이 = (그 칸의 read x 유닛 관측 수) ÷ (그 칸이 묶은 자리 수)
   · read x 유닛 관측 수 = a13 의 줄 수 (= 합친 프로파일의 `유닛수`)
   · 자리 수 = 그 칸에 떨어지는 서로 다른 (개체, 슬롯, array_pos) 의 수
a01 트랙을 다시 안 읽어도 된다 — a13 한 번 읽기로 둘 다 나온다.
칸 규칙은 zg1 과 똑같이 둔다(E=20 · NB=10).
"""
import io
import os
import sys
import gzip
import collections

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
A = _need("RDNA_ACCESS_DIR")
HF = _need("RDNA_FIBERSEQ_DIR")
OUT = os.path.join(HF, "01_latest_results", "REF2_20261012_by_cell_positions_and_depth.tsv")
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

E, NB = 20, 10
NCELL = 2 * E + NB
UNIT, _ = L.load_units()
NMAX = collections.Counter()
for (smp, slot, ap) in UNIT:
    NMAX[(smp, slot)] = max(NMAX[(smp, slot)], int(ap))


def slot_cell(ap, n):
    """zg1 과 같은 규칙."""
    if n < 2 * E + 4:
        return None
    if ap <= E:
        return ap - 1
    if ap > n - E:
        return E + NB + (E - (n - ap) - 1)
    return E + min(NB - 1, int((ap - E - 1) / float(n - 2 * E) * NB))


def lab(c):
    if c < E:
        return "원위%d" % (c + 1)
    if c < E + NB:
        return "가운데 %d/%d" % (c - E + 1, NB)
    return "근위%d" % (2 * E + NB - c)


OBS = collections.Counter()                       # (src, cell) -> read x 유닛 관측 수
POS = collections.defaultdict(set)                # (src, cell) -> {(개체, 슬롯, 자리)}
nline = 0
with gzip.open(os.path.join(A, "a13_unit_state_gmm.tsv.gz"), "rt", encoding="utf8") as f:
    head = f.readline().rstrip(chr(10)).split(chr(9))
    print("a13 열 :", " | ".join(head))
    for l in f:
        x = l.rstrip(chr(10)).split(chr(9))
        nline += 1
        try:
            ap = int(x[4])
        except Exception:
            continue
        smp, src, slot = x[1], x[0], x[3]
        n = NMAX.get((smp, slot), 0)
        c = slot_cell(ap, n)
        if c is None:
            continue
        OBS[(src, c)] += 1
        POS[(src, c)].add((smp, slot, ap))
print("a13 줄 %d · 칸 %d" % (nline, len(OBS)))

rows = []
for src in ("UCSC_ONT_FIBERSEQ", "PACBIO_FIBERSEQ"):
    plat = "ONT" if src.startswith("UCSC") else "PacBio"
    for c in range(NCELL):
        o = OBS.get((src, c), 0)
        k = len(POS.get((src, c), ()))
        if k == 0:
            continue
        rows.append((plat, lab(c), c, o, k, o / float(k)))

with io.open(OUT, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# 칸마다 유닛 자리 수와 자리 하나당 평균 깊이 (2026-10-12)" + chr(10))
    g.write("# 관측수 = a13 의 read x 유닛 줄 수 · 자리수 = 서로 다른 (개체,슬롯,array_pos)"
            + chr(10))
    g.write("# 자리당깊이 = 관측수 / 자리수 — 칸 너비에 안 휘둘리는 자다" + chr(10))
    g.write("# 칸 규칙은 zg1 과 같다 (양 끝 %d 절대 · 가운데 %d 상대)" % (E, NB) + chr(10))
    g.write(chr(9).join(["플랫폼", "자리", "칸번호", "관측수", "자리수", "자리당깊이"]) + chr(10))
    for r in rows:
        g.write(chr(9).join([r[0], r[1], str(r[2]), str(r[3]), str(r[4]),
                             "%.2f" % r[5]]) + chr(10))

print()
print("=== 보기 — 끝 칸은 자리 1 개 · 가운데 칸은 여럿")
print("%-8s %-13s %10s %9s %11s" % ("플랫폼", "자리", "관측수", "자리수", "자리당깊이"))
for r in rows:
    if r[2] in (0, 1, 19, 20, 24, 29, 30, 48, 49):
        print("%-8s %-13s %10d %9d %11.2f" % (r[0], r[1], r[3], r[4], r[5]))
print()
print("적었다 ->", os.path.basename(OUT))
