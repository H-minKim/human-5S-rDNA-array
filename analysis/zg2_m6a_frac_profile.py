# -*- coding: utf-8 -*-
"""지표 ① · ④ — 분모가 깊이인 m6A (2026-10-10).

내가 틀렸던 것 : *"못 만든다"* 고 적었는데, 노드의 `M6A_SITE_*.tsv.gz` 가 바로 그 양이다.
   열 : sample | data_source | slot | array_pos | unit_pos | mod_strand | n | n_m6a | frac
        n     = 그 염기를 덮은 분자 수 (깊이)
        n_m6a = 그 염기가 m6A 인 분자 수
        frac  = n_m6a / n  = "그 염기에서 m6A 인 분자의 비율"
   ⇒ 레퍼런스 정의 그대로다 : *"counting total occurrences across all molecules at a bp
      and dividing by the coverage at that bp"*

우리가 전에 쓰던 `m6A 농도` 는 분모가 A 자리 수였다. 다른 양이다.

내는 것
  ① 자리별(array 안) m6A 비율 — 원위 절대 15 · 가운데 상대 10 · 근위 절대 15
  ④ 유닛 안 자리별(unit_pos) 프로파일 — 구역마다
합치는 법 : pooled (sum n_m6a / sum n) — 분자를 다 합친다
"""
import io
import os
import sys
import gzip
import glob
import collections
import numpy as np

# The data and output roots are environment variables with no default (CONVENTIONS, "the deposited copy wins").
#  변수 표는 `README.md` §2 에 있다.

def _need(k):
    try:
        return os.environ[k]
    except KeyError:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s. "
            "The deposited copy has no default; see README.md section 2."
            % (k, os.path.basename(__file__)))

W = _need("RDNA_FHMM_DIR")
SRCD = os.path.join(W, "m6a_site")
R = os.path.join(_need("RDNA_V12_CANON_DIR"), "02_canonical")
#  쓰는 폴더다 — 발행 입력 폴더를 가리키면 그 표를 덮는다.
HF = _need("RDNA_FIBERSEQ_DIR")
O1 = os.path.join(HF, "01_latest_results", "REF2_20261010d_m6a_depth_denominator_by_position.tsv")
O2 = os.path.join(HF, "01_latest_results", "REF2_20261010e_m6a_depth_denominator_in_unit.tsv")
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
E, NB = 20, 10      # : 끝은 20 까지
NCELL = 2 * E + NB

NMAX = collections.Counter()
with io.open(os.path.join(R, "MAP_V11.tsv"), encoding="utf8") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip(chr(10)).split(chr(9)))}
    for l in f:
        x = l.rstrip(chr(10)).split(chr(9))
        if x[ix["seg_type"]] == "UNIT":
            k = (x[ix["sample"]], x[ix["slot"]])
            NMAX[k] = max(NMAX[k], int(x[ix["array_pos"]]))


def cell(ap, n):
    if n < 2 * E + 4:
        return None
    if ap <= E:
        return ap - 1
    if ap > n - E:
        return E + NB + (E - (n - ap) - 1)
    return E + min(NB - 1, int((ap - E - 1) / float(n - 2 * E) * NB))


def zone(c):
    return "원위" if c < E else ("가운데" if c < E + NB else "근위")


nm = collections.Counter(); nn = collections.Counter()           # ① (plat, cell)
um = collections.Counter(); un = collections.Counter()           # ④ (plat, zone, unit_pos)
rows_read = 0
for p in sorted(glob.glob(os.path.join(SRCD, "*.tsv.gz"))):
    plat = "ONT" if "fsmeth" in os.path.basename(p) else "PacBio"
    with gzip.open(p, "rt", encoding="utf8") as g:
        h = {k: i for i, k in enumerate(g.readline().rstrip(chr(10)).split(chr(9)))}
        for l in g:
            x = l.rstrip(chr(10)).split(chr(9))
            if len(x) != len(h):
                continue
            rows_read += 1
            try:
                ap = int(x[h["array_pos"]]); up = int(x[h["unit_pos"]])
                n = int(x[h["n"]]); k = int(x[h["n_m6a"]])
            except Exception:
                continue
            if n <= 0:
                continue
            N = NMAX.get((x[h["sample"]], x[h["slot"]]), 0)
            c = cell(ap, N)
            if c is None:
                continue
            nm[(plat, c)] += k; nn[(plat, c)] += n
            z = zone(c)
            um[(plat, z, up)] += k; un[(plat, z, up)] += n
print("읽은 행 %d · 자리칸 %d" % (rows_read, len(nn)))


def lab(c):
    if c < E:
        return "원위%d" % (c + 1)
    if c < E + NB:
        return "가운데 %d/%d" % (c - E + 1, NB)
    return "근위%d" % (2 * E + NB - c)


with io.open(O1, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# 지표 ① : 그 염기에서 m6A 인 분자의 비율 (분모 = 깊이) — 레퍼런스 정의" + chr(10))
    g.write("# pooled : sum(n_m6a) / sum(n) · 원위 절대 %d · 가운데 상대 %d · 근위 절대 %d"
            % (E, NB, E) + chr(10))
    g.write(chr(9).join(["플랫폼", "자리", "칸번호", "m6a_분자비율_%", "합_n_m6a", "합_n"]) + chr(10))
    for plat in ("ONT", "PacBio"):
        for c in range(NCELL):
            if nn[(plat, c)] < 1000:
                continue
            g.write(chr(9).join([plat, lab(c), str(c),
                                 "%.4f" % (100.0 * nm[(plat, c)] / nn[(plat, c)]),
                                 str(nm[(plat, c)]), str(nn[(plat, c)])]) + chr(10))

print()
print("=== 지표 ① : m6A 인 분자 비율 (%) — 분모가 깊이")
for plat in ("ONT", "PacBio"):
    v = [(lab(c), 100.0 * nm[(plat, c)] / nn[(plat, c)]) for c in range(NCELL)
         if nn[(plat, c)] >= 1000]
    if not v:
        continue
    print("   [%s]" % plat)
    for i, (L_, y) in enumerate(v):
        if i < 5 or (E <= i < E + NB and (i - E) % 3 == 0) or i >= NCELL - 5:
            print("      %-12s %6.3f" % (L_, y))
    d = np.mean([y for L_, y in v[:E]]); m = np.mean([y for L_, y in v[E:E + NB]])
    p_ = np.mean([y for L_, y in v[E + NB:]])
    print("      ---- 원위 평균 %.3f · 가운데 %.3f · 근위 %.3f  (원위/가운데 %.2f 배 · 근위/가운데 %.2f 배)"
          % (d, m, p_, d / m, p_ / m))

with io.open(O2, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# 지표 ④ : 유닛 안 자리(unit_pos)별 m6A 분자 비율 — 구역마다 · 분모 = 깊이" + chr(10))
    g.write(chr(9).join(["플랫폼", "구역", "unit_pos", "m6a_분자비율_%", "합_n"]) + chr(10))
    for plat in ("ONT", "PacBio"):
        for z in ("원위", "가운데", "근위"):
            for up in sorted({u for (p_, z_, u) in un if p_ == plat and z_ == z}):
                n = un[(plat, z, up)]
                if n < 200:
                    continue
                g.write(chr(9).join([plat, z, str(up),
                                     "%.4f" % (100.0 * um[(plat, z, up)] / n),
                                     str(n)]) + chr(10))
print()
print("적었다 ->", os.path.basename(O1), "·", os.path.basename(O2))
