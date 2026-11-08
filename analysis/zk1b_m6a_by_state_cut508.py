# -*- coding: utf-8 -*-
"""[자동생성 — `_make_cut508_variants.py` 가 zk1_m6a_by_state.py 에서 만든다]
   저·고메틸을 `beta <= 0.508` 로 가른다 (U 없음).
   유닛 안 자리별 m6A — 저메틸 / 고메틸로 가른다 (2026-10-11).

왜 다시 하나 : 앞 그림 (b) 는 구역(원위/가운데/근위)으로 갈랐는데
   사용자가 저메틸/고메틸로 보고 싶다고 했다.
왜 전에는 못 했나 : `M6A_SITE` 는 read 를 이미 합친 표라 상태를 붙일 수 없다.
이번에 쓴 것 : `M6A_READUNIT` — read_id 가 있고 `calls` 에 자리별 ML 점수가 있다.
       열 : sample|data_source|cfg|kind|slot|array_pos|read_id|strand|mod_strand|n_a|n_m6a|frac|in_final|calls
       calls = "유닛자리:ML(0~255)" 목록 · m6A 판정 = ML >= 128
       검산 : 합계가 표의 n_m6a 와 맞는지 본다

분모는 레퍼런스대로 깊이다 — 그 자리를 덮은 분자 수.
상태는 a13 GMM (read x 유닛 단위) 로 붙인다.
"""
import io
import os
import sys
import gzip
import glob
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

SRCD = os.path.join(_need("RDNA_FHMM_DIR"), "m6a_readunit")
#  쓰는 폴더다 — 발행 입력 폴더를 가리키면 그 표를 덮는다.
HF = _need("RDNA_FIBERSEQ_DIR")
OUT = os.path.join(HF, "01_latest_results", "REF2_20261011_m6a_in_unit_by_state_cut508.tsv")
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402


# ─────────────────────────────────────────────────────────────────
# Loader that splits low and high methylation at `beta <= 0.508`
#   · a13 과 같은 키를 돌려준다 : (sample, src, read_id, slot, array_pos) -> 'L'/'H'
#   · U(보류)가 없다 — 모든 유닛이 둘 중 하나로 불린다
#   · 강건성 근거 : CHECK_20261012m_기준3종_강건성.tsv (18 조합 전부 같은 부호)
CUT508 = 0.508


def _load_state_cut508():
    import gzip as _gz
    S = {}
    with _gz.open(os.path.join(L.D, "a13_unit_state_gmm.tsv.gz"), "rt",
                  encoding="utf8") as _f:
        _h = {k: i for i, k in enumerate(_f.readline().rstrip(chr(10)).split(chr(9)))}
        for _l in _f:
            _x = _l.rstrip(chr(10)).split(chr(9))
            try:
                _b = float(_x[_h["beta"]])
            except Exception:
                continue
            S[(_x[_h["sample"]], _x[_h["src"]], _x[_h["read_id"]],
               _x[_h["slot"]], _x[_h["array_pos"]])] = "L" if _b <= CUT508 else "H"
    return S
# ─────────────────────────────────────────────────────────────────

GMM = _load_state_cut508()
print("GMM %d" % len(GMM))
ML = 128

num = collections.Counter()      # (plat, state, unit_pos) -> m6A 인 분자 수
den = collections.Counter()      # (plat, state, unit_pos) -> 그 자리를 가진 분자 수
chk_m = chk_t = 0
nrow = 0
for p in sorted(glob.glob(os.path.join(SRCD, "*.tsv.gz"))):
    plat = "ONT" if "fsmeth" in os.path.basename(p) else "PacBio"
    src = "UCSC_ONT_FIBERSEQ" if plat == "ONT" else "PACBIO_FIBERSEQ"
    with gzip.open(p, "rt", encoding="utf8") as g:
        h = {k: i for i, k in enumerate(g.readline().rstrip(chr(10)).split(chr(9)))}
        for l in g:
            x = l.rstrip(chr(10)).split(chr(9))
            if len(x) != len(h):
                continue
            if "in_final" in h and x[h["in_final"]] != "Y":
                continue
            st = GMM.get((x[h["sample"]], src, x[h["read_id"]],
                          x[h["slot"]], x[h["array_pos"]]), "NA")
            if st not in ("L", "H"):
                continue
            nrow += 1
            c = x[h["calls"]]
            if not c:
                continue
            k_ = 0
            for it in c.split(","):
                if ":" not in it:
                    continue
                a_, b_ = it.split(":", 1)
                try:
                    up = int(a_); v = int(b_)
                except Exception:
                    continue
                den[(plat, st, up)] += 1
                if v >= ML:
                    num[(plat, st, up)] += 1
                    k_ += 1
            try:
                chk_m += abs(k_ - int(x[h["n_m6a"]])); chk_t += int(x[h["n_m6a"]])
            except Exception:
                pass
print("쓴 read x 유닛 %d" % nrow)
print("검산 : ML>=%d 로 센 m6A 와 표의 n_m6a 차이 합 %d / %d (%.3f%%)"
      % (ML, chk_m, chk_t, 100.0 * chk_m / max(1, chk_t)))

rows = []
for (plat, st, up), d in sorted(den.items()):
    if d < 200:
        continue
    rows.append((plat, st, up, 100.0 * num[(plat, st, up)] / d, d))
with io.open(OUT, "w", encoding="utf8", newline=chr(10)) as g:
    g.write("# 유닛 안 자리별 m6A — 분모는 깊이(그 자리를 덮은 분자 수) · 상태별" + chr(10))
    g.write("# m6A 판정 = ML >= %d · 검산 : 표의 n_m6a 와 차이 %.3f%%"
            % (ML, 100.0 * chk_m / max(1, chk_t)) + chr(10))
    g.write(chr(9).join(["플랫폼", "상태", "unit_pos", "m6a_분자비율_%", "깊이"]) + chr(10))
    for r in rows:
        g.write(chr(9).join([r[0], r[1], str(r[2]), "%.4f" % r[3], str(r[4])]) + chr(10))
print()
for plat in ("ONT", "PacBio"):
    for st in ("L", "H"):
        v = np.array([r[3] for r in rows if r[0] == plat and r[1] == st])
        if len(v):
            print("   %-7s %s : 자리 %4d · 평균 %6.2f%%"
                  % (plat, "저메틸" if st == "L" else "고메틸", len(v), v.mean()))
print()
print("적었다 ->", os.path.basename(OUT))
