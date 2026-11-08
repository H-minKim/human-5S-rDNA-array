# -*- coding: utf-8 -*-
r"""FiberHMM 판 표를 한 패스로 다 낸다.

 : *"hmm 으로 해도 괜찮은 내용을 철저히 따져보고, 가능한 것들은 hmm 으로 바꿔줘"*

따진 결과 — 칸마다 답이 다르다
   본문 (a)③ · (b) · (d)(e) · 보충 (a)  : 호출기를 안 쓴다 (m6A·CpG·깊이) -> 해당 없음
   본문 (c) 자리별 large MSP            : 바꾼다 (염기 피복 비율)
   본문 (f) large MSP 짝                : 바꾼다
   본문 (g)(h) 점유 짝                  : 바꾼다 (이미 HMM 이었으나 안 자른 판이었다)
   보충 (c) 날것 vs 짝 (점유)           : 바꾼다
   보충 (d) 구역별 large MSP            : 바꾼다
   보충 (b) 유닛당 뉴클레오솜 개수      : 안 바꾼다 — 아래 까닭
     FiberHMM 분절은 이웃을 묶는다(평균 283 · 95분위 786 · 최대 58,204 bp).
     중점 배정이라 한 분절이 26 유닛을 가로지르며 한 유닛에 통째로 들어간다.
     `ft` 는 분절 중앙 143 bp 로 뉴클레오솜 하나(147)와 맞는다.
     ⇒ "개수" 는 ft 분절에서만 물리적 뜻이 선다. 그 칸만 ft 로 두고 legend 에 적는다.

규칙은 발행 코드를 그대로 베낀다
   칸 배정 `slot_cell` (E=20 절대 · NB=10 상대) · LARGE=150 · DINUC=210      <- zg1
   large MSP 마스크 : 양 옆이 호출된 >150bp                                  <- zr9
   구역 : far(-1500..-1000) mid(-900..-400) up(-300..-1) gene(0..118)        <- zn1b
   개체 단위 Wilcoxon · 구역마다 상태마다 호출 염기 >= 500                    <- zn1b
   유닛 여과 : 호출 염기 >= 1800                                              <- zr7/zr9
상태 : beta <= 0.508
원본 표는 덮지 않는다 — 새 이름(`HMM_*`)으로 낸다.
"""
import io
import os
import sys
import gzip
import glob
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

MY = _need("RDNA_FHMM_DIR")
FH = os.path.join(MY, "data", "fhmm")
#  쓰는 폴더다 — 발행 입력 폴더를 가리키면 그 표를 덮는다.
HF = _need("RDNA_FIBERSEQ_DIR")
R = os.path.join(HF, "01_latest_results")
NL = chr(10)
CUT, MINCOV, LARGE, DINUC = 0.508, 1800, 150, 210
E, NB = 20, 10
NCELL = 2 * E + NB
REG = [("far", -1500, -1000), ("mid", -900, -400), ("up", -300, -1), ("gene", 0, 118)]
sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

SHIFT = {("HG02559", "2"), ("HG03672", "1")}
UNIT, _ = L.load_units()

NMAX = collections.Counter()
STATE = {}
with gzip.open(os.path.join(L.D, "a13_unit_state_gmm.tsv.gz"), "rt",
               encoding="utf8") as f:
    h = {k: i for i, k in enumerate(f.readline().rstrip(NL).split(chr(9)))}
    for l in f:
        x = l.rstrip(NL).split(chr(9))
        smp, slot, ap = x[h["sample"]], x[h["slot"]], x[h["array_pos"]]
        NMAX[(smp, slot)] = max(NMAX[(smp, slot)], int(ap))
        try:
            b = float(x[h["beta"]])
        except Exception:
            continue
        STATE[(smp, x[h["src"]], x[h["read_id"]], slot, ap)] = "L" if b <= CUT else "H"
print("상태 %s · 슬롯 %s" % (format(len(STATE), ","), format(len(NMAX), ",")))


def slot_cell(ap, n):
    if n < 2 * E + 4:
        return None
    if ap <= E:
        return ap - 1
    if ap > n - E:
        return E + NB + (E - (n - ap) - 1)
    r = (ap - E - 1) / float(n - 2 * E)
    return E + min(NB - 1, int(r * NB))


SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FH, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        h = {k: i for i, k in enumerate(f.readline().rstrip(NL).split(chr(9)))}
        for l in f:
            x = l.rstrip(NL).split(chr(9))
            if len(x) != len(h) or x[h["cls"]] == "tf":
                continue
            SEG[(x[h["sample"]], x[h["read_id"]][:30], x[h["slot"]])].append(
                (x[h["array_pos"]], x[h["cls"]],
                 int(x[h["rel_start"]]), int(x[h["rel_end"]])))
print("FHMM 묶음 %s" % format(len(SEG), ","))

din_bp = collections.Counter(); lmsp_bp = collections.Counter()
cov_bp = collections.Counter(); n_mol = collections.Counter()
n_din = collections.Counter(); n_patch = collections.Counter()
nunit_c = collections.Counter()
PAIRL = collections.defaultdict(lambda: {"L": [], "H": []})   # (plat,smp,rid) -> largeMSP
PAIRO = collections.defaultdict(lambda: {"L": [], "H": []})   # (plat,smp,rid) -> occ
ALLL = collections.defaultdict(list)       # (plat,state) -> largeMSP (날것)
ALLO = collections.defaultdict(list)       # (plat,state) -> occ (날것)
RG = collections.defaultdict(lambda: np.zeros(2))   # (plat,smp,state,reg)->[big,cov]
nread = nohmm = 0

for smp, src, rid, slot, arr, lo, us in L.iter_reads(UNIT, 0, 0):
    if src == "PACBIO_FIBERSEQ" and (smp in L.PB_EXCL or (smp, slot) in SHIFT):
        continue
    nread += 1
    plat = "ONT" if src.startswith("UCSC") else "PacBio"
    sg = SEG.get((smp, rid[:30], slot))
    if not sg:
        nohmm += 1
        continue
    cod = np.zeros(len(arr), dtype=np.int8)
    for ap, c, a0, b0 in sg:
        u = UNIT.get((smp, slot, ap))
        if u is None:
            continue
        i0 = (u[2] - 1) + a0 - lo
        i1 = (u[2] - 1) + b0 + 1 - lo
        if i1 <= 0 or i0 >= len(arr):
            continue
        cod[max(0, i0):min(len(arr), i1)] = 1 if c == "nuc" else 2
    own = np.full(len(arr), -1, dtype=np.int32)
    tss = {}
    for ap, _t in us:
        s_, e_, t_ = UNIT[(smp, slot, ap)]
        own[s_ - 1 - lo:e_ - lo] = int(ap)
        tss[int(ap)] = t_ - 1 - lo
    ch = np.flatnonzero(np.diff(cod)) + 1
    bb = np.concatenate(([0], ch, [len(cod)]))
    segs = [(bb[i], bb[i + 1], int(cod[bb[i]])) for i in range(len(bb) - 1)]
    big = np.zeros(len(cod), dtype=np.int8)
    din = np.zeros(len(cod), dtype=np.int8)
    n = NMAX.get((smp, slot), 0)
    for i in range(1, len(segs) - 1):
        a0, e0, c = segs[i]
        if c == 0 or segs[i - 1][2] == 0 or segs[i + 1][2] == 0:
            continue
        ln = e0 - a0
        if c == 2 and ln > LARGE:
            big[a0:e0] = 1
        elif c == 1 and ln > DINUC:
            din[a0:e0] = 1
        else:
            continue
        cm = own[(a0 + e0) // 2]
        if cm >= 0:
            cc = slot_cell(int(cm), n)
            if cc is not None:
                (n_din if c == 1 else n_patch)[(src, cc)] += 1
    # ── 칸별 누적 (본문 (c))
    if n >= 2 * E + 4:
        ok = (cod > 0) & (own >= 0)
        idx = np.flatnonzero(ok)
        if len(idx):
            cells = np.array([slot_cell(int(own[p]), n) if
                              slot_cell(int(own[p]), n) is not None else -1
                              for p in idx], dtype=np.int32)
            m2 = cells >= 0
            for p, cc in zip(idx[m2], cells[m2]):
                cov_bp[(src, int(cc))] += 1
                if big[p]:
                    lmsp_bp[(src, int(cc))] += 1
                if din[p]:
                    din_bp[(src, int(cc))] += 1
            for cc in set(int(c) for c in cells[m2]):
                n_mol[(src, cc)] += 1
    # ── 유닛별
    for ap, _t in us:
        s_, e_, _t2 = UNIT[(smp, slot, str(ap))]
        a0, b0 = s_ - 1 - lo, e_ - lo
        sub = cod[a0:b0]
        m = sub > 0
        cov = int(m.sum())
        if cov < MINCOV:
            continue
        st = STATE.get((smp, src, rid, slot, str(ap)))
        if st not in ("L", "H"):
            continue
        nunit_c[(plat, st)] += 1
        occ = float((sub == 1).sum()) / cov
        lm = float((big[a0:b0][m] == 1).sum()) / cov
        PAIRO[(plat, smp, rid[:30])][st].append(occ)
        PAIRL[(plat, smp, rid[:30])][st].append(lm)
        ALLO[(plat, st)].append(occ); ALLL[(plat, st)].append(lm)
        # 구역 (TSS 기준)
        t0 = tss[int(ap)]
        for nm, ra, rb in REG:
            i0, i1 = t0 + ra, t0 + rb + 1
            if i0 < 0 or i1 > len(cod):
                continue
            s2 = cod[i0:i1]
            mm = s2 > 0
            cv = int(mm.sum())
            if cv < 1:
                continue
            RG[(plat, smp, st, nm)] += (float((big[i0:i1][mm] == 1).sum()), cv)
print("read %s · FHMM 없어 건너뜀 %s · 유닛 %s"
      % (format(nread, ","), format(nohmm, ","), dict(nunit_c)))


def pair(D, plat):
    dl, dh, keys = [], [], []
    for k, v in D.items():
        if k[0] != plat or not (v["L"] and v["H"]):
            continue
        dl.append(np.mean(v["L"])); dh.append(np.mean(v["H"])); keys.append(k)
    if len(dl) < 20:
        return None
    dl = np.array(dl); dh = np.array(dh)
    try:
        p = wilcoxon(dl - dh)[1]
    except Exception:
        p = float("nan")
    return dl, dh, p, keys


# ══ ① 본문 (c) 칸별 프로파일 — 원본 스키마 그대로, beta 는 원본에서 가져온다
SRCP = os.path.join(R, "REF2_20261010c_merged_profile.tsv")
BETA = {}
ZON = {}
if os.path.exists(SRCP):
    hh = None
    for l in io.open(SRCP, encoding="utf8"):
        if l.startswith("#"):
            continue
        x = l.rstrip(NL).split(chr(9))
        if hh is None:
            hh = {k: i for i, k in enumerate(x)}
            continue
        BETA[(x[hh["플랫폼"]], int(x[hh["칸번호"]]))] = x[hh["beta"]]
        ZON[(x[hh["플랫폼"]], int(x[hh["칸번호"]]))] = x[hh["자리"]]
O1 = os.path.join(R, "HMM_20261014_merged_profile.tsv")
with io.open(O1, "w", encoding="utf8", newline=NL) as g:
    g.write("# FiberHMM 2.16.8 판 — zg1 의 칸 규칙(E=20 절대 · NB=10 상대) 그대로" + NL)
    g.write("# beta 는 CpG 라 호출기와 무관 — 원본 REF2_20261010c 에서 가져왔다" + NL)
    g.write(chr(9).join(["플랫폼", "자리", "칸번호", "dinuc_덮은%", "patch_덮은%",
                         "dinuc_개수_per_kb", "patch_개수_per_kb", "beta",
                         "분자수", "호출염기", "유닛수"]) + NL)
    for src, plat in (("UCSC_ONT_FIBERSEQ", "ONT"), ("PACBIO_FIBERSEQ", "PacBio")):
        for c in range(NCELL):
            cv = cov_bp.get((src, c), 0)
            if cv < 1:
                continue
            kb = cv / 1000.0
            g.write(chr(9).join([
                plat, ZON.get((plat, c), str(c)), str(c),
                "%.3f" % (100.0 * din_bp.get((src, c), 0) / cv),
                "%.3f" % (100.0 * lmsp_bp.get((src, c), 0) / cv),
                "%.4f" % (n_din.get((src, c), 0) / kb),
                "%.4f" % (n_patch.get((src, c), 0) / kb),
                BETA.get((plat, c), "-"), str(n_mol.get((src, c), 0)),
                str(cv), str(n_mol.get((src, c), 0))]) + NL)
print("① 적었다 ->", os.path.basename(O1))

# ══ ② 본문 (f) large MSP 짝
O2 = os.path.join(R, "HMM_20261014_largeMSP_pair_comparison.tsv")
with io.open(O2, "w", encoding="utf8", newline=NL) as g:
    g.write("# FiberHMM 판 — large MSP(>150bp) 덮은 염기 / 호출 염기 · 같은 분자 안 짝" + NL)
    g.write(chr(9).join(["플랫폼", "짝read", "짝_저메틸", "ci95_저", "짝_고메틸",
                         "ci95_고", "짝차이", "wilcoxon_p", "유닛수_저",
                         "유닛평균_저", "유닛수_고", "유닛평균_고"]) + NL)
    for plat in ("ONT", "PacBio"):
        r = pair(PAIRL, plat)
        if r is None:
            continue
        dl, dh, p, _ = r
        uL = np.array(ALLL[(plat, "L")]); uH = np.array(ALLL[(plat, "H")])
        g.write(chr(9).join([plat, str(len(dl)),
                             "%.6f" % dl.mean(),
                             "%.6f" % (1.96 * dl.std(ddof=1) / np.sqrt(len(dl))),
                             "%.6f" % dh.mean(),
                             "%.6f" % (1.96 * dh.std(ddof=1) / np.sqrt(len(dh))),
                             "%.6f" % (dl - dh).mean(), "%.3e" % p,
                             str(len(uL)), "%.6f" % uL.mean(),
                             str(len(uH)), "%.6f" % uH.mean()]) + NL)
print("② 적었다 ->", os.path.basename(O2))

# ══ ③ 본문 (g)(h) 점 — read 마다 (저,고)
O3 = os.path.join(R, "HMM_20261014_occupancy_pair_points.tsv")
with io.open(O3, "w", encoding="utf8", newline=NL) as g:
    g.write("# FiberHMM 판 — read 마다 저메틸 유닛 점유 평균 · 고메틸 유닛 점유 평균" + NL)
    g.write("# 점유 = nuc_bp/(nuc_bp+msp_bp) · 유닛 여과 호출염기 >= 1800 · beta <= 0.508" + NL)
    g.write(chr(9).join(["플랫폼", "개체", "read", "저메틸", "고메틸"]) + NL)
    for plat in ("ONT", "PacBio"):
        r = pair(PAIRO, plat)
        if r is None:
            continue
        dl, dh, p, keys = r
        print("   %-8s 짝 read %s · 저 %.4f · 고 %.4f · 차 %+.4f · p %.2e"
              % (plat, format(len(dl), ","), dl.mean(), dh.mean(),
                 (dl - dh).mean(), p))
        for (pl, smp, rid), vl, vh in zip(keys, dl, dh):
            g.write(chr(9).join([pl, smp, rid, "%.6f" % vl, "%.6f" % vh]) + NL)
print("③ 적었다 ->", os.path.basename(O3))

# ══ ④ 보충 (c) 날것 vs 짝
O4 = os.path.join(R, "HMM_20261014_raw_vs_paired.tsv")
with io.open(O4, "w", encoding="utf8", newline=NL) as g:
    g.write("# FiberHMM 판 — 점유 저메틸-고메틸 · 짝 안 지은 판 vs 짝 지은 판" + NL)
    g.write(chr(9).join(["플랫폼", "날것", "짝", "짝read", "유닛_저", "유닛_고",
                         "wilcoxon_p"]) + NL)
    for plat in ("ONT", "PacBio"):
        r = pair(PAIRO, plat)
        if r is None:
            continue
        dl, dh, p, _ = r
        aL = np.array(ALLO[(plat, "L")]); aH = np.array(ALLO[(plat, "H")])
        g.write(chr(9).join([plat, "%.6f" % (aL.mean() - aH.mean()),
                             "%.6f" % (dl - dh).mean(), str(len(dl)),
                             str(len(aL)), str(len(aH)), "%.3e" % p]) + NL)
print("④ 적었다 ->", os.path.basename(O4))

# ══ ⑤ 보충 (d) 구역별 — 개체 단위 Wilcoxon
O5 = os.path.join(R, "HMM_20261014_by_zone.tsv")
with io.open(O5, "w", encoding="utf8", newline=NL) as g:
    g.write("# FiberHMM 판 — 구역별 large MSP(>150bp) 덮은 염기 · 개체 단위 Wilcoxon" + NL)
    g.write("# 구역·여과는 zn1b 와 같다 (구역마다 상태마다 호출 염기 >= 500)" + NL)
    g.write(chr(9).join(["kind", "platform", "key", "n", "v_L", "v_H",
                         "diff", "p"]) + NL)
    for plat in ("ONT", "PacBio"):
        for nm, _a, _b in REG:
            aL, aH = [], []
            smps = {s for (p_, s, st, r_) in RG if p_ == plat and r_ == nm}
            for s in sorted(smps):
                vL = RG.get((plat, s, "L", nm)); vH = RG.get((plat, s, "H", nm))
                if vL is None or vH is None or vL[1] < 500 or vH[1] < 500:
                    continue
                aL.append(vL[0] / vL[1]); aH.append(vH[0] / vH[1])
            if len(aL) < 5:
                continue
            aL = np.array(aL); aH = np.array(aH)
            p = wilcoxon(aL, aH)[1]
            g.write(chr(9).join(["region", plat, nm, str(len(aL)),
                                 "%.5f" % aL.mean(), "%.5f" % aH.mean(),
                                 "%.5f" % (aL.mean() - aH.mean()),
                                 "%.3g" % p]) + NL)
            print("   %-8s %-5s n %2d · 저 %.5f · 고 %.5f · 차 %+.5f · p %.2e"
                  % (plat, nm, len(aL), aL.mean(), aH.mean(),
                     aL.mean() - aH.mean(), p))
print("⑤ 적었다 ->", os.path.basename(O5))
