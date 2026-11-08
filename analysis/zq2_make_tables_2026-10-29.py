# -*- coding: utf-8 -*-
r"""그림이 읽을 표 셋을 한 패스로 낸다.

① `HMM_20261029_alu_by_zone.tsv` — 보충 (c) 가 읽을 표. 창 다섯.
     `far`·`mid`·`up`·`gene` 는 발행 표와 같아야 한다(검산) ·
     `alu` 는 `alu_own`(제 유닛 염기만)이다 — 나머지 넷이 `own_frac 1.0000` 이므로
     그래야 다섯 창이 같은 자 위에 선다.
② `HMM_20261029_paired_molecule_fraction.tsv` — 본문 (g) 범례의 `31.8% / 20.9%` 를
     FiberHMM 판에서 다시 잰다 (적힌 값은 ft 판 `CHECK_20261012s` 였다).
③ `HMM_20261029_label_rate_ratio.tsv` — 범례의 `2.80 배` 를 갈음한다.
     `2.80` 은 어느 표에서도 재현되지 않았다(생성기 다섯 군데에 글자로만).
     발행 표 `REF2_20261011_m6a_in_unit_by_state_cut508.tsv`(= 본문 (d)(e) 가 그리는 표)
     에서 깊이가중 평균으로 재면 ONT 2.19 · PacBio 1.45 로 플랫폼마다 다르다
     ⇒ 한 수로 적을 수 없다. 표로 내고 범례가 읽게 한다.

자·관문은 `zy7_make_hmm_tables.py` 를 그대로 베낀다 — 바꾸면 발행 넷이 안 맞는다.
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
#  저장소 안 파일로 가는 길은 `__file__` 에서 상대로 푼다 — `a_lib` 이 이 폴더에 같이 예치돼 있다.
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
#  산출도 그 폴더에 둔다 (as-ran 은 캠페인 폴더였다 — `zy7` 과 같은 자리로 모았다).
OUTD = R
NL = chr(10)
TB = chr(9)
ST = chr(9733)
CUT, MINCOV, LARGE = 0.508, 1800, 150

REG = [("alu", -2122, -1805), ("far", -1500, -1000), ("mid", -900, -400),
       ("up", -300, -1), ("gene", 0, 118)]
OWNONLY = {"alu"}                     # 이 창만 제 유닛 염기로 센다

sys.path.insert(0, _HERE)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8',
                              errors='replace', line_buffering=True)
import a_lib as L                                      # noqa: E402

SHIFT = {("HG02559", "2"), ("HG03672", "1")}
UNIT, _ = L.load_units()

STATE = {}
with gzip.open(os.path.join(L.D, "a13_unit_state_gmm.tsv.gz"), "rt",
               encoding="utf8") as f:
    h = {k: i for i, k in enumerate(f.readline().rstrip(NL).split(TB))}
    for l in f:
        x = l.rstrip(NL).split(TB)
        try:
            b = float(x[h["beta"]])
        except Exception:                                # noqa: BLE001
            continue
        STATE[(x[h["sample"]], x[h["src"]], x[h["read_id"]],
               x[h["slot"]], x[h["array_pos"]])] = "L" if b <= CUT else "H"

SEG = collections.defaultdict(list)
for p in sorted(glob.glob(os.path.join(FH, "FHMM_UNITS_*.tsv"))):
    with io.open(p, encoding="utf8") as f:
        h = {k: i for i, k in enumerate(f.readline().rstrip(NL).split(TB))}
        for l in f:
            x = l.rstrip(NL).split(TB)
            if len(x) != len(h) or x[h["cls"]] == "tf":
                continue
            SEG[(x[h["sample"]], x[h["read_id"]][:30], x[h["slot"]])].append(
                (x[h["array_pos"]], x[h["cls"]],
                 int(x[h["rel_start"]]), int(x[h["rel_end"]])))
print("상태 칸 %s · FHMM 묶음 %s" % (format(len(STATE), ","), format(len(SEG), ",")))

RG = collections.defaultdict(lambda: np.zeros(4))   # big, cov, own_cov, own_big
MOL = collections.defaultdict(set)                  # (plat,smp,rid30) -> {L,H}
nread = 0

for smp, src, rid, slot, arr, lo, us in L.iter_reads(UNIT, 0, 0):
    if src == "PACBIO_FIBERSEQ" and (smp in L.PB_EXCL or (smp, slot) in SHIFT):
        continue
    nread += 1
    plat = "ONT" if src.startswith("UCSC") else "PacBio"
    sg = SEG.get((smp, rid[:30], slot))
    if not sg:
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
    for i in range(1, len(segs) - 1):
        a0, e0, c = segs[i]
        if c == 0 or segs[i - 1][2] == 0 or segs[i + 1][2] == 0:
            continue
        if c == 2 and (e0 - a0) > LARGE:
            big[a0:e0] = 1
    for ap, _t in us:
        s_, e_, _t2 = UNIT[(smp, slot, str(ap))]
        a0, b0 = s_ - 1 - lo, e_ - lo
        sub = cod[a0:b0]
        if int((sub > 0).sum()) < MINCOV:
            continue
        st = STATE.get((smp, src, rid, slot, str(ap)))
        if st not in ("L", "H"):
            continue
        MOL[(plat, smp, rid[:30])].add(st)
        t0 = tss[int(ap)]
        for nm, ra, rb in REG:
            i0, i1 = t0 + ra, t0 + rb + 1
            if i0 < 0 or i1 > len(cod):
                continue
            s2 = cod[i0:i1]
            mm = s2 > 0
            if not mm.any():
                continue
            om = mm & (own[i0:i1] == int(ap))
            RG[(plat, smp, st, nm)] += (
                float((big[i0:i1][mm] == 1).sum()), int(mm.sum()),
                float(om.sum()), float((big[i0:i1][om] == 1).sum()))
print("read %s" % format(nread, ","))


def tally(nm):
    own = nm in OWNONLY
    ci, bi = (2, 3) if own else (1, 0)
    out = {}
    for plat in ("ONT", "PacBio"):
        aL, aH = [], []
        for s in sorted({k[1] for k in RG if k[0] == plat and k[3] == nm}):
            vL, vH = RG.get((plat, s, "L", nm)), RG.get((plat, s, "H", nm))
            if vL is None or vH is None or vL[ci] < 500 or vH[ci] < 500:
                continue
            aL.append(vL[bi] / vL[ci])
            aH.append(vH[bi] / vH[ci])
        if len(aL) < 5:
            out[plat] = None
            continue
        aL, aH = np.array(aL), np.array(aH)
        out[plat] = (len(aL), aL.mean(), aH.mean(), aL.mean() - aH.mean(),
                     wilcoxon(aL, aH)[1])
    return out


RES = collections.OrderedDict((nm, tally(nm)) for nm, _a, _b in REG)

# ── ① 구역별 (창 다섯)
O1 = os.path.join(OUTD, "HMM_20261029_alu_by_zone.tsv")
with io.open(O1, "w", encoding="utf8", newline=NL) as g:
    g.write("# FiberHMM 판 — 구역별 large MSP(>150bp) 덮은 염기 · 개체 단위 Wilcoxon" + NL)
    g.write("# zy7_make_hmm_tables.py 의 구역 절을 그대로 베끼고 alu 창을 더했다 (2026-10-29)" + NL)
    g.write("# 창(TSS 기준) alu -2122..-1805 · far -1500..-1000 · mid -900..-400 "
            "· up -300..-1 · gene 0..118" + NL)
    g.write("# alu 는 제 유닛 염기만 센다(own) — 나머지 넷은 own_frac 1.0000 이라 "
            "같은 자다. alu 를 그냥 세면 1.1%가 앞 유닛에서 온다" + NL)
    g.write(TB.join(["kind", "platform", "key", "n", "v_L", "v_H", "diff", "p"]) + NL)
    for nm, _a, _b in REG:
        for plat in ("ONT", "PacBio"):
            r = RES[nm][plat]
            if r is None:
                continue
            g.write(TB.join(["region", plat, nm, str(r[0]), "%.5f" % r[1],
                             "%.5f" % r[2], "%.5f" % r[3], "%.3g" % r[4]]) + NL)

# ── ② 짝 분자 비율 (FiberHMM 판)
O2 = os.path.join(OUTD, "HMM_20261029_paired_molecule_fraction.tsv")
PAIRPCT = {}
with io.open(O2, "w", encoding="utf8", newline=NL) as g:
    g.write("# 저메틸·고메틸 유닛을 둘 다 지나간 분자의 비율 — FiberHMM 판 (2026-10-29)" + NL)
    g.write("# 분모 = 상태가 불린 유닛을 하나 이상 가진 분자 (유닛 여과 호출염기>=1800 뒤)" + NL)
    g.write("# 분자 키 = (개체, read_id[:30], 슬롯) — zy7 과 같은 키" + NL)
    g.write("# 적혀 있던 31.8 / 20.9 는 ft 판(CHECK_20261012s · 2335 / 686)이었다" + NL)
    g.write(TB.join(["platform", "mol_all", "mol_pair", "pair_pct"]) + NL)
    for plat in ("ONT", "PacBio"):
        ks = [v for k, v in MOL.items() if k[0] == plat]
        a = len(ks)
        b = sum(1 for v in ks if len(v) == 2)
        PAIRPCT[plat] = 100.0 * b / a if a else float("nan")
        g.write(TB.join([plat, str(a), str(b), "%.4f" % PAIRPCT[plat]]) + NL)

# ── ③ 표지율 배수 (발행 표에서)
O3 = os.path.join(OUTD, "HMM_20261029_label_rate_ratio.tsv")
P3 = os.path.join(R, "REF2_20261011_m6a_in_unit_by_state_cut508.tsv")
rows3 = []
hdr3 = None
for l in io.open(P3, encoding="utf8"):
    if l.startswith("#"):
        continue
    x = l.rstrip(NL).split(TB)
    if hdr3 is None:
        hdr3 = {k: i for i, k in enumerate(x)}
        continue
    if len(x) == len(hdr3):
        rows3.append(x)
LAB = {}
with io.open(O3, "w", encoding="utf8", newline=NL) as g:
    g.write("# m6A 표지율 — 저메틸 유닛 ÷ 고메틸 유닛 (2026-10-29)" + NL)
    g.write("# 자 = 그 자리에서 m6A 인 분자 비율(%) · 깊이로 가중한 유닛 평균" + NL)
    g.write("# 읽은 표 = REF2_20261011_m6a_in_unit_by_state_cut508.tsv (본문 (d)(e) 가 그리는 표)" + NL)
    g.write("# 범례에 글자로 박혀 있던 '2.80 배' 는 어느 표에서도 재현되지 않았다 — "
            "플랫폼마다 다르므로 한 수로 적을 수 없다" + NL)
    g.write(TB.join(["platform", "m6a_L_pct", "m6a_H_pct", "ratio_L_over_H",
                     "n_pos_L", "n_pos_H"]) + NL)
    for plat in ("ONT", "PacBio"):
        o = {}
        for st in ("L", "H"):
            v = [(float(r[hdr3["m6a_분자비율_%"]]), float(r[hdr3["깊이"]]))
                 for r in rows3 if r[hdr3["플랫폼"]] == plat and r[hdr3["상태"]] == st]
            a = np.array([x for x, _ in v])
            w = np.array([y for _, y in v])
            o[st] = (float(np.average(a, weights=w)), len(a))
        LAB[plat] = o["L"][0] / o["H"][0]
        g.write(TB.join([plat, "%.4f" % o["L"][0], "%.4f" % o["H"][0],
                         "%.4f" % LAB[plat], str(o["L"][1]), str(o["H"][1])]) + NL)

# ───────────────────────────────── 검산
print()
print(ST + ST + " 검산 ① 발행 네 창이 그대로 재현되는가 (alu 는 새 창)")
pub = {}
for l in io.open(os.path.join(R, "HMM_20261014_by_zone.tsv"), encoding="utf8"):
    if l.startswith("#") or l.startswith("kind"):
        continue
    x = l.rstrip(NL).split(TB)
    if len(x) == 8:
        pub[(x[1], x[2])] = x
ok = bad = 0
for nm, _a, _b in REG:
    if nm == "alu":
        continue
    for plat in ("ONT", "PacBio"):
        r = RES[nm][plat]
        w = pub[(plat, nm)]
        got = "%s %s %s %s %s" % (r[0], "%.5f" % r[1], "%.5f" % r[2],
                                  "%.5f" % r[3], "%.3g" % r[4])
        want = "%s %s %s %s %s" % (w[3], w[4], w[5], w[6], w[7])
        ok += got == want
        bad += got != want
        print("   %-7s %-5s %-30s %s" % (plat, nm, want,
                                         "ok" if got == want else ST + "다름 " + got))
print("   ⇒ " + ST + "같음 %d · 어긋남 %d" % (ok, bad))
if bad:
    raise SystemExit(ST + ST + "발행 넷이 안 맞는다 — 멈춘다.")

print()
print(ST + ST + " 검산 ② alu 가 alu_own 인가 (앞 글 zq1 과 맞대기)")
for plat, want in (("ONT", "0.02913"), ("PacBio", "0.03454")):
    got = "%.5f" % RES["alu"][plat][3]
    print("   %-7s zq1 의 alu_own %s ↔ 이번 alu %s %s"
          % (plat, want, got, "ok" if got == want else ST + "다름"))

print()
print(ST + ST + " 검산 ③ 새 표 값")
for plat in ("ONT", "PacBio"):
    print("   %-7s 짝분자 %.2f%% (분모 %s · 짝 %s) · 표지율 배수 %.2f"
          % (plat, PAIRPCT[plat],
             format(len([1 for k in MOL if k[0] == plat]), ","),
             format(len([1 for k, v in MOL.items() if k[0] == plat and len(v) == 2]), ","),
             LAB[plat]))
print()
for p in (O1, O2, O3):
    print("낸 것 : %s" % p)
