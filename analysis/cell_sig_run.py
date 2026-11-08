# -*- coding: utf-8 -*-
"""세 세포주 — 두 자리가 귀무에 대해 유의한가 (블라인드 + 호출기 없음) · 그림용 곡선 (2026-10-12)

usage: python cell_sig_run.py HepG2|K562|GM12878
"""
import io
import os
import sys
import json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cell_sig_lib as S

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
DS = {"HepG2": (H("~/hepg2"), "HEPG2"),
      "K562": (H("~/cells"), "K562r"),
      "GM12878": (H("~/cells"), "GM12878FS")}
NAME = sys.argv[1]
DIR, LAB = DS[NAME]
NB_BLIND = 1000
NB_CONF = 1000 if NAME == "HepG2" else 2000
NB_CURVE = 200
CUT = 0.508
HALF, PER, NP = S.HALF, S.PER, S.NP
rng = np.random.default_rng(20261012)

# ── 적재
keys, beta = [], []
with io.open(os.path.join(DIR, "UNIT_%s.tsv" % LAB), encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(NL).split(TB)
        keys.append(x[2])
        beta.append(float(x[6]) if (x[6] and int(x[4]) >= 10) else np.nan)
beta = np.array(beta)
N = len(keys)
ix = dict((k, i) for i, k in enumerate(keys))
assert len(ix) == N
mi, ce, sz = [], [], []
with io.open(os.path.join(DIR, "FHMM_UNITS_%s.tsv" % LAB), encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(NL).split(TB)
        if x[4] != "tf":
            continue
        c = (int(x[5]) + int(x[6])) // 2
        if -HALF <= c < HALF and x[2] in ix:
            mi.append(ix[x[2]]); ce.append(c + HALF); sz.append(int(x[7]))
mi, ce, sz = np.array(mi), np.array(ce), np.array(sz)
A = np.zeros((N, NP), np.int16)
M = np.zeros((N, NP), np.int16)
with io.open(os.path.join(DIR, "SITE_M6A_%s.tsv" % LAB), encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(NL).split(TB)
        i = ix.get(x[1])
        if i is None:
            continue
        for v in (x[3].split(",") if len(x) > 3 and x[3] else []):
            r = int(v)
            if S.LO <= r <= S.HI:
                A[i, r - S.LO] += 1
        for v in (x[4].split(",") if len(x) > 4 and x[4] else []):
            r = int(v)
            if S.LO <= r <= S.HI:
                M[i, r - S.LO] += 1
print("%s : 유닛 %d · tf 사건 %d · β 있는 유닛 %d" % (NAME, N, len(ce), int(np.isfinite(beta).sum())))

R = {"name": NAME, "n_units": int(N), "n_tf": int(len(ce)), "n_beta": int(np.isfinite(beta).sum())}

# ── β 이봉성
b = beta[np.isfinite(beta)]
hs, ed = np.histogram(b, bins=40, range=(0, 1))
lo_i, hi_i = np.searchsorted(ed, 0.30), np.searchsorted(ed, 0.70)
j = lo_i + int(np.argmin(hs[lo_i:hi_i]))
R["beta"] = dict(median=float(np.median(b)), trough_bin=[float(ed[j]), float(ed[j + 1])],
                 trough_n=int(hs[j]), frac_045_055=float(((b >= .45) & (b < .55)).mean()),
                 low_frac_at_cut=float((b < CUT).mean()), hist=hs.tolist())
bimodal = R["beta"]["frac_045_055"] < 0.05
R["beta"]["bimodal"] = bool(bimodal)
print("  β 중앙 %.3f · 골 %.3f~%.3f (n=%d) · 0.45~0.55 %.2f%% · 컷 0.508 -> 저 %.1f%% · 양분 %s"
      % (R["beta"]["median"], ed[j], ed[j + 1], hs[j], 100 * R["beta"]["frac_045_055"],
         100 * R["beta"]["low_frac_at_cut"], "가능" if bimodal else "불가(단봉)"))

# ── ① 블라인드 격자 스캔 (전체 분자)
bl5 = S.blind(mi, ce, sz, N, NB_BLIND, rng, minev=30, sbin=5)
bl = S.blind(mi, ce, sz, N, NB_BLIND, rng, minev=30, sbin=10)
R["blind5"] = (None if bl5 is None else dict(
    thr=bl5["thr"], nk=bl5["nk"], NB=bl5["NB"], sites=bl5["sites"],
    peaks=[list(p[:2]) + list(p[2]) for p in bl5["peaks"][:12]]))
R["blind"] = dict(thr=bl["thr"], nk=bl["nk"], NB=bl["NB"],
                  bins=[list(x) for x in bl["bins"]],
                  sites=bl["sites"], peaks=[list(p[:2]) + list(p[2]) for p in bl["peaks"][:12]])
print("  [블라인드 10bp] 쓴 칸 %d/15 · 격자 문턱(95%%) %.2f" % (bl["nk"], bl["thr"]))
if bl5 is not None:
    print("  [블라인드 5bp ] 쓴 칸 %d/29 · 격자 문턱 %.2f · " % (bl5["nk"], bl5["thr"]) + " · ".join(
        "%s p_창 %.4f p_격자 %.4f" % (n_, r_["p_win"], r_["p_grid"]) for n_, r_ in bl5["sites"].items()))
for nm, r in bl["sites"].items():
    print("    %-7s 관측 %.1f× (칸 %s · 자리 %+d) · p_창 %.4f · p_격자 %.4f · 사전칸 p %s"
          % (nm, r["obs"], "%d~%d" % r["bin"], r["pos"], r["p_win"], r["p_grid"],
             ("%.4f" % r["p_win_bin"]) if "p_win_bin" in r else "-"))
print("    상위 봉우리 : " + " · ".join("%+d(%d~%d, %.0f×)" % (p[1], p[2][0], p[2][1], p[0])
                                      for p in bl["peaks"][:8]))

# ── ② 호출기 없는 확인검정 (전체 분자)
okm = A.sum(1) >= 20
cf = S.confirm(A[okm], M[okm], NB_CONF, rng)
R["confirm"] = cf
R["confirm_n"] = int(okm.sum())
print("  [확인·호출기 없음] 분자 %d · 치환 %d" % (okm.sum(), NB_CONF))
for nm, r in cf.items():
    print("    %-7s 관측 %.3f (자리 %+d · 폭 %d) · 귀무 중앙 %.3f · 95%% %.3f · p %.4f"
          % (nm, r["obs"], r["pos"], r["width"], r["null_med"], r["null95"], r["p"]))

# ── 그림용 곡선
arr = {"rel": S.REL, "bl_ENR": bl["ENR"], "bl_bins": np.array(bl["bins"]), "bl_thr": bl["thr"]}
# (c) 10 bp 크기 무리 전체 밀도·배수
CL = [(20, 29), (30, 39), (40, 49), (50, 59)]
for k, (a0, a1) in enumerate(CL):
    s = (sz >= a0) & (sz <= a1)
    o, bs = S.dens_fold(mi[s], ce[s], N, NB_CURVE, rng)
    arr["c%d_obs" % k] = o
    arr["c%d_base" % k] = bs
# (d)(e)(f) 메틸 양분
LAYERS = []
if bimodal:
    lay = np.full(N, -1)
    lay[np.isfinite(beta) & (beta < CUT)] = 0
    lay[np.isfinite(beta) & (beta >= CUT)] = 1
    LAYERS = [("low", 0), ("high", 1)]
else:
    lay = np.where(np.isfinite(beta), 2, -1)   # 단봉 — 양분 불가, 전체 한 층
    LAYERS = [("all", 2)]
R["layers"] = {}
for nmL, li in LAYERS:
    idx = np.where(lay == li)[0]
    n_l = len(idx)
    remap = -np.ones(N, int)
    remap[idx] = np.arange(n_l)
    R["layers"][nmL] = dict(n=int(n_l))
    for bnm, (a0, a1), cpos in (("b35", (30, 39), -11), ("b40", (40, 49), 62)):
        s = (sz >= a0) & (sz <= a1) & (remap[mi] >= 0)
        o, bs = S.dens_fold(remap[mi[s]], ce[s], n_l, NB_CURVE, rng)
        arr["d_%s_%s_obs" % (nmL, bnm)] = o
        arr["d_%s_%s_base" % (nmL, bnm)] = bs
        j = int(np.where(S.REL == cpos)[0][0])
        flo = 1000.0 / (n_l * S.SM)
        fold = float(o[j] / max(bs[j], flo))
        nev = int(((ce[s] - HALF >= cpos - 5) & (ce[s] - HALF <= cpos + 5)).sum())
        R["layers"][nmL][bnm] = dict(center=cpos, fold=fold, n_events=nev)
    sA = A[idx].sum(0).astype(float)
    sM = M[idx].sum(0).astype(float)
    R["layers"][nmL]["dub_m11"] = S.dub_paper(sA, sM, -11)
    R["layers"][nmL]["dub_p58"] = S.dub_paper(sA, sM, 58)
    R["layers"][nmL]["dub_p62"] = S.dub_paper(sA, sM, 62)
    print("  [층 %s n=%d] TFIIIB 35~39 @-11 %.1f× (사건 %d) · ICR 40~44 @+62 %.1f× (사건 %d) · "
          "flank/core TFIIIB-11 %.3f · ICR+58 %.3f"
          % (nmL, n_l, R["layers"][nmL]["b35"]["fold"], R["layers"][nmL]["b35"]["n_events"],
             R["layers"][nmL]["b40"]["fold"], R["layers"][nmL]["b40"]["n_events"],
             R["layers"][nmL]["dub_m11"], R["layers"][nmL]["dub_p58"]))

out_npz = os.path.join(DIR, "sig_%s.npz" % NAME)
np.savez_compressed(out_npz, **arr)
io.open(os.path.join(DIR, "sig_%s.json" % NAME), "w", encoding="utf8").write(
    json.dumps(R, ensure_ascii=False, indent=1, default=lambda o: float(o) if hasattr(o, "item") else str(o)))
print("저장 ->", out_npz)
