# -*- coding: utf-8 -*-
"""뉴클레오솜을 걷어낸 판 — 세포주 전체 풀 · 유전자부 저/고 (2026-10-14)
The definitions are those used for the nucleosome-free HMM figure.
 (나) 자리 단위 : 구역창에 nuc 가 겹치지 않는 유닛만 · hit = 중심이 창 안 + 크기 범위(TFIIIB 30~59 · ICR 20~49)
 (다) 구간 대칭 : flank/core 는 core + flank 를 다 덮는 구간(c0-35..c0+27)에 nuc 가 없는 유닛만 · 유닛째 제외
 f = 그 깨끗한 유닛의 HMM 점유 · 천장 1/(1-f)
Also : whether the nucleosome-free occupancy exceeds a background window (same width, same size range, units without a nucleosome), and the low versus high difference (Fisher plus a read-cluster bootstrap)
usage: python cell_nucfree.py HepG2|K562all|GM12878all
"""
import io
import os
import sys
import json
import numpy as np
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cell_sig_lib as S

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
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


DS = {"HepG2": (_need("RDNA_HEPG2_DIR"), "HEPG2"), "K562all": (_need("RDNA_CELLS_DIR"), "K562all"),
      "GM12878all": (_need("RDNA_CELLS_DIR"), "GM12878all")}
NAME = sys.argv[1]
DIR, LAB = DS[NAME]
LO, HI, NP = S.LO, S.HI, S.NP
SITE = {"TFIIIB": dict(win=(-28, -10), size=(30, 59), c0=-11), "ICR": dict(win=(57, 68), size=(20, 49), c0=58)}
rng = np.random.default_rng(20261014)
NBOOT = 2000

keys, reads = [], []
with io.open(os.path.join(DIR, "UNIT_%s.tsv" % LAB), encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(NL).split(TB)
        keys.append(x[2]); reads.append(x[3])
N = len(keys)
ix = dict((k, i) for i, k in enumerate(keys))
rid = {}
read_i = np.array([rid.setdefault(r, len(rid)) for r in reads])
NR = len(rid)
NUC = np.zeros((N, NP), np.uint8)
TF = {s: np.zeros((N, NP), np.uint8) for s in SITE}
with io.open(os.path.join(DIR, "FHMM_UNITS_%s.tsv" % LAB), encoding="utf8") as f:
    f.readline()
    for l in f:
        x = l.rstrip(NL).split(TB)
        i = ix.get(x[2])
        if i is None:
            continue
        s0, e0, sz = int(x[5]), int(x[6]), int(x[7])
        if x[4] == "nuc":
            a, b = max(LO, s0) - LO, min(HI, e0) - LO
            if b >= a:
                NUC[i, a:b + 1] = 1
        elif x[4] == "tf":
            c = (s0 + e0) // 2
            if LO <= c <= HI:
                for s, d in SITE.items():
                    if d["size"][0] <= sz <= d["size"][1]:
                        TF[s][i, c - LO] = 1
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
            if LO <= r <= HI:
                A[i, r - LO] += 1
        for v in (x[4].split(",") if len(x) > 4 and x[4] else []):
            r = int(v)
            if LO <= r <= HI:
                M[i, r - LO] += 1
lab = np.array(["" for _ in range(N)], dtype=object)
with io.open(os.path.join(DIR, "strata", "%s_gene.tsv" % LAB if LAB != "HEPG2" else "HepG2_gene.tsv"), encoding="utf8") as f:
    f.readline()
    for l in f:
        k, v = l.rstrip(NL).split(TB)
        if k in ix:
            lab[ix[k]] = v
SETS = {"전체": np.ones(N, bool), "low": lab == "low", "high": lab == "high"}
CN = np.concatenate([np.zeros((N, 1), np.int32), np.cumsum(NUC, 1, dtype=np.int32)], 1)
CT = {s: np.concatenate([np.zeros((N, 1), np.int32), np.cumsum(TF[s], 1, dtype=np.int32)], 1) for s in SITE}


def window_arrays(s_idx, W, site):
    """각 유닛에서 창(s_idx..s_idx+W-1) : (뉴클 없음, hit)"""
    nonuc = (CN[:, s_idx + W] - CN[:, s_idx]) == 0
    hit = (CT[site][:, s_idx + W] - CT[site][:, s_idx]) > 0
    return nonuc, hit


ARR = {}
R = {"name": NAME, "n_units": int(N), "n_reads": int(NR), "sets": {k: int(v.sum()) for k, v in SETS.items()}, "sites": {}}
print("%s : 유닛 %d · read %d · 층 %s" % (NAME, N, NR, R["sets"]))
for site, d in SITE.items():
    a, b = d["win"]
    W = b - a + 1
    s0 = a - LO
    nonuc, hit = window_arrays(s0, W, site)
    hit_nn = hit & nonuc
    res = {"W": W, "sets": {}}
    # 배경 창 : 같은 폭 · 같은 크기 범위 · 뉴클 없는 유닛 · 두 자리 근처는 뺀다
    bg_s = [s for s in range(0, NP - W + 1)
            if not any((s + LO <= bb + 10 and s + LO + W - 1 >= aa - 10) for aa, bb in (SITE["TFIIIB"]["win"], SITE["ICR"]["win"]))]
    BG = {k: [] for k in SETS}
    for s in bg_s:
        nn, ht = window_arrays(s, W, site)
        for k, m in SETS.items():
            n_nn = int((nn & m).sum())
            BG[k].append(((ht & nn & m).sum() / n_nn) if n_nn >= 30 else np.nan)
    for k, m in SETS.items():
        n_all = int(m.sum()); n_nn = int((nonuc & m).sum()); h_all = int((hit & m).sum()); h_nn = int((hit_nn & m).sum())
        occ_raw = h_all / max(n_all, 1)
        occ_nn = h_nn / max(n_nn, 1)
        bg = np.array([v for v in BG[k] if v == v])
        p_bg = (float((bg >= occ_nn).sum()) + 1) / (len(bg) + 1) if len(bg) else None
        # (다) 구간 대칭 : flank/core
        c0 = d["c0"]
        lo_i, hi_i = max(LO, c0 - 35) - LO, min(HI, c0 + 27) - LO
        clean = ((CN[:, hi_i + 1] - CN[:, lo_i]) == 0) & m
        n_clean = int(clean.sum())
        h_clean = int((hit & clean).sum())
        f = h_clean / max(n_clean, 1)
        sA = A[clean].sum(0).astype(float); sM = M[clean].sum(0).astype(float)
        sA_all = A[m].sum(0).astype(float); sM_all = M[m].sum(0).astype(float)
        res["sets"][k] = dict(n_all=n_all, n_nonuc=n_nn, nonuc_frac=n_nn / max(n_all, 1), hit_all=h_all, hit_nn=h_nn,
                              occ_raw=occ_raw, occ_nn=occ_nn, bg_n=int(len(bg)), bg_median=float(np.median(bg)) if len(bg) else None,
                              bg_p95=float(np.percentile(bg, 95)) if len(bg) else None,
                              fold_vs_bg=(occ_nn / float(np.median(bg))) if len(bg) and np.median(bg) > 0 else None, p_bg=p_bg,
                              n_clean=n_clean, hit_clean=h_clean, f=f, ceil=1.0 / (1.0 - f) if f < 1 else None,
                              dub_clean=S.dub_paper(sA, sM, c0), dub_raw=S.dub_paper(sA_all, sM_all, c0))
    # 저/고 차이
    L, Hh = SETS["low"], SETS["high"]
    out = {}
    for tag, cond in (("raw", None), ("nn", nonuc)):
        if cond is None:
            nl, hl, nh, hh = int(L.sum()), int((hit & L).sum()), int(Hh.sum()), int((hit & Hh).sum())
            per_read = np.zeros((NR, 4))
            for arr, col in ((L, 0), (hit & L, 1), (Hh, 2), (hit & Hh, 3)):
                np.add.at(per_read[:, col], read_i[arr], 1)
        else:
            nl, hl, nh, hh = int((cond & L).sum()), int((hit & cond & L).sum()), int((cond & Hh).sum()), int((hit & cond & Hh).sum())
            per_read = np.zeros((NR, 4))
            for arr, col in ((cond & L, 0), (hit & cond & L, 1), (cond & Hh, 2), (hit & cond & Hh, 3)):
                np.add.at(per_read[:, col], read_i[arr], 1)
        ratio = (hl / max(nl, 1)) / max(hh / max(nh, 1), 1e-12) if hh else float("inf")
        pF = stats.fisher_exact([[hl, nl - hl], [hh, nh - hh]])[1]
        bs = []
        for _ in range(NBOOT):
            cnt = np.bincount(rng.integers(0, NR, NR), minlength=NR)
            t = cnt @ per_read
            if t[0] > 0 and t[2] > 0 and t[3] > 0:
                bs.append((t[1] / t[0]) / (t[3] / t[2]))
            elif t[0] > 0 and t[2] > 0 and t[1] > 0:
                bs.append(float("inf"))
        bs = np.array(bs)
        bsv = np.sort(bs)
        out[tag] = dict(n_low=nl, hit_low=hl, n_high=nh, hit_high=hh, ratio=ratio, fisher_p=float(pF),
                        boot_lo=float(np.quantile(bsv, 0.025, method="lower")) if len(bsv) else None,
                        boot_hi=float(np.quantile(bsv, 0.975, method="higher")) if len(bsv) else None,
                        boot_p_le1=(float((bsv <= 1.0).sum()) + 1) / (len(bsv) + 1) if len(bsv) else None,
                        boot_n=int(len(bsv)))
    res["diff"] = out
    # 깨끗한 유닛만의 호출기 없음 확인검정 (자리별 · 전체 유닛)
    c0 = d["c0"]
    lo_i, hi_i = max(LO, c0 - 35) - LO, min(HI, c0 + 27) - LO
    cl = (CN[:, hi_i + 1] - CN[:, lo_i]) == 0
    okm = cl & (A.sum(1) >= 20)
    if okm.sum() >= 100:
        cf = S.confirm(A[okm], M[okm], 1000, rng)
        res["confirm_clean"] = dict(n=int(okm.sum()), **cf[site])
    # 군집 치환 검정용 : read 안에서 저/고 라벨만 섞는다 (read 효과를 보존)
    labi = np.where(lab == "low", 0, np.where(lab == "high", 1, -1)).astype(np.int8)
    ARR[site + "_nonuc"] = nonuc; ARR[site + "_hit"] = hit
    ARR["lab"] = labi; ARR["read"] = read_i
    R["sites"][site] = res
    print("  [%s] 창 %d bp · 전체 유닛 중 뉴클 없는 비율 %.1f%%" % (site, W, 100 * res["sets"]["전체"]["nonuc_frac"]))
    for k in ("전체", "low", "high"):
        r = res["sets"][k]
        print("    %-3s 유닛 %6d · 뉴클없음 %6d (%.0f%%) · 점유 생판 %.4f → 뉴클제거 %.4f (hit %d/%d) · 배경 중앙 %s · p_bg %s · 깨끗 %d f=%.3f 천장 %.3f flank/core %.3f(생판 %.3f)"
              % (k, r["n_all"], r["n_nonuc"], 100 * r["nonuc_frac"], r["occ_raw"], r["occ_nn"], r["hit_nn"], r["n_nonuc"],
                 ("%.4f" % r["bg_median"]) if r["bg_median"] is not None else "-", ("%.4f" % r["p_bg"]) if r["p_bg"] is not None else "-",
                 r["n_clean"], r["f"], r["ceil"] or 0, r["dub_clean"], r["dub_raw"]))
    for tag in ("raw", "nn"):
        o = out[tag]
        print("    저÷고 %-3s : 저 %d/%d · 고 %d/%d · 비 %s [군집부트 %s~%s · 비<=1 비율 %s] · 피셔 p %.1e" % (
            tag, o["hit_low"], o["n_low"], o["hit_high"], o["n_high"], ("%.2f" % o["ratio"]) if o["ratio"] != float("inf") else "∞(고 hit 0)",
            ("%.2f" % o["boot_lo"]) if o["boot_lo"] is not None else "-", ("%.2f" % o["boot_hi"]) if o["boot_hi"] is not None else "-",
            ("%.4f" % o["boot_p_le1"]) if o["boot_p_le1"] is not None else "-", o["fisher_p"]))
    if "confirm_clean" in res:
        c = res["confirm_clean"]
        print("    깨끗한 유닛만 확인검정(호출기 없음) : n %d · 최대 %.3f (자리 %+d 폭 %d) · p %.4f" % (c["n"], c["obs"], c["pos"], c["width"], c["p"]))
np.savez_compressed(os.path.join(DIR, "nucfree_arrays_%s.npz" % NAME), **ARR)
out = os.path.join(DIR, "nucfree_%s.json" % NAME)
io.open(out, "w", encoding="utf8").write(json.dumps(R, ensure_ascii=False, indent=1, default=lambda o: float(o) if hasattr(o, "item") else str(o)))
print("저장 ->", out)
