# -*- coding: utf-8 -*-
# v142b (2026-10-19) — 창 축을 하드콜 자 로도 낸다
#
#   왜
#     `Fig 3 (A)` 는 한 그림에 두 선을 겹쳐 그린다 :
#       · 주황 유닛 축  = `FIG1_2_UNITS.tsv` — 하드콜 비율 (Σ(b_call × n_canon) / Σ n_canon)
#       · 검정 창 축    = `FIG1_V5_A_bp.tsv` ← `JUNCTION_FINE_METH.tsv` — 확률 평균 (Σ sum_p / Σ n_call)
#     둘은 같은 자가 아니었다 (2026-10-18 실측). 원위에서 창 축이 체계적으로 +0.103 높고,
#     그 가운데 자 차이가 +0.093 · 모으는 법 차이는 +0.007 뿐이다.
#     ⇒ 창 축을 하드콜로 다시 내면 두 선이 같은 자 위에 선다. (2026-10-19)
#
#   무엇이 다른가 — v142 와 딱 한 가지
#     원래 : 창마다 `n_call`(호출 수)과 `sum_p`(확률의 합)만 쌓았다.
#     여기 : `sum_h`(ML ≥ 128 인 호출의 수) 를 함께 쌓는다. 나머지 파싱·창·좌표는 글자 그대로 같다.
#     ⇒ `beta`(확률 평균)와 `beta_hc`(하드콜 비율) 를 한 표에 둔다. 옛 열은 그대로 남는다.
#
#   옛 표를 덮지 않는다 — 새 이름으로 낸다 : `02_tables/JUNCTION_FINE_METH_HC.tsv`
#     (`JUNCTION_FINE_METH.tsv` 는 다른 그림 묶음도 읽는다. 옛 판 보관 : `02_tables/00_previous_2026-10-19_before_hardcall_window/`)
#
#   검산 : 새 표의 `n_call`·`sum_p` 가 옛 표와 같아야 한다 (같은 자료를 같은 길로 읽었으므로).
#            다르면 파싱이 달라진 것이니 멈추고 본다.
import io, os, glob, sys
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
# 자료 뿌리는 환경변수 하나로 (2026-10-21 · 저장소 공개 준비)
#   안 주면 지금까지 쓰던 자리를 그대로 쓴다 — 돌아가는 것이 달라지지 않는다.
DATA_ROOT = os.environ["DATA_ROOT"]
V_ROOT = os.environ["V_ROOT"]
OUT = DATA_ROOT; T = OUT + "/02_tables"
SAMD = V_ROOT + "/13_junction50k/final_SAM"
OFF = V_ROOT + "/13_junction50k/mapping_table/FLANK_OFFSET.tsv"
LOG = io.open(OUT + "/04_logs/V142B_FINE_HC.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
OUTMAX, INMAX = 50000, 45000
HC = 128.0 / 255.0            # 하드콜 컷 — ML ≥ 128 (공통)


def binof(d):
    a = abs(d)
    if a < 2000: w = 50
    elif a < 10000: w = 200
    elif a <= 50000: w = 1000
    else: return None
    lo = (a // w) * w
    return lo if d > 0 else -lo


O = pd.read_csv(OFF, sep="\t", dtype={"slot": str})
REF = {}
for r in O.itertuples(index=False):
    REF[r.ref_name] = (r.sample, r.slot, int(r.ref_front), int(r.ref_front) + int(r.canon_len))
P("[0] 참조 %d · 하드콜 컷 ML>=128 (p >= %.6f)" % (len(REF), HC))


def ref_span(cig):
    n = 0; num = ""
    for ch in cig:
        if ch.isdigit(): num += ch
        elif ch in "MDN=X": n += int(num); num = ""
        else: num = ""
    return n


def qr_pairs(cig, pos):
    out = []; q = 0; r = pos; num = ""
    for ch in cig:
        if ch.isdigit(): num += ch; continue
        n = int(num); num = ""
        if ch in "M=X": out.append((q, r, n)); q += n; r += n
        elif ch in "IS": q += n
        elif ch in "DN": r += n
        elif ch == "H": pass
    return out


def read_calls(flag, seq, cig, pos, mm, ml):
    rev = bool(flag & 16)
    orig = seq.translate(COMP)[::-1] if rev else seq
    probs_off = 0; target = None
    for part in mm.rstrip(";").split(";"):
        f = part.split(",")
        nskip = len(f) - 1
        if f[0].startswith("C+m"):
            target = [int(x) for x in f[1:]]; break
        probs_off += nskip
    if target is None: return []
    cpos = [i for i, ch in enumerate(orig) if ch == "C"]
    idx = []; k = -1
    for s in target:
        k += s + 1
        if k >= len(cpos): break
        idx.append(cpos[k])
    pr = ml[probs_off:probs_off + len(idx)]
    if len(pr) < len(idx): idx = idx[:len(pr)]
    keep = [j for j, i in enumerate(idx) if i + 1 < len(orig) and orig[i + 1] in "Gg"]
    idx = [idx[j] for j in keep]; pr = [pr[j] for j in keep]
    L = len(seq)
    segs = qr_pairs(cig, pos)
    qi = [(L - 1 - i) if rev else i for i in idx]
    order = np.argsort(qi); qi = np.array(qi)[order]; pv = np.array(pr, float)[order]
    res = []; j = 0
    for q0, r0, n in segs:
        while j < len(qi) and qi[j] < q0: j += 1
        while j < len(qi) and qi[j] < q0 + n:
            res.append((r0 + (qi[j] - q0), pv[j] / 255.0)); j += 1
    return res


rows = []
files = sorted(glob.glob(SAMD + "/*/*/*.sam"))
P("[1] SAM %d 개" % len(files))
import time
t0 = time.time()
for fi, f in enumerate(files):
    smp = os.path.basename(f).split(".")[0]
    acc = {}
    with io.open(f, encoding="utf8", errors="replace") as fh:
        for line in fh:
            if line[0] == "@": continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 11: continue
            rn = t[2]
            if rn not in REF: continue
            flag = int(t[1])
            if flag & 0x900: continue
            pos = int(t[3]) - 1; cig = t[5]; seq = t[9]
            if cig == "*" or seq == "*": continue
            s0, sl, a0, a1 = REF[rn]
            r0 = pos; r1 = pos + ref_span(cig)
            if r1 < a0 - OUTMAX or r0 > a1 + OUTMAX: continue
            mm = ml = None
            for x in t[11:]:
                if x.startswith("MM:Z:"): mm = x[5:]
                elif x.startswith("ML:B:C,"):
                    ml = [int(v) for v in x[7:].split(",") if v != ""]
                    if not ml: ml = None
            if mm is None or ml is None: continue
            for rp, p in read_calls(flag, seq, cig, pos, mm, ml):
                if rp < a0:
                    d, side = a0 - rp, "원위"
                elif rp >= a1:
                    d, side = rp - a1 + 1, "근위"
                else:
                    dd, dp = rp - a0 + 1, a1 - rp
                    if dd <= dp: d, side = -dd, "원위"
                    else: d, side = -dp, "근위"
                    if abs(d) > INMAX: continue
                if d > OUTMAX: continue
                b = binof(d)
                if b is None: continue
                k = (s0, sl, side, b)
                a = acc.setdefault(k, [0, 0.0, 0])
                a[0] += 1; a[1] += p
                if p >= HC: a[2] += 1          # 하드콜
    if acc:
        rows.extend([dict(sample=k[0], slot=k[1], side=k[2], b=k[3],
                          n_call=v[0], sum_p=v[1], sum_h=v[2]) for k, v in acc.items()])
    if (fi + 1) % 20 == 0:
        P("   %d/%d ... %s · 누적 행 %d · %.0f 초" % (fi + 1, len(files), smp, len(rows), time.time() - t0))

F = pd.DataFrame(rows).groupby(["sample", "slot", "side", "b"], as_index=False).agg(
    n_call=("n_call", "sum"), sum_p=("sum_p", "sum"), sum_h=("sum_h", "sum"))
F["beta"] = F.sum_p / F.n_call
F["beta_hc"] = F.sum_h / F.n_call
F.to_csv(T + "/JUNCTION_FINE_METH_HC.tsv", sep="\t", index=False, float_format="%.5f")
P("[2] 저장 02_tables/JUNCTION_FINE_METH_HC.tsv · 행 %d · 하플로 %d · %.0f 초"
  % (len(F), F.groupby(["sample", "slot"]).ngroups, time.time() - t0))

# ---------- 검산 : 옛 표와 n_call·sum_p 가 같아야 한다 ----------
Oldp = T + "/JUNCTION_FINE_METH.tsv"
if os.path.exists(Oldp):
    Old = pd.read_csv(Oldp, sep="\t", dtype={"slot": str})
    F2 = F.copy(); F2["slot"] = F2.slot.astype(str)
    m = Old.merge(F2, on=["sample", "slot", "side", "b"], suffixes=("_old", "_new"), how="outer",
                  indicator=True)
    P()
    P("[3] 옛 표와 대조 — 옛 행 %d · 새 행 %d · 양쪽 %d · 옛쪽만 %d · 새쪽만 %d"
      % (len(Old), len(F), int((m._merge == "both").sum()),
         int((m._merge == "left_only").sum()), int((m._merge == "right_only").sum())))
    b = m[m._merge == "both"]
    dn = (b.n_call_old - b.n_call_new).abs().max()
    dp = (b.sum_p_old - b.sum_p_new).abs().max()
    P("   n_call 최대차 %s · sum_p 최대차 %.5f  ⇒ %s"
      % (dn, dp, "같다 — 파싱이 그대로다" if (dn == 0 and dp < 0.02) else "다르다 — 멈추고 본다"))
    P("   하드콜 ↔ 확률 : array 안 β 중앙 hc %.4f ↔ prob %.4f"
      % (b.beta_hc.median(), b.beta_new.median() if "beta_new" in b else b.beta.median()))
LOG.close()
