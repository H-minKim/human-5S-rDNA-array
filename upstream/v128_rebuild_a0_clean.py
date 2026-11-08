# -*- coding: utf-8 -*-
# v128 (2026-09-23 · ④검증 보강) — 손상된 SAM 2 개를 깨끗한 사본으로 다시 계산해 수정판 표를 재구성
# Background : $V_ROOT/13_junction50k/final_SAM/node3/rj50k_final2/{HG02717,HG02886}.placed50k.final.sam is
#          기록 MD5(MD5_13.txt)와 다르다(92ce16e6→f656f570 · caf46658→37cfaf1b · 크기 동일 · 2026-09-23 디스크 사고).
# the copy in 04_placement/junction50k_v2/node3/rj50k_final2/ matches the recorded MD5, so only those two samples are recomputed from it.
#   방법 : v118(a0=ref_front) 과 똑같은 계산을 두 파일에만 돌리고, 나머지는 v118 의 파일별 중간 결과(04_logs/v118_partial)를 그대로 써서
#          최종 표를 v118 과 같은 방식(groupby 합·beta 평균)으로 다시 만든다.
#   출력 02_tables/FLANK_METH_HAP_A0_CLEAN.tsv · 04_logs/V128_REBUILD_CLEAN.txt
import io, os, glob, hashlib
import numpy as np, pandas as pd
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"; T = OUT + "/02_tables"
OFF = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11/13_junction50k/mapping_table/FLANK_OFFSET.tsv"
CLEAN = OUT + "/04_placement/junction50k_v2/node3/rj50k_final2"
PART = OUT + "/04_logs/v118_partial"
LOG = io.open(OUT + "/04_logs/V128_REBUILD_CLEAN.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
md5 = lambda f: hashlib.md5(open(f, "rb").read()).hexdigest()
WINS = [(0, 2000), (2000, 5000), (5000, 10000), (10000, 20000), (20000, 50000)]
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
BAD = {"HG02717": "92ce16e6df0d9c44f3f8bd8a9cf3ba4c", "HG02886": "caf4665861b72e2c911c31a369fb403e"}

O = pd.read_csv(OFF, sep="\t", dtype={"slot": str})
REF = {r.ref_name: (r.sample, r.slot, int(r.ref_front), int(r.ref_front) + int(r.canon_len)) for r in O.itertuples(index=False)}

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

def one(f):
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
            a0, a1 = REF[rn][2], REF[rn][3]
            r0 = pos; r1 = pos + ref_span(cig)
            if not any((a0 - hi <= r1 and r0 <= a0 - lo) or (a1 + lo <= r1 and r0 <= a1 + hi) for lo, hi in WINS): continue
            mm = ml = None
            for x in t[11:]:
                if x.startswith("MM:Z:"): mm = x[5:]
                elif x.startswith("ML:B:C,"):
                    ml = [int(v) for v in x[7:].split(",") if v != ""]
                    if not ml: ml = None
            if mm is None or ml is None: continue
            for rp, p in read_calls(flag, seq, cig, pos, mm, ml):
                if rp < a0: d = a0 - rp; side = "원위"
                elif rp >= a1: d = rp - a1 + 1; side = "근위"
                else: continue
                for lo, hi in WINS:
                    if lo <= d < hi:
                        k = (REF[rn][0], REF[rn][1], side, "%d-%dk" % (lo // 1000, hi // 1000))
                        a = acc.setdefault(k, [0, 0.0, set()])
                        a[0] += 1; a[1] += p; a[2].add(t[0]); break
    return pd.DataFrame([dict(sample=k[0], slot=k[1], side=k[2], win=k[3], n_call=v[0], beta=v[1] / max(v[0], 1), n_read=len(v[2]))
                         for k, v in acc.items()])

P("[0] 깨끗한 사본 확인")
rows = []
for smp, want in BAD.items():
    f = os.path.join(CLEAN, "%s.placed50k.final.sam" % smp)
    got = md5(f)
    P("   %s : %s (기록 %s) → %s" % (smp, got[:8], want[:8], "일치" if got == want else "불일치"))
    assert got == want, smp
    d = one(f); d["src"] = smp
    rows.append(d)
    P("      다시 계산 : 행 %d · 원위 0-2k n_call %s" % (len(d), d[(d.side == "원위") & (d.win == "0-2k")].n_call.tolist()))
NEW = pd.concat(rows, ignore_index=True)

P()
P("[1] v118 중간 결과에서 손상분을 빼고 깨끗한 값으로 바꿔 표 재구성")
parts = sorted(glob.glob(PART + "/*.tsv"))
P("   중간 파일 %d 개" % len(parts))
keep = []
for p in parts:
    b = os.path.basename(p)
    if b.split("_")[0] in BAD:
        P("   제외(손상 SAM 기반) : %s" % b); continue
    keep.append(pd.read_csv(p, sep="\t", dtype={"slot": str}))
F = pd.concat(keep + [NEW.drop(columns=["src"])], ignore_index=True)
F = F.groupby(["sample", "slot", "side", "win"], as_index=False).agg(
    n_call=("n_call", "sum"), n_read=("n_read", "sum"), beta=("beta", "mean"))
F.to_csv(T + "/FLANK_METH_HAP_A0_CLEAN.tsv", sep="\t", index=False, float_format="%.4f")
P("   저장 02_tables/FLANK_METH_HAP_A0_CLEAN.tsv · 행 %d · 하플로 %d" % (len(F), F.groupby(["sample", "slot"]).ngroups))

P()
P("[2] 참고판(mirror_fix)과 대조")
M = pd.read_csv(__import__("os").environ["REGISTRY_DIR"] + "/verify/flank_orientation/mirror_fix/02_tables/FLANK_METH_HAP.tsv", sep="\t", dtype={"slot": str})
K = ["sample", "slot", "side", "win"]
J = F.merge(M, on=K, how="outer", suffixes=("_a", "_b"), indicator=True)
P("   행 이 판 %d · 참고판 %d · 양쪽 %d · 한쪽만 %d" % (len(F), len(M), (J._merge == "both").sum(), (J._merge != "both").sum()))
P("   β 최대차 %.8f · n_call 최대차 %s · n_read 최대차 %s"
  % ((J.beta_a - J.beta_b).abs().max(), (J.n_call_a - J.n_call_b).abs().max(), (J.n_read_a - J.n_read_b).abs().max()))
d = J[(J.beta_a - J.beta_b).abs() > 1e-9]
P("   다른 행 %d %s" % (len(d), sorted(set(d["sample"])) if len(d) else ""))
P()
P("[3] 손상 SAM 으로 만든 v118 표와의 차이(참고)")
A = pd.read_csv(T + "/FLANK_METH_HAP_A0.tsv", sep="\t", dtype={"slot": str})
J2 = F.merge(A, on=K, suffixes=("_c", "_d"))
d2 = J2[(J2.beta_c - J2.beta_d).abs() > 1e-9]
P("   다른 행 %d · 시료 %s" % (len(d2), sorted(set(d2["sample"]))))
LOG.close()
