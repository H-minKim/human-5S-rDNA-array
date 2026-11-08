# -*- coding: utf-8 -*-
r"""Signal comparison by layer (origin), CpG deamination deficit and methylation profile, produced in one place

 1. For each group, the null z of the three layers and the 'most distant zone' -> separates lineage-type from position-type
 2. Per position (proximal 30 ... middle ... distal 30), C>T|CpG and C>T|non-CpG and their ratio
 3. Methylation beta at the same positions (n_reads>=3, stage 2 mean)
Output: REPRO/VTSTRAT_z.tsv, VTSTRAT_cpg_profile.tsv
"""
import io, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
CODE = __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir))
sys.path.insert(0, CODE)
import v11_io as G
from vgroup import GROUPS, GI, vgroup, iscpg_factory
import origin_filter as OF

OUT = os.path.join((os.environ.get("OUTDIR") or os.environ["RDNA_OUT_DIR"]), "REPRO")
iscpg = iscpg_factory(G.CAN)
# to use the layers of the common cohort rule, the two axes are taken separately (defaults are the same as the old behavior)
#   V11_PRESET   : slot cohort that enters the panel           (primary 564, h1 557, h3 490)
#   V11_FPRESET  : cohort on which origin f is counted (same as V11_PRESET if not given)
#   => stage 2 (methylation only H3) = V11_PRESET=h3, V11_FPRESET=primary, f is not touched
PRESET = os.environ.get("V11_PRESET", "primary")
FPRESET = os.environ.get("V11_FPRESET", PRESET)
U, _, _ = G.load(PRESET, verbose=False)
MT = G.load_meth(min_reads=1, verbose=False)
HK = sorted(U)
NH = len(HK)
if FPRESET == PRESET:
    hapf, kcnt, NS, _ = OF.build(U, HK)
else:
    UF, _, _ = G.load(FPRESET, verbose=False)
    hapf, kcnt, NS, _ = OF.build(UF, sorted(UF))
print("[VTSTRAT] 판 코호트 %s 슬롯 %d · f 코호트 %s (NS %d)" % (PRESET, NH, FPRESET, NS))
OKC = {"전부": {"다개체 공유형", "배열내 다중형", "단일유닛형"},
       "공유형 뺌": {"배열내 다중형", "단일유닛형"},
       "단일유닛형만": {"단일유닛형"}}
LAYERS = ["전부", "공유형 뺌", "단일유닛형만"]

# ---------- 1. collect z of the three layers (reads tables already produced) ----------
SRC = {"전부": "VTZONE2_zones.tsv", "공유형 뺌": "VTZONE3_notshared_zones.tsv",
       "단일유닛형만": "VTZONE3_single_zones.tsv"}
Z = {}
for lay, f in SRC.items():
    L = [l.rstrip("\n").split("\t") for l in io.open(os.path.join(OUT, f), encoding="utf8")]
    H = {c: i for i, c in enumerate(L[0])}
    Z[lay] = {r[0]: (float(r[H["최대_귀무z"]]), r[H["가장_먼_구역"]],
                     [float(x) for x in r[1:5]]) for r in L[1:]}

g = io.open(os.path.join(OUT, "VTSTRAT_z.tsv"), "w", encoding="utf8", newline="\n")
g.write("무리\t" + "\t".join("z_" + l for l in LAYERS) + "\t"
        + "\t".join("구역_" + l for l in LAYERS) + "\t판정\n")
print("■ 1. 층마다 귀무 z 와 가장 먼 구역 — 계통형인가 자리형인가")
print("   %-12s %18s %18s %18s  %s" % ("무리", LAYERS[0], LAYERS[1], LAYERS[2], "판정"))
rows = []
for t in GROUPS:
    if t == "기타":
        continue
    za, zb, zc = (Z[l][t][0] for l in LAYERS)
    ra, rb, rc = (Z[l][t][1] for l in LAYERS)
    if abs(za) >= 10 and abs(zb) < 5:
        v = "계통형 (공유형이 만든 것)"
    elif abs(zb) >= 10 and abs(zc) >= 5 and rb == rc:
        v = "자리형 (드문 층에도 남음)"
    elif abs(zb) >= 5:
        v = "둘 다 있음"
    else:
        v = "약함"
    rows.append((t, za, zb, zc, ra, rb, rc, v))
    g.write("%s\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n" % (t, za, zb, zc, ra, rb, rc, v))
    print("   %-12s %+7.1f %-10s %+7.1f %-10s %+7.1f %-10s  %s"
          % (t, za, ra, zb, rb, zc, rc, v))
g.close()

# ---------- 2 and 3. CpG / non-CpG and methylation per position ----------
REG = ["P%d" % d for d in range(1, 31)] + ["MID"] + ["D%d" % d for d in range(30, 0, -1)]
RI = {r: i for i, r in enumerate(REG)}
NR = len(REG)
KEYS = ["C>T|CpG", "G>A|CpG", "C>T|비CpG", "G>A|비CpG"]
S = {lay: np.zeros((NH, NR, len(KEYS))) for lay in LAYERS}
C = np.zeros((NH, NR))
BS = np.zeros((NH, NR)); BC = np.zeros((NH, NR))
for hi, k in enumerate(HK):
    for row in U[k]:
        ap, ad = row[0], row[1]
        v = {lay: np.zeros(len(KEYS)) for lay in LAYERS}
        for x in row[4]:
            gname = vgroup(x, iscpg)
            if gname not in KEYS:
                continue
            kl = OF.klass(x, k, hapf, kcnt)
            for lay in LAYERS:
                if kl in OKC[lay]:
                    v[lay][KEYS.index(gname)] += 1
        m = MT.get((k[0], k[1], row[6]["array_pos"]))
        b = m["beta_call"] if (m and m["n_reads"] >= 3) else None
        for r in ([("P%d" % ap)] if ap <= 30 else []) + ([("D%d" % ad)] if ad <= 30 else []) \
                + (["MID"] if (ap > 30 and ad > 30) else []):
            i = RI[r]
            C[hi, i] += 1
            for lay in LAYERS:
                S[lay][hi, i] += v[lay]
            if b is not None:
                BS[hi, i] += b; BC[hi, i] += 1
Cn = np.where(C > 0, C, np.nan)
M = {lay: np.nanmean(S[lay] / Cn[:, :, None], 0) for lay in LAYERS}
BETA = np.nanmean(np.where(BC > 0, BS / np.where(BC > 0, BC, np.nan), np.nan), 0)

g = io.open(os.path.join(OUT, "VTSTRAT_cpg_profile.tsv"), "w", encoding="utf8", newline="\n")
g.write("자리\t쪽\t접합부로부터_유닛\t메틸베타\t"
        + "\t".join("%s_%s" % (l, kk) for l in LAYERS for kk in KEYS) + "\n")
for r in REG:
    i = RI[r]
    g.write("%s\t%s\t%d\t%.5f\t%s\n"
            % (r, "MID" if r == "MID" else ("prox" if r[0] == "P" else "dist"),
               -1 if r == "MID" else int(r[1:]), BETA[i],
               "\t".join("%.5f" % M[l][i, j] for l in LAYERS for j in range(len(KEYS)))))
g.close()

print("\n■ 2. CpG 와 비CpG 를 함께 — 공유형 뺀 층 (유닛당)")
print("   %-8s %9s %9s %9s %9s %9s" % ("자리", "메틸β", "C>T|CpG", "C>T|비CpG", "CpG/비CpG", "G>A|CpG"))
for r in ("P1", "P3", "P10", "P20", "MID", "D20", "D10", "D5", "D3", "D1"):
    i = RI[r]
    a = M["공유형 뺌"][i, KEYS.index("C>T|CpG")]
    b = M["공유형 뺌"][i, KEYS.index("C>T|비CpG")]
    c = M["공유형 뺌"][i, KEYS.index("G>A|CpG")]
    print("   %-8s %9.3f %9.4f %9.4f %9.3f %9.4f" % (r, BETA[i], a, b, (a / b) if b else np.nan, c))
print("\n저장 %s/VTSTRAT_z.tsv · VTSTRAT_cpg_profile.tsv" % OUT)

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('VTSTRAT_cpg_profile.tsv',): _taben.to_en_file(_o.path.join(OUT, _f))
