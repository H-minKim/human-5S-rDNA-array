# -*- coding: utf-8 -*-
r"""Fig 2 M and Supplementary Fig 2 G: rebuilds the array type groups on cohort zone windows

   Earlier path: the majority-type path of figx9_12_prep.py (Fig 2 M table FIGX9_12_GROUPS) and answer_x12.py
           variant windows = the 9 fixed-depth windows of v20_nxn_scan (proximal 1-3, 4-7, 8-12, 13-25, middle, and the four symmetric distal ones)
   New version: variant windows = 4 cohort zones (S24 boundaries proximal 3|4, 17|18, distal 10|11 -> proximal 1-3, proximal 4-17, middle, distal 1-10)
           a unit that falls in two zones goes to the shallower one (same rule as zone_of in VTZONE.py)
   Nothing else changes: 40 tokens (top SD of carrier rate across haplotypes, candidates = units 2-98%), at least 2 units in a window,
     Jaccard >=0.90 single linkage, group = 3 or more haplotypes, separating token = within-group >=0.80 and within-group - outside-group >=0.30
   Gate: run with X12_WINDOWS=nine, and the values must equal the three published tables (majority-type rows of X12_TYPE_GROUPS, X12_SIG_TOKENS, X12_HAP_GROUPS)
   Environment variables (no defaults): RDNA_DATA_ROOT (root that holds v11_readvar), X12_WINDOWS (nine|zone), X12_OUT (output folder)
     X12_PUB (folder with the three published tables) only for the gate
   Output: X12_TYPE_GROUPS.tsv, X12_SIG_TOKENS.tsv, X12_HAP_GROUPS.tsv, X12_CHAIN.tsv, X12_WINDOWS.tsv, X12_ZONE_LOG.txt"""
import io, os, sys, collections
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd

ROOT = os.environ["RDNA_DATA_ROOT"].rstrip("/\\") + "/"
T = ROOT + "v11_readvar/02_tables/"
MAJ = ROOT + "v11_readvar/03_figures/latest_2026-10-01_majority_reference/L2_MAJREF_variant_cols_only.parquet"
WIN = os.environ["X12_WINDOWS"]
OUTD = os.environ["X12_OUT"]
os.makedirs(OUTD, exist_ok=True)
LOG = io.open(os.path.join(OUTD, "X12_ZONE_LOG.txt"), "w", encoding="utf8", newline="\n")
CUT, NTOP, MINU = 0.90, 40, 2
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def P(s=""):
    print(s, flush=True); LOG.write(str(s) + "\n")


def win_nine(D):
    """Identical, character for character, to figx9_12_prep.py and answer_x12.py"""
    return [("근1_3", (D.ap <= 3)), ("근4_7", (D.ap >= 4) & (D.ap <= 7)),
            ("근8_12", (D.ap >= 8) & (D.ap <= 12)),
            ("근13_25", (D.ap >= 13) & (D.ap <= 25)),
            ("중앙", (D.ap > 25) & (D.ad > 25)),
            ("원13_25", (D.ad >= 13) & (D.ad <= 25)),
            ("원8_12", (D.ad >= 8) & (D.ad <= 12)),
            ("원4_7", (D.ad >= 4) & (D.ad <= 7)),
            ("원1_3", (D.ad <= 3))]


def win_zone(D):
    """4 cohort zones, on overlap the shallower one (zone_of in VTZONE.py)"""
    ap, ad = D.ap.values, D.ad.values
    pz = np.where(ap <= 3, 0, np.where(ap <= 17, 1, -1))
    dz = np.where(ad <= 10, 3, -1)
    z = np.full(len(ap), 2)
    z = np.where(pz >= 0, pz, z)
    z = np.where(dz >= 0, 3, z)
    z = np.where((pz >= 0) & (dz >= 0), np.where(ap <= ad, pz, 3), z)
    P("   구역 창 : 두 구역에 다 걸린 유닛 %d → 더 얕은 쪽" % int(((pz >= 0) & (dz >= 0)).sum()))
    return [(nm, pd.Series(z == i, index=D.index)) for i, nm in enumerate(["근1_3", "근4_17", "가운데", "원1_10"])]


WINDOWS = {"nine": win_nine, "zone": win_zone}[WIN]


def build_var(D, label):
    D = D.copy(); D["sl"] = D["sample"].astype(str) + "_" + D.slot.astype(str)
    slots = sorted(D.sl.unique())
    tok = D.varlist.fillna("-").astype(str).str.split(",")
    cnt = collections.Counter()
    for tl in tok:
        for t in tl:
            if t and t != "-": cnt[t] += 1
    cand = [t for t, c in cnt.items() if 0.02 * len(D) <= c <= 0.98 * len(D)]
    TOKM = pd.DataFrame({t: tok.map(lambda L, t=t: 1.0 if t in L else 0.0).values for t in cand}, index=D.index)
    sd = TOKM.groupby(D.sl.values).mean().std()
    TOP = sd.sort_values(ascending=False).head(NTOP).index.tolist()
    out = pd.DataFrame(index=slots)
    W = WINDOWS(D)
    for wn, m in W:
        sub = TOKM[m.values]; key = D.sl.values[m.values]
        g = sub.groupby(key); keep = g.size()[g.size() >= MINU].index
        mu = g.mean().reindex(slots); mu.loc[~mu.index.isin(keep), :] = np.nan
        for t in TOP: out["%s|VAR_%s" % (wn, t)] = mu[t].values
        P("   창 %-8s 유닛 %7d · 유닛 %d 개 이상 하플로타입 %d" % (wn, int(m.sum()), MINU, len(keep)))
    P("   [%s · 창 %s %d] 후보 토큰 %d → 쓰는 토큰 %d · 열 %d · 하플로 %d"
      % (label, WIN, len(W), len(cand), len(TOP), out.shape[1], out.shape[0]))
    return out, TOP, len(W)


def pcoa(J):
    Dm = 1.0 - np.nan_to_num(J, nan=0.0); np.fill_diagonal(Dm, 0.0)
    n = len(Dm); Jc = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * Jc @ (Dm ** 2) @ Jc; B = (B + B.T) / 2
    w, v = np.linalg.eigh(B); idx = np.argsort(-w)[:2]
    lam = np.clip(w[idx], 0, None)
    return v[:, idx] * np.sqrt(lam)


# ---------------------------------------------------------------- data (majority-type path of figx9_12_prep.py, unchanged)
L = pd.read_parquet(T + "L2_HAP_UNIT.parquet",
                    columns=["sample", "slot", "array_pos", "ap", "ad", "N", "n_read", "f_hi", "varlist"])
L = L[L.n_read >= 5].copy()
M = pd.read_parquet(MAJ, columns=["sample", "slot", "array_pos", "ap", "ad", "N", "varlist"])
M = M.merge(L[["sample", "slot", "array_pos", "n_read"]], on=["sample", "slot", "array_pos"], how="inner")
P("MAJREF 자리 %s · 하플로 %d" % (format(len(M), ","), (M["sample"].astype(str) + "_" + M.slot.astype(str)).nunique()))
V, TOP, NW = build_var(M, "다수형")

haps = list(V.index)
Xb = (V.fillna(0).values > 0).astype(np.int8)
tot = Xb.sum(1).astype(float)
inter = (Xb.astype(np.int32) @ Xb.astype(np.int32).T).astype(float)
uni = tot[:, None] + tot[None, :] - inter
J = np.where(uni > 0, inter / uni, np.nan)
par = list(range(len(haps)))


def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a


iu = np.triu_indices(len(haps), 1)
for a, b in zip(iu[0][J[iu] >= CUT], iu[1][J[iu] >= CUT]):
    ra, rb = find(a), find(b)
    if ra != rb: par[rb] = ra
lab = np.array([find(i) for i in range(len(haps))])
sz = collections.Counter(lab)
big = [k for k, v in sorted(sz.items(), key=lambda x: -x[1]) if v >= 3]
gid = {k: i + 1 for i, k in enumerate(big)}
P("자카드 중앙 %.3f (10-90%% %.3f-%.3f) · 무리(>=3) %d · 무리 안 하플로 %d · 최대 %d · 하플로 5 개 이상 무리 %d (하플로 %d)"
  % (np.nanmedian(J[iu]), np.nanpercentile(J[iu], 10), np.nanpercentile(J[iu], 90), len(big),
     sum(sz[k] for k in big), max(sz[k] for k in big), sum(1 for k in big if sz[k] >= 5),
     sum(sz[k] for k in big if sz[k] >= 5)))

# ---------------------------------------------------------------- group summary (same columns and same computation as the majority-type rows of FIGX9_12_GROUPS)
LL = L.copy(); LL["hap"] = LL["sample"].astype(str) + "_" + LL.slot.astype(str)
prox = LL[(LL.ap <= LL.ad) & (LL.ap <= 10)].groupby("hap").f_hi.mean()
dist = LL[(LL.ad < LL.ap) & (LL.ad <= 10)].groupby("hap").f_hi.mean()
midd = LL[(LL.ap > 25) & (LL.ad > 25)].groupby("hap").f_hi.mean()
POP = pd.read_csv(T + "POP_COHORT.tsv", sep="\t")
POP["hap"] = POP["sample"].astype(str) + "_" + POP.slot.astype(str)
meta = POP.set_index("hap")[["sp", "pop", "N"]]
dep = LL.groupby("hap").n_read.median()
co = pcoa(J)
EMB = pd.DataFrame([dict(hap=h, sample=h.rsplit("_", 1)[0], pc1=co[i, 0], pc2=co[i, 1], 무리=gid.get(lab[i], 0),
                         sp=meta.sp.get(h, None), N=meta.N.get(h, np.nan), 깊이=dep.get(h, np.nan),
                         근위=prox.get(h, np.nan), 가운데=midd.get(h, np.nan), 원위=dist.get(h, np.nan))
                    for i, h in enumerate(haps)])
grows = []
for g, s in EMB[EMB.무리 > 0].groupby("무리"):
    pc = s.sp.value_counts()
    grows.append(dict(기준="다수형", 무리=int(g), 하플로=len(s), 개체=s["sample"].nunique(),
                      주집단=(pc.index[0] if len(pc) else "?"),
                      주집단몫=(pc.iloc[0] / pc.sum()) if len(pc) else np.nan, 집단수=int(len(pc)),
                      **{("n_" + p): int(pc.get(p, 0)) for p in POPS},
                      CN중앙=s.N.median(), 깊이중앙=s.깊이.median(), 근위중앙=s.근위.median(),
                      가운데중앙=s.가운데.median(), 원위중앙=s.원위.median(),
                      pc1중앙=s.pc1.median(), pc2중앙=s.pc2.median()))
GRP = pd.DataFrame(grows).sort_values(["기준", "하플로"], ascending=[True, False])

# ---------------------------------------------------------------- membership and separating tokens (answers 11 and 10 of answer_x12.py, unchanged)
MEMD = pd.DataFrame([dict(무리=gid.get(lab[i], 0), hap=h, sample=h.rsplit("_", 1)[0], slot=h.rsplit("_", 1)[1],
                          무리크기=int(sz[lab[i]]), 보유토큰수=int(Xb[i].sum())) for i, h in enumerate(haps)]
                    ).sort_values(["무리", "hap"])
VF = V.fillna(0); hasz = (VF > 0); drows = []
for k in big:
    m = np.array([lab[i] == k for i in range(len(haps))]); o = ~m
    for c in VF.columns:
        v = hasz[c].values; pin = v[m].mean(); pout = v[o].mean()
        if pin >= 0.80 and (pin - pout) >= 0.30:
            z, t = c.split("|VAR_")
            drows.append(dict(무리=gid[k], 구역=z, 토큰=t, 무리안_보유율=round(float(pin), 4),
                              무리밖_보유율=round(float(pout), 4), 차=round(float(pin - pout), 4)))
DT = pd.DataFrame(drows).sort_values(["무리", "차"], ascending=[True, False])
P("가르는 토큰 줄 %d · 무리 %d / %d" % (len(DT), DT.무리.nunique(), len(big)))

# pairwise Jaccard within groups, single-linkage chain diagnosis (answer 9)
CH = []
for k in big:
    ii = [i for i in range(len(haps)) if lab[i] == k]
    pv = J[np.ix_(ii, ii)][np.triu_indices(len(ii), 1)]
    CH.append(dict(무리=gid[k], 하플로=len(ii), 쌍=len(pv), 자카드최소=round(float(pv.min()), 4),
                   자카드중앙=round(float(np.median(pv)), 4), 문턱이상_몫=round(float((pv >= CUT).mean()), 4)))
CH = pd.DataFrame(CH)
P("무리 안 자카드 (하플로 5 개 이상) :\n" + CH[CH.하플로 >= 5].to_string(index=False))

# the two haplotypes of the same individual (answer 12)
smp = pd.Series([h.rsplit("_", 1)[0] for h in haps]); cc = smp.value_counts()
bothin = same = 0
for s_ in [s for s, c in cc.items() if c == 2]:
    ii = [i for i in range(len(haps)) if haps[i].rsplit("_", 1)[0] == s_]
    if all(lab[i] in big for i in ii):
        bothin += 1; same += int(lab[ii[0]] == lab[ii[1]])
P("개체 %d · 두 하플로타입이 둘 다 무리에 든 개체 %d · 그 중 같은 무리 %d" % (cc.size, bothin, same))

GRP.to_csv(os.path.join(OUTD, "X12_TYPE_GROUPS.tsv"), sep="\t", index=False, float_format="%.6g")
DT.to_csv(os.path.join(OUTD, "X12_SIG_TOKENS.tsv"), sep="\t", index=False)
MEMD.to_csv(os.path.join(OUTD, "X12_HAP_GROUPS.tsv"), sep="\t", index=False)
CH.to_csv(os.path.join(OUTD, "X12_CHAIN.tsv"), sep="\t", index=False)
io.open(os.path.join(OUTD, "X12_WINDOWS.tsv"), "w", encoding="utf8", newline="\n").write(
    "창\t창수\n%s\t%d\n" % (WIN, NW))

# ---------------------------------------------------------------- gate (only for nine)
if WIN == "nine" and os.environ.get("X12_PUB"):
    PB = os.environ["X12_PUB"]
    ok = True
    for nm, new, filt in (("X12_TYPE_GROUPS.tsv", GRP, lambda d: d[d.기준 == "다수형"]),
                          ("X12_SIG_TOKENS.tsv", DT, lambda d: d), ("X12_HAP_GROUPS.tsv", MEMD, lambda d: d)):
        a = filt(pd.read_csv(os.path.join(PB, nm), sep="\t")).reset_index(drop=True)
        b = pd.read_csv(os.path.join(OUTD, nm), sep="\t").reset_index(drop=True)
        same_ = a.shape == b.shape and list(a.columns) == list(b.columns)
        if same_:
            for c in a.columns:
                if a[c].dtype.kind in "fi" and b[c].dtype.kind in "fi":
                    same_ &= bool(np.allclose(a[c].values.astype(float), b[c].values.astype(float), equal_nan=True, rtol=1e-5))
                else:
                    same_ &= bool((a[c].astype(str).values == b[c].astype(str).values).all())
        P("관문 %-22s 발행 %s ↔ 다시 %s · 같음 %s" % (nm, a.shape, b.shape, same_))
        ok &= same_
    P("관문 %s" % ("통과" if ok else "실패"))
LOG.close()

# the deposited table(s) carry English column names and values (common/fig02_figS2/taben.py)
import sys as _s, os as _o
_s.path.insert(0, _o.environ.get("V11_CODE", _o.path.dirname(_o.path.dirname(_o.path.abspath(__file__)))))
import taben as _taben
for _f in ('X12_TYPE_GROUPS.tsv', 'X12_SIG_TOKENS.tsv', 'X12_HAP_GROUPS.tsv', 'X12_WINDOWS.tsv'): _taben.to_en_file(_o.path.join(OUTD, _f))
