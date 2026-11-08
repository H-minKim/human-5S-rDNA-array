# -*- coding: utf-8 -*-
# Data for the new versions of FigX9 / FigX12, 2026-10-08
#   (1) redo the array-type groups on both the reference (CHM13) basis and the majority-allele basis
#   (2) produce the 2D haplotype layout (PCoA) on both bases, so similar types come close
#   (3) group <-> population / proximal, middle and distal methylation / CN / depth associations + label shuffling null
#   (4) CN-controlled version of FigX12 C
#
#   Why produce both
#     The variant tokens in `HAPFEAT.parquet` (02 tables folder) used by FigX9, X11 and X12 are on the CHM13 reference basis.
#     At 11 positions CHM13 carries the minor allele, so almost every haplotype "has" that token (563 0.99, 437 0.99 ...).
#     Such tokens raise the Jaccard of every pair together and change the Jaccard 0.90 grouping.
#     The majority-allele version was already made on 2026-10-01 (03 figures folder, latest_2026-10-01 majority-basis subfolder).
#
#   Gate (1): first reproduce the variant columns of `HAPFEAT.parquet` exactly, and stop if that fails.
#   Never write to the canonical files, all outputs are written only under new names below the paper Fig subfolder of the 03 figures folder.
import io, os, sys, hashlib, collections
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
T = (__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/02_tables/")
MAJ = (__import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/03_figures/latest_2026-10-01_majority_reference/L2_MAJREF_variant_cols_only.parquet")
md5 = lambda p: hashlib.md5(io.open(p, "rb").read()).hexdigest()
LOG = io.open(os.path.join(HERE, "FIGX9_12_PREP.txt"), "w", encoding="utf8", newline="\n")


def P(s=""):
    print(s, flush=True)
    LOG.write(str(s) + "\n")
    LOG.flush()


CUT = 0.90
NTOP = 40
MINU = 2
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

# ---------------- common: variant window (identical to v20_nxn_scan.py down to the characters) ----------------
def var_windows(D):
    return [("근1_3", (D.ap <= 3)), ("근4_7", (D.ap >= 4) & (D.ap <= 7)),
            ("근8_12", (D.ap >= 8) & (D.ap <= 12)),
            ("근13_25", (D.ap >= 13) & (D.ap <= 25)),
            ("중앙", (D.ap > 25) & (D.ad > 25)),
            ("원13_25", (D.ad >= 13) & (D.ad <= 25)),
            ("원8_12", (D.ad >= 8) & (D.ad <= 12)),
            ("원4_7", (D.ad >= 4) & (D.ad <= 7)),
            ("원1_3", (D.ad <= 3))]


def build_var(D, label):
    """Follow the variant column construction of v20_nxn_scan.py exactly"""
    D = D.copy()
    D["sl"] = D["sample"].astype(str) + "_" + D.slot.astype(str)
    slots = sorted(D.sl.unique())
    tok_lists = D.varlist.fillna("-").astype(str).str.split(",")
    cnt = collections.Counter()
    for tl in tok_lists:
        for t in tl:
            if t and t != "-":
                cnt[t] += 1
    cand = [t for t, c in cnt.items() if 0.02 * len(D) <= c <= 0.98 * len(D)]
    TOKM = pd.DataFrame({t: tok_lists.map(lambda L, t=t: 1.0 if t in L else 0.0).values
                         for t in cand}, index=D.index)
    sd = TOKM.groupby(D.sl.values).mean().std()
    TOP = sd.sort_values(ascending=False).head(NTOP).index.tolist()
    out = pd.DataFrame(index=slots)
    for wn, m in var_windows(D):
        sub = TOKM[m.values]
        key = D.sl.values[m.values]
        g = sub.groupby(key)
        keep = g.size()[g.size() >= MINU].index
        mu = g.mean().reindex(slots)
        mu.loc[~mu.index.isin(keep), :] = np.nan
        for t in TOP:
            out["%s|VAR_%s" % (wn, t)] = mu[t].values
    P("   [%s] 후보 토큰 %d → 쓰는 토큰 %d · 열 %d · 하플로 %d"
      % (label, len(cand), len(TOP), out.shape[1], out.shape[0]))
    return out, TOP


def jaccard(X):
    tot = X.sum(1).astype(float)
    inter = (X.astype(np.int32) @ X.astype(np.int32).T).astype(float)
    uni = tot[:, None] + tot[None, :] - inter
    return np.where(uni > 0, inter / uni, np.nan)


def union_find(J):
    par = list(range(len(J)))

    def find(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    iu = np.triu_indices(len(J), 1)
    for a, b in zip(iu[0][J[iu] >= CUT], iu[1][J[iu] >= CUT]):
        ra, rb = find(a), find(b)
        if ra != rb:
            par[rb] = ra
    return np.array([find(i) for i in range(len(J))])


def pcoa(J):
    """Classical multidimensional scaling (PCoA), 2 axes, with 1 - Jaccard as the distance"""
    Dm = 1.0 - np.nan_to_num(J, nan=0.0)
    np.fill_diagonal(Dm, 0.0)
    n = len(Dm)
    Jc = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * Jc @ (Dm ** 2) @ Jc
    B = (B + B.T) / 2
    w, v = np.linalg.eigh(B)
    idx = np.argsort(-w)[:2]
    lam = np.clip(w[idx], 0, None)
    co = v[:, idx] * np.sqrt(lam)
    var = lam / np.clip(np.sum(np.clip(w, 0, None)), 1e-12, None)
    return co, var


# ================================================================= [1] reference basis
P("=" * 78)
P("[1] 참조(CHM13) 기준 — 지금 FigX9·X11·X12 가 쓰는 것")
L = pd.read_parquet(T + "L2_HAP_UNIT.parquet",
                    columns=["sample", "slot", "array_pos", "ap", "ad", "N", "n_read",
                             "f_hi", "varlist"])
L = L[L.n_read >= 5].copy()
P("   L2 자리 %s · 하플로 %d · 개체 %d"
  % (format(len(L), ","), (L["sample"].astype(str) + "_" + L.slot.astype(str)).nunique(),
     L["sample"].nunique()))
VREF, TOPREF = build_var(L, "참조")

# ---------------- gate (1) ----------------
HF = pd.read_parquet(T + "HAPFEAT.parquet")
VC = [c for c in HF.columns if "|VAR_" in c]
common = sorted(set(VC) & set(VREF.columns))
P()
P("관문 ① HAPFEAT 의 변이 열 재현")
P("   HAPFEAT 변이 열 %d · 새로 만든 열 %d · 겹치는 열 %d" % (len(VC), VREF.shape[1], len(common)))
a = HF.loc[VREF.index, common].values.astype(float)
b = VREF[common].values.astype(float)
both = ~(np.isnan(a) | np.isnan(b))
dmax = np.abs(a[both] - b[both]).max() if both.any() else np.nan
nanmis = int((np.isnan(a) ^ np.isnan(b)).sum())
P("   값 최대 차이 %.8f · 결측 자리 어긋남 %d" % (dmax, nanmis))
if len(common) < 0.9 * len(VC) or dmax > 1e-9 or nanmis > 0:
    P("   재현 실패 — 멈춘다")
    sys.exit(1)
P("   ⇒ 그대로 재현된다")

# ================================================================= [2] majority-allele basis
P()
P("=" * 78)
P("[2] 다수형 기준 — 2026-10-01 에 만들어 둔 판본")
M = pd.read_parquet(MAJ, columns=["sample", "slot", "array_pos", "ap", "ad", "N", "varlist"])
nr = L[["sample", "slot", "array_pos", "n_read"]].copy()
M = M.merge(nr, on=["sample", "slot", "array_pos"], how="inner")
P("   MAJREF 자리 %s (n_read>=5 로 맞춘 뒤) · 하플로 %d"
  % (format(len(M), ","), (M["sample"].astype(str) + "_" + M.slot.astype(str)).nunique()))
VMAJ, TOPMAJ = build_var(M, "다수형")
P("   두 판에서 쓰는 토큰이 겹치는 수 %d / %d" % (len(set(TOPREF) & set(TOPMAJ)), NTOP))

# ================================================================= [3] groups and layout
P()
P("=" * 78)
P("[3] 자카드 %.2f 무리짓기 + 2 차원 배치 — 두 기준" % CUT)
OUT = {}
for nm, V in (("참조", VREF), ("다수형", VMAJ)):
    haps = list(V.index)
    X = (V.fillna(0).values > 0).astype(np.int8)
    J = jaccard(X)
    iu = np.triu_indices(len(haps), 1)
    lab = union_find(J)
    sz = collections.Counter(lab)
    big = [k for k, v in sorted(sz.items(), key=lambda x: -x[1]) if v >= 3]
    gid = {k: i + 1 for i, k in enumerate(big)}
    co, var = pcoa(J)
    OUT[nm] = dict(haps=haps, J=J, lab=lab, big=big, gid=gid, co=co, var=var, sz=sz)
    P("   [%s] 자카드 중앙 %.3f (10-90%% %.3f-%.3f) · 무리(>=3) %d · 무리 안 하플로 %d · 최대 무리 %d"
      % (nm, np.nanmedian(J[iu]), np.nanpercentile(J[iu], 10), np.nanpercentile(J[iu], 90),
         len(big), sum(sz[k] for k in big), max(sz[k] for k in big)))
    P("        PCoA 축 1 %.1f%% · 축 2 %.1f%%" % (100 * var[0], 100 * var[1]))

try:
    from sklearn.metrics import adjusted_rand_score as ari
    same = sorted(set(OUT["참조"]["haps"]) & set(OUT["다수형"]["haps"]))
    i1 = [OUT["참조"]["haps"].index(h) for h in same]
    i2 = [OUT["다수형"]["haps"].index(h) for h in same]
    A = ari(OUT["참조"]["lab"][i1], OUT["다수형"]["lab"][i2])
    P()
    P("   두 기준의 무리 나눔이 같은가 — 조정 랜드 지수 %.3f (하플로 %d)" % (A, len(same)))
    P("   (1 이면 같은 나눔 · 0 이면 우연 수준)")
except Exception as e:
    A = np.nan
    P("   (sklearn 없음 : %s)" % e)

# ================================================================= [4] associations (middle included)
P()
P("=" * 78)
P("[4] 무리 ↔ 무엇 — 집단 · 근위/가운데/원위 메틸 · CN · 깊이")
LL = L.copy()
LL["hap"] = LL["sample"].astype(str) + "_" + LL.slot.astype(str)
prox = LL[(LL.ap <= LL.ad) & (LL.ap <= 10)].groupby("hap").f_hi.mean()
dist = LL[(LL.ad < LL.ap) & (LL.ad <= 10)].groupby("hap").f_hi.mean()
midd = LL[(LL.ap > 25) & (LL.ad > 25)].groupby("hap").f_hi.mean()      # newly made
POP = pd.read_csv(T + "POP_COHORT.tsv", sep="\t")
POP["hap"] = POP["sample"].astype(str) + "_" + POP.slot.astype(str)
meta = POP.set_index("hap")[["sp", "pop", "N"]]
dep = LL.groupby("hap").n_read.median()
P("   근위1-10 %d · 가운데(양끝 25 밖) %d · 원위1-10 %d 하플로"
  % (len(prox), len(midd), len(dist)))


def cramers_v(x, y):
    ct = pd.crosstab(x, y)
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return np.nan
    chi2 = stats.chi2_contingency(ct)[0]
    n = ct.values.sum()
    return float(np.sqrt(chi2 / (n * (min(ct.shape) - 1))))


def eta2(groups, vals):
    d = pd.DataFrame(dict(g=groups, v=vals)).dropna()
    if d.g.nunique() < 2:
        return np.nan
    gm = d.v.mean()
    ssb = sum(len(s) * (s.v.mean() - gm) ** 2 for _, s in d.groupby("g"))
    sst = ((d.v - gm) ** 2).sum()
    return float(ssb / sst) if sst > 0 else np.nan


rng = np.random.default_rng(20261008)
rows = []
for nm in ("참조", "다수형"):
    o = OUT[nm]
    mem = [(h, o["gid"][k]) for h, k in zip(o["haps"], o["lab"]) if k in o["gid"]]
    hs = [h for h, _ in mem]
    gs = pd.Series([g for _, g in mem], index=hs)
    for key, kind, ser in (("집단(인종)", "cat", meta.sp.reindex(hs)),
                           ("세부집단", "cat", meta["pop"].reindex(hs)),
                           ("근위 메틸 f_hi", "num", prox.reindex(hs)),
                           ("가운데 메틸 f_hi", "num", midd.reindex(hs)),
                           ("원위 메틸 f_hi", "num", dist.reindex(hs)),
                           ("array 길이 CN", "num", meta.N.reindex(hs)),
                           ("읽기 깊이", "num", dep.reindex(hs))):
        f = cramers_v if kind == "cat" else eta2
        obs = f(gs.values, ser.values)
        nul = [f(rng.permutation(gs.values), ser.values) for _ in range(200)]
        rows.append(dict(기준=nm, 항목=key, 종류=kind, n=int(ser.notna().sum()),
                         관측=obs, 귀무중앙=float(np.nanmedian(nul)),
                         순=obs - float(np.nanmedian(nul)),
                         p=float((np.sum(np.array(nul) >= obs) + 1) / (len(nul) + 1))))
        P("   [%s] %-14s 관측 %.3f · 귀무 %.3f · 순 %+.3f · p %.3f (하플로 %d)"
          % (nm, key, obs, np.nanmedian(nul), obs - np.nanmedian(nul), rows[-1]["p"], rows[-1]["n"]))
ASSOC = pd.DataFrame(rows)

# ================================================================= [5] CN control of FigX12 C
P()
P("=" * 78)
P("[5] FigX12 C — 무리마다 근위 메틸이 다른가, CN 을 통제하면")
crows = []
for nm in ("참조", "다수형"):
    o = OUT[nm]
    mem = [(h, o["gid"][k]) for h, k in zip(o["haps"], o["lab"]) if k in o["gid"]]
    d = pd.DataFrame(mem, columns=["hap", "무리"]).set_index("hap")
    d["근위"] = prox.reindex(d.index)
    d["가운데"] = midd.reindex(d.index)
    d["원위"] = dist.reindex(d.index)
    d["CN"] = meta.N.reindex(d.index)
    d["sp"] = meta.sp.reindex(d.index)
    d = d.dropna(subset=["근위", "CN"])
    big = d.무리.value_counts()
    d = d[d.무리.isin(big[big >= 5].index)]          # Groups with 5 or more only (so the test is valid)
    for col in ("근위", "가운데", "원위"):
        v = d.dropna(subset=[col])
        if v.무리.nunique() < 2:
            continue
        raw = eta2(v.무리.values, v[col].values)
        # CN control: the same test on residuals after regressing on CN
        ok = v.dropna(subset=["CN"])
        z = np.polyfit(ok.CN.values, ok[col].values, 1)
        res = ok[col].values - np.polyval(z, ok.CN.values)
        adj = eta2(ok.무리.values, res)
        kw = stats.kruskal(*[g[col].values for _, g in v.groupby("무리")])
        kwa = stats.kruskal(*[res[ok.무리.values == g] for g in sorted(set(ok.무리))])
        crows.append(dict(기준=nm, 지표=col, 무리수=int(v.무리.nunique()), n=len(v),
                          eta2_raw=raw, eta2_CN통제=adj,
                          KW_p_raw=float(kw.pvalue), KW_p_CN통제=float(kwa.pvalue),
                          rho_CN=float(stats.spearmanr(ok.CN, ok[col]).statistic)))
        P("   [%s] %-4s 무리 %d · n %d · eta2 %.3f → CN 통제 %.3f · KW p %.4g → %.4g (CN↔메틸 rho %+.2f)"
          % (nm, col, v.무리.nunique(), len(v), raw, adj, kw.pvalue, kwa.pvalue,
             stats.spearmanr(ok.CN, ok[col]).statistic))
CTRL = pd.DataFrame(crows)

# ================================================================= save
emb = []
for nm in ("참조", "다수형"):
    o = OUT[nm]
    for i, h in enumerate(o["haps"]):
        emb.append(dict(기준=nm, hap=h, sample=h.rsplit("_", 1)[0],
                        pc1=o["co"][i, 0], pc2=o["co"][i, 1],
                        무리=o["gid"].get(o["lab"][i], 0),
                        무리크기=int(o["sz"][o["lab"][i]]),
                        sp=meta.sp.get(h, None), pop=meta["pop"].get(h, None),
                        N=meta.N.get(h, np.nan), 깊이=dep.get(h, np.nan),
                        근위=prox.get(h, np.nan), 가운데=midd.get(h, np.nan), 원위=dist.get(h, np.nan)))
EMB = pd.DataFrame(emb)

# ---------------- group summary (for FigX12 A, B, C, D) ----------------
grows = []
for nm in ("참조", "다수형"):
    d = EMB[(EMB.기준 == nm) & (EMB.무리 > 0)]
    for g, s in d.groupby("무리"):
        pc = s.sp.value_counts()
        grows.append(dict(기준=nm, 무리=int(g), 하플로=len(s), 개체=s["sample"].nunique(),
                          주집단=(pc.index[0] if len(pc) else "?"),
                          주집단몫=(pc.iloc[0] / pc.sum()) if len(pc) else np.nan,
                          집단수=int(len(pc)),
                          **{("n_" + p): int(pc.get(p, 0)) for p in POPS},
                          CN중앙=s.N.median(), 깊이중앙=s.깊이.median(),
                          근위중앙=s.근위.median(), 가운데중앙=s.가운데.median(),
                          원위중앙=s.원위.median(),
                          pc1중앙=s.pc1.median(), pc2중앙=s.pc2.median()))
GRP = pd.DataFrame(grows).sort_values(["기준", "하플로"], ascending=[True, False])
P()
P("[무리 요약] 참조 %d · 다수형 %d" % ((GRP.기준 == "참조").sum(), (GRP.기준 == "다수형").sum()))
P(GRP[GRP.기준 == "다수형"].head(8)[["무리", "하플로", "개체", "주집단", "주집단몫", "집단수",
                                   "CN중앙", "근위중앙", "가운데중앙", "원위중앙"]].round(3).to_string(index=False))

# ---------------- depth control, since depth became significant on the majority-allele basis ----------------
P()
P("=" * 78)
P("[6] 읽기 깊이 통제 — 다수형에서 깊이가 유의해졌다 (순 +0.092 · p 0.015)")
drows = []
for nm in ("참조", "다수형"):
    d = EMB[(EMB.기준 == nm) & (EMB.무리 > 0)].dropna(subset=["깊이"]).copy()
    big = d.무리.value_counts()
    d = d[d.무리.isin(big[big >= 5].index)]
    for col in ("근위", "가운데", "원위"):
        v = d.dropna(subset=[col])
        if v.무리.nunique() < 2:
            continue
        raw = eta2(v.무리.values, v[col].values)
        z = np.polyfit(np.log10(v.깊이.values), v[col].values, 1)
        res = v[col].values - np.polyval(z, np.log10(v.깊이.values))
        adj = eta2(v.무리.values, res)
        # Control CN and depth together
        A2 = np.column_stack([np.log10(v.깊이.values), v.N.values, np.ones(len(v))])
        coef, *_ = np.linalg.lstsq(A2, v[col].values, rcond=None)
        res2 = v[col].values - A2 @ coef
        adj2 = eta2(v.무리.values, res2)
        drows.append(dict(기준=nm, 지표=col, n=len(v), 무리수=int(v.무리.nunique()),
                          eta2_raw=raw, eta2_깊이통제=adj, eta2_깊이CN통제=adj2,
                          rho_깊이=float(stats.spearmanr(v.깊이, v[col]).statistic)))
        P("   [%s] %-4s n %d · eta2 %.3f → 깊이 통제 %.3f → 깊이+CN 통제 %.3f (깊이↔메틸 rho %+.2f)"
          % (nm, col, len(v), raw, adj, adj2, drows[-1]["rho_깊이"]))
DEPC = pd.DataFrame(drows)
VAR = pd.DataFrame([dict(기준=nm, 축1=OUT[nm]["var"][0], 축2=OUT[nm]["var"][1],
                         무리수=len(OUT[nm]["big"]),
                         무리안하플로=int(sum(OUT[nm]["sz"][k] for k in OUT[nm]["big"])),
                         최대무리=int(max(OUT[nm]["sz"][k] for k in OUT[nm]["big"])),
                         ARI=A) for nm in ("참조", "다수형")])

for nm, d in (("FIGX9_12_EMBED.tsv", EMB), ("FIGX9_12_ASSOC.tsv", ASSOC),
              ("FIGX9_12_CNCTRL.tsv", CTRL), ("FIGX9_12_META.tsv", VAR),
              ("FIGX9_12_GROUPS.tsv", GRP), ("FIGX9_12_DEPTHCTRL.tsv", DEPC)):
    d.to_csv(os.path.join(HERE, nm), sep="\t", index=False, float_format="%.6g")
P()
P("[저장] " + " · ".join("%s (md5 %s)" % (n, md5(os.path.join(HERE, n))[:12])
                        for n in ("FIGX9_12_EMBED.tsv", "FIGX9_12_ASSOC.tsv",
                                  "FIGX9_12_CNCTRL.tsv", "FIGX9_12_META.tsv",
                                  "FIGX9_12_GROUPS.tsv", "FIGX9_12_DEPTHCTRL.tsv")))
LOG.close()
