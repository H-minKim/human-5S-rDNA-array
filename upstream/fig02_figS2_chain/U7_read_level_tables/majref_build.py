# -*- coding: utf-8 -*-
# Build the majority-allele-reference variant table. Rewrites the L2 (haplotype x unit) variants against the majority allele
#
#   Why : the reference unit unit_canon.fa is a single unit of the single individual CHM13, so it carries the minor allele at 11 positions
#        (alert : the PROMPT_5S_ reference-unit CHM13-single-unit alert .md file in $RDNA_DATA_ROOT\rDNA_registry, 2026-09-29)
#   Rule : exactly the rule of the hook (hook_majref_sitecustomize.py), flip only positions where the majority allele exceeds half,
#          and at such a position a unit carrying the majority allele has no variant, a unit carrying the CHM13 allele gets the inverse token
#   Pilot test passed : the 11 positions computed independently from my data match the alert, variants per unit 15.516→6.966
#
#   Version separation (user rule) : nothing is deleted
#     previous  the 2026-10-01 CHM13-reference folder under $RDNA_DATA_ROOT\v11_readvar (figure folder 03)
#     latest    the 2026-10-01 majority-reference folder under the same figure folder (NEW below)
#   Nothing is written to the canonical table folder (T2 below), the new version goes only to the latest folder
import io, os, re, sys, collections, datetime
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
T2 = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/02_tables"
NEW = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/03_figures/latest_2026-10-01_majority_reference"
os.makedirs(NEW, exist_ok=True)
LOG = io.open(os.path.join(NEW, "MAJREF_BUILD.txt"), "w", encoding="utf8", newline="\n")


def P(s=""):
    print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()


COMMON = 0.5            # same definition as in v13_L1_read_unit.py
ELEMS = ["AluY", "NTS1", "GArich", "NTS2", "coding", "NA"]


def _pos(t):
    m = re.match(r"^(\d+)", str(t))
    return int(m.group(1)) if m else None


VG = pd.read_csv(os.path.join(T2, "VAR_GROUP_MAP.tsv"), sep="\t")
REFB = {str(r.tok): (str(r.ref) if isinstance(r.ref, str) else "") for _, r in VG.iterrows()}
CL = dict(zip(VG.tok.astype(str), VG["class"].astype(str)))
EL = dict(zip(VG.tok.astype(str), VG["element"].fillna("NA").astype(str)))
CPG = dict(zip(VG.tok.astype(str), (VG.is_cpg.astype(str) == "Y").astype(int)))
COD = dict(zip(VG.tok.astype(str), (VG.is_coding.astype(str) == "Y").astype(int)))
# unmatched as defined in v13 = tokens without an annotation in the group table (matched==0, unit_freq missing)
MATCHED = set(VG.loc[VG.matched == 1, "tok"].astype(str))


def _inv(t):
    """The variant carried by a 'CHM13 allele' unit under the majority-allele reference (L2 format)"""
    m = re.match(r"^(\d+)([ACGT])$", str(t))
    if m:
        rb = REFB.get(str(t), "")
        if not rb or rb == "nan":
            raise ValueError("ref 염기 모름: %r" % t)
        return "%s%s" % (m.group(1), rb)
    m = re.match(r"^(\d+)D(\d+)$", str(t))
    if m:
        return "%sI%s" % (m.group(1), m.group(2))
    m = re.match(r"^(\d+)I(\d+)$", str(t))
    if m:
        return "%sD%s" % (m.group(1), m.group(2))
    raise ValueError("역토큰 못 만듦: %r" % t)


P("[다수형 기준 변이표] %s" % datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
KEEP = ["sample", "slot", "array_pos", "ap", "ad", "side", "varlist", "n_var", "n_var_gene",
        "n_var_common", "n_var_rare", "n_var_unmatched", "n_snv", "n_ins", "n_del",
        "n_var_cpg", "n_var_coding", "vsig", "vsig_n", "vsig_minfreq",
        "data_source", "pop", "super_pop", "N", "both_eligible"] + ["nv_" + e for e in ELEMS]
L2 = pd.read_parquet(os.path.join(T2, "L2_HAP_UNIT.parquet"))
P("   입력 L2_HAP_UNIT.parquet : 행 %s · 하플로 %d · 개체 %d"
  % (format(len(L2), ","), L2.groupby(["sample", "slot"]).ngroups, L2["sample"].nunique()))

# ---------- 1. majority allele per position ----------
alle = collections.defaultdict(collections.Counter)
per_unit = []
for s in L2.varlist.values:
    dd = {}
    if isinstance(s, str) and s not in ("-", ""):
        for t in s.split(","):
            p = _pos(t)
            if p is None:
                continue
            dd[p] = min(dd[p], t) if p in dd else t
    per_unit.append(dd)
    for p, t in dd.items():
        alle[p][t] += 1
nu = len(L2)
for p in alle:
    alle[p]["REF"] = nu - sum(alle[p].values())
FLIP = {}
for p in alle:
    t, c = alle[p].most_common(1)[0]
    if t != "REF" and c / nu > 0.5:
        FLIP[p] = t
P("   변이 불린 자리 %d · 뒤집는 자리 %d : %s"
  % (len(alle), len(FLIP), " ".join("%d(%s→%s)" % (p, FLIP[p], _inv(FLIP[p]))
                                    for p in sorted(FLIP))))
pd.DataFrame([dict(bp=p, 다수형토큰=FLIP[p], 역토큰=_inv(FLIP[p]),
                   다수형비율=alle[p][FLIP[p]] / nu,
                   그자리_대립형질수=len([k for k in alle[p] if k != "REF"]))
              for p in sorted(FLIP)]).to_csv(
    os.path.join(NEW, "MAJREF_FLIPS_this_run.tsv"), sep="\t", index=False, float_format="%.6f")

# ---------- 2. flip varlist ----------
newlist = []
for dd in per_unit:
    s2 = set()
    for p, t in dd.items():
        if p not in FLIP or t != FLIP[p]:
            s2.add(t)
    for p, m in FLIP.items():
        if p not in dd:
            s2.add(_inv(m))
    newlist.append(",".join(sorted(s2)))
L2["varlist_maj"] = newlist
P("   유닛당 변이 수 : %.3f → %.3f (%.1f%%)"
  % (L2.n_var.mean(), np.mean([len(x.split(",")) if x else 0 for x in newlist]),
     100 * (np.mean([len(x.split(",")) if x else 0 for x in newlist]) / L2.n_var.mean() - 1)))

# ---------- 3. recompute derived columns (same definitions as v13, frequencies are recounted) ----------
cnt = collections.Counter()
for x in newlist:
    if x:
        for t in x.split(","):
            cnt[t] += 1
UFQ = {t: c / nu for t, c in cnt.items()}
P("   다수형 기준 토큰 종류 %d (전 기준 %d)" % (len(cnt), len(VG)))
P("   역토큰은 그룹표에 없다 → class·element·CpG 는 짝 토큰(같은 자리 다수형)의 것을 물려준다")
INV2FWD = {_inv(m): m for m in FLIP.values()}
rec = []
for x in newlist:
    if not x:
        rec.append((0, 0, 0, 0, 0, 0, 0, 0, [0] * 6, "", 0, np.nan)); continue
    tk = x.split(",")
    ncom = nrare = nun = nsnv = nins = ndel = ncpg = ncod = 0
    ec = [0] * 6
    sig = []
    fmin = 1.0
    for t in tk:
        # Keep the definition : unmatched = a token not in the group table (same meaning as the previous version).
        #   Frequencies were recounted against the majority allele, so every token has a frequency →
        #   judging unmatched by frequency would change its meaning.
        look0 = INV2FWD.get(t, t)
        if look0 not in MATCHED:
            nun += 1
        f = UFQ.get(t, np.nan)
        if False:
            pass
        elif f >= COMMON:
            ncom += 1
        else:
            nrare += 1; sig.append(t); fmin = min(fmin, f)
        look = INV2FWD.get(t, t)          # for an inverse token, use the annotation of its partner token
        c = CL.get(look, "")
        if c == "snv":
            nsnv += 1
        elif c == "ins":
            nins += 1
        elif c == "del":
            ndel += 1
        # an inverse token swaps del↔ins
        if t in INV2FWD:
            if c == "del":
                ndel -= 1; nins += 1
            elif c == "ins":
                nins -= 1; ndel += 1
        ncpg += CPG.get(look, 0); ncod += COD.get(look, 0)
        e = EL.get(look, "NA")
        ec[ELEMS.index(e) if e in ELEMS else 5] += 1
    rec.append((ncom, nrare, nun, nsnv, nins, ndel, ncpg, ncod, ec,
                ";".join(sorted(sig)), len(sig), (fmin if sig else np.nan)))
M = pd.DataFrame(dict(
    sample=L2["sample"], slot=L2.slot, array_pos=L2.array_pos, ap=L2.ap, ad=L2.ad, side=L2.side,
    N=L2.N, data_source=L2.data_source, pop=L2["pop"], super_pop=L2.super_pop,
    # L2.pop is not a column but the DataFrame.pop method, always use L2["pop"]
    both_eligible=L2.both_eligible,
    varlist=L2.varlist_maj,
    n_var=[len(x.split(",")) if x else 0 for x in newlist],
    n_var_common=[r[0] for r in rec], n_var_rare=[r[1] for r in rec],
    n_var_unmatched=[r[2] for r in rec], n_snv=[r[3] for r in rec],
    n_ins=[r[4] for r in rec], n_del=[r[5] for r in rec],
    n_var_cpg=[r[6] for r in rec], n_var_coding=[r[7] for r in rec],
    vsig=[r[9] for r in rec], vsig_n=[r[10] for r in rec], vsig_minfreq=[r[11] for r in rec]))
for i, e in enumerate(ELEMS):
    M["nv_" + e] = [r[8][i] for r in rec]
# Version that keeps every L2 column and swaps in only the variant columns, so code can run by changing only the path
FULL = L2.drop(columns=["varlist_maj"]).copy()
for c in M.columns:
    if c in ("sample", "slot", "array_pos", "ap", "ad", "side", "N",
             "data_source", "pop", "super_pop", "both_eligible"):
        continue
    FULL[c] = M[c].values
FULL.to_parquet(os.path.join(NEW, "L2_HAP_UNIT_MAJREF.parquet"), index=False)
M.to_parquet(os.path.join(NEW, "L2_MAJREF_variant_cols_only.parquet"), index=False)
P("   갈아끼운 변이 열 %d 개 · 전체 열 %d 개 유지"
  % (len([c for c in M.columns if c not in ("sample", "slot", "array_pos", "ap", "ad", "side",
                                            "N", "data_source", "pop", "super_pop", "both_eligible")]),
     FULL.shape[1]))

# ---------- 4. before/after comparison table ----------
BA = []
for col in ["n_var", "n_var_common", "n_var_rare", "n_var_unmatched",
            "n_snv", "n_ins", "n_del", "n_var_cpg", "n_var_coding", "vsig_n"] + ["nv_" + e for e in ELEMS]:
    if col not in L2.columns:
        continue
    BA.append(dict(열=col, 전_평균=L2[col].mean(), 후_평균=M[col].mean(),
                   전_중앙=L2[col].median(), 후_중앙=M[col].median(),
                   차=M[col].mean() - L2[col].mean()))
B = pd.DataFrame(BA)
B.to_csv(os.path.join(NEW, "BEFORE_AFTER_unit_variant_cols.tsv"), sep="\t", index=False, float_format="%.4f")
P()
P("[전/후 — 유닛 변이 열]")
P("   %-16s %10s %10s %10s" % ("열", "전 평균", "후 평균", "차"))
for _, r in B.iterrows():
    P("   %-16s %10.3f %10.3f %+10.3f" % (r.열, r.전_평균, r.후_평균, r.차))
# new token frequency table
TOKN = pd.DataFrame({"tok": list(cnt.keys()), "n_unit_seen": list(cnt.values())})
TOKN["unit_freq"] = TOKN.n_unit_seen / nu
TOKN["역토큰인가"] = TOKN.tok.isin(INV2FWD).map({True: "Y", False: ""})
TOKN.sort_values("n_unit_seen", ascending=False).to_csv(
    os.path.join(NEW, "VAR_GROUP_MAP_MAJREF.tsv"), sep="\t", index=False, float_format="%.6f")
P()
P("   절반 넘는 토큰 : 전 %d 개 → 후 %d 개"
  % (int((VG.unit_freq.fillna(0) > 0.5).sum()), int((TOKN.unit_freq > 0.5).sum())))
P()
P("[저장] %s" % NEW.replace("/", "\\"))
for f in ("L2_HAP_UNIT_MAJREF.parquet", "MAJREF_FLIPS_this_run.tsv",
          "BEFORE_AFTER_unit_variant_cols.tsv", "VAR_GROUP_MAP_MAJREF.tsv", "MAJREF_BUILD.txt"):
    P("   %s" % f)
LOG.close()
