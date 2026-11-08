# -*- coding: utf-8 -*-
# 타 생물종 5S array — array 하나하나의 메틸을 한 표로
#
#   자료 뿌리(읽기만) : $ANALYSIS_ROOT\2026-09-11_5S_other_species\data_v2
# Output is written only to this folder.
#
#   출처
#     제브라피시 31 : 00_catalog/ZEBRAFISH_ALL_ARRAYS.tsv  (array·flank·차 이미 있음)
#     그 외        : 09_array_vs_flank/AVF.tsv            (array / flank_left / flank_right)
#     애기장대 CG   : 05_bedmethyl/.../allC.bed.gz 에서 CG 문맥만 다시 셈
#                    (AVF 의 애기장대 값은 전문맥이라 다른 종과 못 비교한다)
#     생쥐          : 10_mouse_5S/README_MOUSE.md 의 내부대조 수치
#     붉은털(Tier2) : 11_rhesus_tier2/rhesus_tier2_5S_meth.tsv (공개 자리별 메틸, read 없음)
#
#  matplotlib 안 씀
import io, os, gzip, collections
import numpy as np
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


R = _need("RDNA_DATA_ROOT")
OUT = _need("RDNA_PREV_OUT_DIR")
LOG = io.open(OUT + "/X1_BUILD.txt", "w", encoding="utf8", newline="\n")


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


P("=" * 104)
P("타 생물종 5S array 메틸 — array 별 종합표")
P("=" * 104)
P()

# Species metadata
META = {
    "zebrafish_GRCz12tu":      ("제브라피시", "Danio rerio", "성체 뇌", "R10.4.1 guppy6.4.2", 6, 1),
    "zfinch_v03_dip":          ("금화조", "Taeniopygia guttata", "혈액/조직", "R10.4.1 dorado2.1.2", 1, 1),
    "arabidopsis_ColCC":       ("애기장대", "Arabidopsis thaliana", "유식물", "R10.4.1 dorado2.1.2", 1, 1),
    "bonobo_mPanPan1":         ("보노보", "Pan paniscus", "세포주", "R9.4.1 guppy6.3.7+6.3.8", 1, 1),
    "chimp_mPanTro3":          ("침팬지", "Pan troglodytes", "세포주", "R9.4.1 guppy6.3.8", 1, 1),
    "gorilla_mGorGor1":        ("고릴라", "Gorilla gorilla", "세포주", "R9.4.1 guppy6.3.8", 1, 1),
    "orang_bornean_mPonPyg2":  ("보르네오오랑", "Pongo pygmaeus", "세포주", "R9.4.1 guppy6.3.8", 1, 1),
    "orang_sumatran_mPonAbe1": ("수마트라오랑", "Pongo abelii", "세포주", "R9.4.1 guppy6.3.8", 1, 1),
    "siamang_mSymSyn1":        ("시아망", "Symphalangus syndactylus", "세포주", "R9.4.1 dorado0.9.1", 1, 1),
    "mouse_mm39":              ("생쥐", "Mus musculus", "뇌", "R10.4.1 dorado v5.2r1", 1, 1),
    "rhesus_T2T_MMU8v2":       ("붉은털원숭이", "Macaca mulatta", "공개 자리별", "(read 없음)", 0, 2),
}

rows = []

# ---------------------------------------------------------------- AVF
A = pd.read_csv(R + "/09_array_vs_flank/AVF.tsv", sep="\t")
A["frac_5mC"] = pd.to_numeric(A.frac_5mC, errors="coerce")
P("  AVF.tsv %d 행 · 종 %d · array %d"
  % (len(A), A.species_key.nunique(), A.array_id.nunique()))


def pooled(g):
    """호출 수로 가중한 pooled 메틸 (단순 평균 금지)"""
    n = g.n_c.sum()
    return (g.n_c * g.frac_5mC).sum() / n if n > 0 else np.nan, int(n), int(g.reads.sum())


for (sp, aid), g in A.groupby(["species_key", "array_id"]):
    ar = g[g.region == "array"]
    fl = g[g.region != "array"]
    if not len(ar):
        continue
    m_a, n_a, r_a = pooled(ar)
    m_f, n_f, r_f = pooled(fl) if len(fl) else (np.nan, 0, 0)
    rows.append(dict(species_key=sp, array_id=aid, chrom=ar.chrom.iloc[0],
                     array_start=int(ar.array_start.iloc[0]),
                     array_len=int(ar.array_len.iloc[0]),
                     meth=m_a, flank=m_f, n_call=n_a, n_call_flank=n_f,
                     reads=r_a, n_sample=ar["sample"].nunique(), src="AVF"))
T = pd.DataFrame(rows)
P("  -> array %d 개 (종 %d)" % (len(T), T.species_key.nunique()))
P()

# ---------------------------------------------------------------- 애기장대 CG 재계산
P("  애기장대 : AVF 값은 전문맥이라 CG 만 다시 센다")
bp = R + "/05_bedmethyl/arabidopsis_ColCC/ColCC_ATRun3.5S.allC.bed.gz"
ara = T[T.species_key == "arabidopsis_ColCC"].copy()
iv = [(r.chrom, r.array_start, r.array_start + r.array_len, r.array_id)
      for r in ara.itertuples(index=False)]
acc = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
with gzip.open(bp, "rt") as f:
    for l in f:
        t = l.rstrip("\n").split("\t")
        code = t[3]
        if not code.startswith("m,"):
            continue                        # 5mC 만 (h 는 따로)
        ctx = code.split(",")[1]
        c, s = t[0], int(t[1])
        nmod, ncan = int(t[11]), int(t[12])
        for cc, s0, s1, aid in iv:
            if c == cc and s0 <= s < s1:
                acc[aid][ctx][0] += nmod
                acc[aid][ctx][1] += nmod + ncan
                break
P("   array                                   CG        CHG       CHH      CG 호출수")
for aid in ara.array_id:
    v = acc[aid]
    def fr(k):
        a, b = v[k]
        return a / b if b else np.nan
    P("   %-38s %.4f    %.4f    %.4f   %s"
      % (aid.replace("arabidopsis_ColCC__", ""), fr("CG"), fr("CHG"), fr("CHH"),
         format(v["CG"][1], ",")))
    T.loc[T.array_id == aid, "meth_allctx"] = T.loc[T.array_id == aid, "meth"]
    T.loc[T.array_id == aid, "meth"] = fr("CG")
    T.loc[T.array_id == aid, "CHG"] = fr("CHG")
    T.loc[T.array_id == aid, "CHH"] = fr("CHH")
    T.loc[T.array_id == aid, "flank"] = np.nan          # flank 는 bed 범위 밖
P("   flank 는 bedmethyl 범위 밖이라 애기장대 Δ 는 못 낸다 (전문맥 Δ 만 있음)")
P()

# ---------------------------------------------------------------- 제브라피시 유닛수
Z = pd.read_csv(R + "/00_catalog/ZEBRAFISH_ALL_ARRAYS.tsv", sep="\t")
Z = Z[Z.meth_mean.astype(str) != "-"].copy()
zmap = {r.array_id: (r.refunit_len, r.g8_units, r.locus_type_lit, r.chrom_name)
        for r in Z.itertuples(index=False)}
P("  제브라피시 유닛 정보 %d array" % len(zmap))

# ---------------------------------------------------------------- 유닛 정보 일반
U = pd.read_csv(R + "/00_catalog/UNIT_REF_TABLE.tsv", sep="\t")
P("  UNIT_REF_TABLE %d 행 · 열 %s" % (len(U), list(U.columns)[:8]))
ucol = [c for c in U.columns if "unit" in c.lower() or "len" in c.lower()]
umap = {}
for r in U.itertuples(index=False):
    d = r._asdict()
    umap[d.get("array_id")] = d

for i, r in T.iterrows():
    aid = r.array_id
    if aid in zmap:
        T.loc[i, "unit_len"] = zmap[aid][0]
        T.loc[i, "n_unit"] = zmap[aid][1]
        T.loc[i, "note"] = zmap[aid][2]
        T.loc[i, "chrom_name"] = zmap[aid][3]
    elif aid in umap:
        d = umap[aid]
        for k in d:
            if k.lower() in ("refunit_len", "unit_len", "ref_len"):
                T.loc[i, "unit_len"] = d[k]
            if k.lower() in ("n_units", "g8_units", "units"):
                T.loc[i, "n_unit"] = d[k]

# ---------------------------------------------------------------- 생쥐 (README 실측)
T = pd.concat([T, pd.DataFrame([dict(
    species_key="mouse_mm39", array_id="mouse_mm39__chr8_124264035_124351175",
    chrom="chr8", array_start=124264035, array_len=87140,
    meth=0.114, flank=(0.585 * 98830 + 0.536 * 32635 + 0.519 * 38618 + 0.550 * 82800)
                 / (98830 + 32635 + 38618 + 82800),
    n_call=67911, n_call_flank=252883, reads=127, n_sample=1, src="README_MOUSE",
    unit_len=1701, n_unit=52, note="mm39 축약판(T2T 는 134kb/80copy)")])],
    ignore_index=True)

# ---------------------------------------------------------------- 붉은털 (Tier 2)
RH = pd.read_csv(R + "/11_rhesus_tier2/rhesus_tier2_5S_meth.tsv", sep="\t")
for c in ("mean_pct", "coverage_pct", "span_bp", "refunit_len", "g8_units",
          "cpg_sites_with_value"):
    if c in RH.columns:
        RH[c] = pd.to_numeric(RH[c], errors="coerce")
RH = RH.dropna(subset=["mean_pct"])
#  RH 파일에는 array 행과 flank 행이 섞여 있다. 갈라서 flank 는 대조로만 쓴다.
RH["is_flank"] = RH.array_status.astype(str).str.contains("FLANK")
fl_mean = ((RH[RH.is_flank].mean_pct * RH[RH.is_flank].cpg_sites_with_value).sum()
           / RH[RH.is_flank].cpg_sites_with_value.sum()) / 100.0
P("  붉은털 Tier2 : array %d 행 · flank %d 행 · flank 가중평균 %.3f"
  % ((~RH.is_flank).sum(), RH.is_flank.sum(), fl_mean))
for r in RH[~RH.is_flank].itertuples(index=False):
    if "NON_5S" in str(r.array_status):
        continue
    T = pd.concat([T, pd.DataFrame([dict(
        species_key="rhesus_T2T_MMU8v2", array_id=r.array_id, chrom="chr1",
        array_start=int(str(r.region).split(":")[1].split("-")[0]),
        array_len=int(r.span_bp), meth=r.mean_pct / 100.0, flank=fl_mean,
        n_call=int(r.cpg_sites_with_value), n_call_flank=0, reads=0, n_sample=0,
        src="Tier2 공개", unit_len=r.refunit_len, n_unit=r.g8_units,
        note="자리 피복 %.1f%% · %s" % (r.coverage_pct, r.array_status))])],
        ignore_index=True)

# ---------------------------------------------------------------- 메타 붙이기
T["종"] = T.species_key.map(lambda k: META.get(k, ("?",))[0])
T["학명"] = T.species_key.map(lambda k: META.get(k, ("", "?"))[1])
T["시료"] = T.species_key.map(lambda k: META.get(k, ("", "", "?"))[2])
T["호출기"] = T.species_key.map(lambda k: META.get(k, ("", "", "", "?"))[3])
T["tier"] = T.species_key.map(lambda k: META.get(k, ("", "", "", "", 0, 1))[5])
T["delta"] = T.meth - T.flank
T = T.sort_values(["species_key", "meth"])
T.to_csv(OUT + "/tables/OTHERSP_ARRAY_METH.tsv", sep="\t", index=False,
         float_format="%.5f")
P()
P("=" * 104)
P("종합표 — array %d 개 · 종 %d" % (len(T), T.species_key.nunique()))
P("=" * 104)
P()
P("   종            array  메틸 범위          flank 범위        Δ 범위           유닛수 범위")
for sp, g in T.groupby("species_key", sort=False):
    nm = META.get(sp, ("?",))[0]
    dd = g.delta.dropna()
    uu = g.n_unit.dropna()
    P("   %-12s %4d  %.3f-%.3f     %s     %s     %s"
      % (nm, len(g), g.meth.min(), g.meth.max(),
         ("%.3f-%.3f" % (g.flank.min(), g.flank.max())) if g.flank.notna().any() else "  (없음)   ",
         ("%+.3f~%+.3f" % (dd.min(), dd.max())) if len(dd) else "  (없음)   ",
         ("%d-%d" % (uu.min(), uu.max())) if len(uu) else "(없음)"))
P()
P("  저장 : 10_other_species/tables/OTHERSP_ARRAY_METH.tsv")
LOG.close()
print("done")
