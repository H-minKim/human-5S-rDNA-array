# -*- coding: utf-8 -*-
# — array − flank 를 정수 n_5mC 로 다시 합산 (x1 은 반올림된 frac_5mC x n_c 였음)
#   원본 x1 표(OTHERSP_ARRAY_METH.tsv)는 그대로 두고 새 표를 낸다.
#   두 정의를 같이 낸다 :
# (a) per array, summed over samples (x1 definition; for zebra finch LIG and UL are summed) (b) per sample
#   입력 : $ANALYSIS_ROOT/2026-09-11_5S_other_species/data_v2/09_array_vs_flank/AVF.tsv (읽기만)
#   출력 : tables/OTHERSP_ARRAY_METH_n5mC.tsv · X16_AVF_N5MC.txt
import io
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
LOG = io.open(OUT + "/X16_AVF_N5MC.txt", "w", encoding="utf8", newline="\n")


def P(s=""):
    print(s)
    LOG.write(s + "\n")


A = pd.read_csv(R + "/09_array_vs_flank/AVF.tsv", sep="\t")
for c in ("n_c", "n_5mC", "frac_5mC", "reads"):
    A[c] = pd.to_numeric(A[c], errors="coerce")
OLD = pd.read_csv(OUT + "/tables/OTHERSP_ARRAY_METH.tsv", sep="\t")
OLD = OLD[OLD.src == "AVF"].set_index("array_id")


def pool(g, how):
    n = g.n_c.sum()
    if n <= 0:
        return np.nan
    return (g.n_5mC.sum() if how == "int" else (g.n_c * g.frac_5mC).sum()) / n


rows = []
for (sp, aid), g in A.groupby(["species_key", "array_id"]):
    ar, fl = g[g.region == "array"], g[g.region != "array"]
    r = dict(species_key=sp, array_id=aid,
             meth_frac=pool(ar, "frac"), flank_frac=pool(fl, "frac"),
             meth_int=pool(ar, "int"), flank_int=pool(fl, "int"),
             n_call=int(ar.n_c.sum()), n_call_flank=int(fl.n_c.sum()), n_sample=ar["sample"].nunique())
    r["delta_frac"] = r["meth_frac"] - r["flank_frac"]
    r["delta_int"] = r["meth_int"] - r["flank_int"]
    r["delta_x1_table"] = (OLD.loc[aid, "meth"] - OLD.loc[aid, "flank"]) if aid in OLD.index else np.nan
    # (나) 시료별
    per = []
    for s, gs in g.groupby("sample"):
        a2, f2 = gs[gs.region == "array"], gs[gs.region != "array"]
        if len(a2) and len(f2):
            per.append((s, pool(a2, "int") - pool(f2, "int"), int(a2.reads.sum()), int(f2.reads.sum())))
    r["delta_int_per_sample"] = "; ".join("%s %+.3f (array read %d · flank read %d)" % p for p in per)
    rows.append(r)
T = pd.DataFrame(rows)
T.to_csv(OUT + "/tables/OTHERSP_ARRAY_METH_n5mC.tsv", sep="\t", index=False, float_format="%.5f")

P("x1 표 재현 검산 : delta_frac vs x1 표 delta 최대 차 %.2e" % np.nanmax(np.abs(T.delta_frac - T.delta_x1_table)))
P("정수 n_5mC 로 바꾼 효과 : |delta_int − delta_frac| 최대 %.4f · 중앙 %.5f" %
  (np.nanmax(np.abs(T.delta_int - T.delta_frac)), np.nanmedian(np.abs(T.delta_int - T.delta_frac))))
P()
GRP = {"유인원 5 (시아망 제외)": ["bonobo_mPanPan1", "chimp_mPanTro3", "gorilla_mGorGor1", "orang_bornean_mPonPyg2", "orang_sumatran_mPonAbe1"],
       "시아망": ["siamang_mSymSyn1"]}
for nm, keys in GRP.items():
    g = T[T.species_key.isin(keys)]
    P("%-22s x1(frac) %+.3f ~ %+.3f  →  정수 n_5mC %+.3f ~ %+.3f" % (nm, g.delta_frac.min(), g.delta_frac.max(), g.delta_int.min(), g.delta_int.max()))
    for r in g.itertuples():
        P("   %-45s frac %+.4f · int %+.4f" % (r.array_id.split("__")[1], r.delta_frac, r.delta_int))
f = T[(T.species_key == "zfinch_v03_dip") & ~T.array_id.str.contains("CM109798|CM109759")]
P()
P("금화조 chr2·chr4 (6 array)")
P("  (가) array 마다 LIG+UL 합산 : frac %+.3f ~ %+.3f · int %+.3f ~ %+.3f" % (f.delta_frac.min(), f.delta_frac.max(), f.delta_int.min(), f.delta_int.max()))
allps = []
for r in f.itertuples():
    P("   %-38s 합산 int %+.3f | 시료별 %s" % (r.array_id.split("__")[1], r.delta_int, r.delta_int_per_sample))
    for part in r.delta_int_per_sample.split("; "):
        allps.append(float(part.split()[1]))
P("  (나) 시료별 : %+.3f ~ %+.3f (%d 행)" % (min(allps), max(allps), len(allps)))
LOG.close()
