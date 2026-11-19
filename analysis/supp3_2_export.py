# -*- coding: utf-8 -*-
# Write the Source Data table for supplementary panel 3-2.
#   The drawing script plots straight from L1_READ_UNIT with no intermediate table,
#          so there is no table for the Source Data workbook. This writes exactly what the figure shows.
#   The figure is not redrawn; only the table is written, under the same filter and the same ordering.
#   Check: 52 molecules, 1,200 covered bins and 83 positions, which must match the figure's own output.
import io, os, sys, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
# The data roots come from environment variables.
#   They have no default: the script stops with a message if one is not set.
DATA_ROOT = os.environ["DATA_ROOT"]

S = DATA_ROOT + "/"
# Where the tables are read and written.
# Write the deposited table's column names and values in English.
#   The renaming happens at a single point just before writing.
_KO2EN_COL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": {"구역": "zone", "유닛안자리": "pos_in_unit",
                                    "유전자기준": "pos_gene_ref", "호출수": "n_call",
                                    "메틸호출수": "n_meth_call"},
    "FIG3_B_POSBETA.tsv": {"층": "stratum", "축": "axis", "자리": "pos"},
    "FIG3_C_LOWCOUNT.tsv": {"컷": "cut", "층": "stratum", "N중앙": "N_median",
                            "저메틸수": "n_low_meth", "저메틸수_sd": "n_low_meth_sd",
                            "비례했다면": "expected_if_proportional",
                            "저메틸몫": "low_meth_fraction"},
    "SUPP3_2_POSSD.tsv": {"분자수": "n_molecules", "평균_b_call": "mean_b_call",
                          "분자사이_sd": "sd_between_molecules"},
    "SUPP3_B_HAP3.tsv": {"하플로타입": "haplotype", "길이층": "length_class",
                         "유닛수": "n_units", "분자순번": "molecule_index",
                         "유닛자리": "unit_position"},
}
_KO2EN_VAL = {
    "FIG3_A_UNITPROFILE_ZONE.tsv": ("구역", {"근위": "proximal", "내부": "interior",
                                            "원위": "distal", "전체": "all"}),
    "FIG3_B_POSBETA.tsv": ("층", {"전체": "all"}),
    "FIG3_C_LOWCOUNT.tsv": ("컷", {"0.508(여기서 잴 골)": "0.508 (trough measured in this work)"}),
    "SUPP3_B_HAP3.tsv": ("길이층", {"짧은 array": "short array", "중간 길이": "medium",
                                   "긴 array": "long array"}),
}


def _to_english(n, df):
    b = os.path.basename(n)
    v = _KO2EN_VAL.get(b)
    if v and v[0] in df.columns:
        df = df.copy()
        df[v[0]] = df[v[0]].map(lambda x: v[1].get(x, x))
    c = _KO2EN_COL.get(b)
    if c:
        df = df.rename(columns=c)
    return df

HERE = os.environ.get("TABLE_DIR", os.path.dirname(os.path.abspath(__file__)))
HAPS, HAPSLOT = "NA19338", 2

L = pd.read_parquet(S + "02_tables/L1_READ_UNIT.parquet",
                    columns=["sample", "slot", "read_id", "ap", "ad", "N", "gate", "b_call"])
L = L[(L.gate == 1) & (L["sample"].astype(str) == HAPS) & (L.slot.astype(int) == HAPSLOT)].copy()
NN = int(L.N.iloc[0])
order = list(L.groupby("read_id").b_call.mean().sort_values().index)
pos = {r: i for i, r in enumerate(order)}

# ---- A and B: molecule by unit, in long form ----
A = L[["read_id", "ap", "ad", "b_call"]].copy()
A["분자순번"] = A.read_id.map(pos)          # in increasing mean methylation, the vertical order of the figure
A["하플로타입"] = "%s|%d" % (HAPS, HAPSLOT)
A["array_길이_유닛"] = NN
A = A.sort_values(["분자순번", "ap"])[
    ["하플로타입", "array_길이_유닛", "read_id", "분자순번", "ap", "ad", "b_call"]]

# ---- C: the spread between molecules, per position ----
g = L.groupby("ap").b_call.agg(["mean", "std", "count"]).reset_index()
g.columns = ["ap", "평균_b_call", "분자사이_sd", "분자수"]
g["ad"] = NN - g.ap + 1
C = g[g.분자수 >= 5][["ap", "ad", "분자수", "평균_b_call", "분자사이_sd"]]

for df, n in ((A, "SUPP3_2_MOLECULES.tsv"), (C, "SUPP3_2_POSSD.tsv")):
    p = os.path.join(HERE, n)
    _to_english(n, df).to_csv(p, sep="\t", index=False, float_format="%.5f", lineterminator="\n")
    print("saved %-28s rows %6d, md5 %s"
          % (n, len(df), hashlib.md5(io.open(p, "rb").read()).hexdigest()[:12]))

print()
print("check: does it match the figure output")
print("   molecules (read) %d        (figure: 52)" % A.read_id.nunique())
print("   covered cells %d           (figure: 1200)" % len(A))
print("   array units %d        (figure: 83)" % NN)
print("   positions with read>=5 %d" % len(C))
for nm, m in (("근위 1-10", C.ap <= 10), ("가운데 21-55", (C.ap > 20) & (C.ap <= 55)),
              ("끝 71-83", C.ap > 70)):
    s = C[m]
    print("   %-12s mean β %.3f, mean between-molecule sd %.3f (max %.3f)"
          % (nm, s.평균_b_call.mean(), s.분자사이_sd.mean(), s.분자사이_sd.max()))
print("   (figure legend: 0.093 / 0.059 / 0.026)")
