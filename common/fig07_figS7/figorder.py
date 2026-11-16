# -*- coding: utf-8 -*-
# The panel order is by chromosome number, then by coordinate.
#   Within a species: increasing chromosome number, then increasing array start coordinate, counted per locus, with the maternal and paternal copies merged.
import re
import pandas as pd

_CACHE = {}


def panel_order(NEW):
    """{group name : (chromosome number, start coordinate)}; smaller values come first"""
    if NEW in _CACHE:
        return _CACHE[NEW]
    MT = pd.read_csv(NEW + "/tables/PERARRAY_META.tsv", sep="\t")
    CN = pd.read_csv(NEW + "/tables/CHROM_NAMES.tsv", sep="\t")
    G = pd.read_csv(NEW + "/tables/UNITBETA_PERARRAY_v13.tsv.gz", sep="\t",
                    usecols=["array_id", "묶음"]).drop_duplicates()
    chno = dict(zip(CN.접근번호.astype(str), CN.염색체번호.astype(str)))
    k = dict(zip(G.array_id, G.묶음))
    out = {}
    for t in MT.itertuples(index=False):
        key = k.get(t.array_id)
        if key is None:
            continue
        lab = str(t.이름표).replace("_mat", "").replace("_pat", "")
        nm = lab if lab.lower().startswith("chr") else chno.get(str(t.chrom).replace("_mat", "").replace("_pat", ""),
                                                                str(t.chrom))
        m = re.search(r"(\d+)", str(nm))
        n = int(m.group(1)) if m else 999
        s = float(t.array시작)
        if key not in out or (n, s) < out[key]:
            out[key] = (n, s)
    _CACHE[NEW] = out
    return out
