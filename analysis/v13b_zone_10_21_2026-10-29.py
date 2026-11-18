# -*- coding: utf-8 -*-
# Methylation by position within the unit, accumulated per zone (proximal / interior / distal).
#
#   The zones are the common zones used by the other panels of this figure.
#
#   Zones, grouped from the common five; no new boundary is introduced:
#     proximal = ap <= ad and ap <= 17 (the common "proximal 1-3" plus "proximal 4-17")
#     distal = ad < ap and ad <= 22 (the common "distal 11-22" plus "distal 1-10")
#     interior = everything else (the common "interior"). The three partition the whole set.
#     (ad = array_pos, ap = N - array_pos + 1)
#
#   Ruler and gates, the same as the other panels:
#     hard call at ML >= 128; each read x unit needs at least 150 canonical CpG calls and coverage of at least 0.9.
#     To keep the cohort identical to the other panels of this figure, it is further restricted to the 45,544 units of L2_HAP_UNIT with n_read >= 5.
#     A bin's value is sum(methylated calls) / sum(calls), pooled over read x unit; it is not a mean over haplotypes.
#
#   To run: python v13_unitprofile_zone_2026-10-18.py (about 4 minutes, over the 6 raw files in order).
import io, os, re, gzip, glob, time, collections, hashlib
import numpy as np
import pandas as pd
# The data roots come from environment variables.
#   They have no default: the script stops with a message if one is not set.
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

DATA_ROOT = os.environ["DATA_ROOT"]
V_ROOT = os.environ["V_ROOT"]
CANON_FA = os.environ["CANON_FA"]

V = V_ROOT
OUT = DATA_ROOT
CANON = CANON_FA
MAP = os.path.join(V, r"03_maps/g8canonv10/CONSENSUS_MAP_G8_full.tsv")
RAW = os.path.join(V, r"08_methylation/raw")
DST = os.path.join(OUT, r"03_figures/paper_Fig/FIG3_A_UNITPROFILE_ZONE.tsv")
LOG = io.open(os.path.join(OUT, r"04_logs/V13_UNITPROFILE_ZONE.txt"), "w",
              encoding="utf8", newline="\n")
# The script reproduces FIG3_A_UNITPROFILE_ZONE.tsv byte for byte.
MLCUT = 128


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


# ---------- canonical CpG sites ----------
canon = "".join(x.strip() for x in io.open(CANON, encoding="utf8") if x[0] != ">")
CPG = sorted(i + 1 for i in range(len(canon) - 1) if canon[i:i + 2] == "CG")
CPGSET = set(CPG)

# ---------- the unit set this figure uses, the same cohort as the other panels ----------
L2 = pd.read_parquet(os.path.join(OUT, r"02_tables/L2_HAP_UNIT.parquet"),
                     columns=["sample", "slot", "ap", "ad", "n_read"])
L2 = L2[L2.n_read >= 5]
KEEP = set(zip(L2["sample"].astype(str), L2.slot.astype(str), L2.ad.astype(int)))
NHAP = L2.assign(h=L2["sample"].astype(str) + "|" + L2.slot.astype(str)).h.nunique()

# ---------- unit coordinates, shifted by the variant list ----------
IN = re.compile(r"^(\d+)I(\d+)$")
DE = re.compile(r"^(\d+)D(\d+)$")
UMAP, UN = {}, {}
with io.open(MAP, encoding="utf8") as h:
    hd = h.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hd)}
    for line in h:
        p = line.rstrip("\n").split("\t")
        if p[ix["seg_type"]] != "UNIT":
            continue
        smp = p[ix["sample"]]; sl = p[ix["slot"]]
        ap = int(p[ix["array_pos"]]); ul = int(p[ix["unit_len"]])
        UN[(smp, sl)] = int(p[ix["N"]])
        vl = p[ix["varlist"]]
        ins, dele = {}, {}
        if vl and vl != "-":
            for v in vl.split(","):
                m = IN.match(v)
                if m:
                    ins[int(m.group(1))] = int(m.group(2)); continue
                m = DE.match(v)
                if m:
                    dele[int(m.group(1))] = int(m.group(2))
        arr = np.zeros(ul + 2, dtype=np.int16)
        sp = 1; cp = 1; skip = 0
        while cp <= 2240 and sp <= ul:
            if skip > 0:
                skip -= 1; cp += 1; continue
            if cp in dele:
                skip = dele[cp] - 1; cp += 1; continue
            arr[sp] = cp; sp += 1
            if cp in ins:
                for _ in range(ins[cp]):
                    if sp <= ul:
                        arr[sp] = 0; sp += 1
            cp += 1
        UMAP[(smp, sl, ap)] = arr

P("=" * 104)
P("Methylation by position within the unit, by zone (Fig3 (A))")
P("=" * 104)
P()
P("  canonical CpG %d, MLCUT %d" % (len(CPG), MLCUT))
P("  zones: proximal (ap<=ad and ap<=10), interior (the rest), distal (ad<ap and ad<=21)")
P("  cohort: units with n_read>=5 in L2_HAP_UNIT %s, haplotypes %d" % (format(len(KEEP), ","), NHAP))
P()

ACC = collections.defaultdict(lambda: [0, 0])       # (zone, cpg) -> [n, n_meth]
NU = collections.Counter()                          # number of read x unit observations per zone
USEEN = collections.defaultdict(set)                # distinct units per zone
n_gate_out = 0
files = sorted(glob.glob(os.path.join(RAW, "METH_READUNIT_*.tsv.gz")))
for f in files:
    t1, nk = time.time(), 0
    with gzip.open(f, "rt", encoding="utf8") as h:
        h.readline()
        for line in h:
            p = line.rstrip("\n").split("\t")
            if p[1] != "HPRC_ONT_submissions":
                continue
            smp, sl = p[0], p[4]
            ad = int(p[5])                           # array_pos is the index from the distal end
            if (smp, sl, ad) not in KEEP:            # only the units used by the other panels
                continue
            cs = p[11]
            if not cs or cs == ".":
                continue
            arr = UMAP.get((smp, sl, ad))
            if arr is None:
                continue
            N = UN.get((smp, sl), 0)
            if N <= 0:
                continue
            ap = N - ad + 1                          # index from the proximal end
            L = len(arr) - 1
            hit, nraw = [], 0
            for it in cs.split(","):
                a, v = it.split(":")
                a = int(a); v = int(v)
                nraw += 1
                if a < 1 or a > L:
                    continue
                cp = int(arr[a])
                if cp == 0 or cp not in CPGSET:
                    continue
                hit.append((cp, 1 if v >= MLCUT else 0))
            if len(hit) < 150 or len(hit) / max(nraw, 1) < 0.9:
                n_gate_out += 1
                continue
            if ap <= ad and ap <= 10:
                zone = "근위"
            elif ad < ap and ad <= 21:
                zone = "원위"
            else:
                zone = "내부"
            NU[zone] += 1
            USEEN[zone].add((smp, sl, ad))
            for cp, m in hit:
                a_ = ACC[(zone, cp)]
                a_[0] += 1; a_[1] += m
                a2 = ACC[("전체", cp)]
                a2[0] += 1; a2[1] += m
            nk += 1
    P("  %-28s read×unit %9s, %.0f s"
      % (os.path.basename(f)[14:-7], format(nk, ","), time.time() - t1))

rows = [(z, cp, cp - 2121, n, m, m / n) for (z, cp), (n, m) in ACC.items()]
D = _to_english("FIG3_A_UNITPROFILE_ZONE.tsv",
                pd.DataFrame(rows, columns=["구역", "유닛안자리", "유전자기준",
                                            "호출수", "메틸호출수", "beta"]))
D = D.sort_values(["zone", "pos_in_unit"])
D.to_csv(DST, sep="\t", index=False, float_format="%.5f")

P()
P("  rows %d (zones 4 × positions %d), saved %s" % (len(D), D.pos_in_unit.nunique(), os.path.basename(DST)))
P("  md5 %s" % hashlib.md5(io.open(DST, "rb").read()).hexdigest()[:12])
P()
P("  read×unit count and unit count by zone")
tot_u = 0
for z in ("근위", "내부", "원위"):
    P("    %-4s read×unit %9s, units %7s" % (z, format(NU[z], ","), format(len(USEEN[z]), ",")))
    tot_u += len(USEEN[z])
P("    %-4s %s (= %s + %s + %s, must equal 45,544)"
  % ("sum", format(tot_u, ","), *[format(len(USEEN[z]), ",") for z in ("근위", "내부", "원위")]))
P("  read×unit dropped at the gate %s" % format(n_gate_out, ","))
P()
P("  whole-unit mean β by zone (all positions pooled)")
for z in ("all", "proximal", "interior", "distal"):
    s = D[D.zone == z]
    P("    %-4s β %.4f (calls %s)" % (z, s.n_meth_call.sum() / s.n_call.sum(), format(int(s.n_call.sum()), ",")))
