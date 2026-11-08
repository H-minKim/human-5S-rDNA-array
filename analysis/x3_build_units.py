# -*- coding: utf-8 -*-
# 타 생물종 — 유닛 하나하나의 구조와 메틸을 한 표로
#
#   유닛 좌표 : 06_units_from_canonical/<array_id>.units.tsv
#               (array_pos · chrom_start/end · unit_len · strand)
#   메틸      : 05_bedmethyl/<species>/<시료>.5S.cpg.bed  (동물, CpG)
#               arabidopsis 는 allC.bed.gz 에서 CG 만
#   변이      : 07_consensus_from_reads/*.unit_consensus.tsv
#               (array_pos 별 n_var_cons · n_sub/ins/del · n_var_strict · depth)
#
#   '-' 가닥 array 는 array_pos 순서와 조립본 좌표 순서가 반대다 — 검산해서 확인한다.
# Output is written only to this folder.
import io, os, glob, gzip
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
LOG = io.open(OUT + "/X3_UNITS.txt", "w", encoding="utf8", newline="\n")


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


P("=" * 104)
P("유닛 구조 + 유닛별 메틸")
P("=" * 104)
P()

# ---------------------------------------------------------------- 유닛 좌표
fs = sorted(glob.glob(R + "/06_units_from_canonical/*.units.tsv"))
U = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
U = U[U.verify.astype(str) == "OK"] if "verify" in U.columns else U
P("  유닛 좌표 %s 행 · array %d · 종 %d"
  % (format(len(U), ","), U.array_id.nunique(), U.species_key.nunique()))

# 가닥 검산
P()
P("  '-' 가닥 array 에서 array_pos 와 좌표 방향이 반대인가")
chk = []
for aid, g in U.groupby("array_id"):
    if len(g) < 3:
        continue
    g = g.sort_values("array_pos")
    rho = np.corrcoef(g.array_pos, g.chrom_start)[0, 1]
    chk.append((aid, g.strand.iloc[0], rho))
C = pd.DataFrame(chk, columns=["array_id", "strand", "rho"])
for s, g in C.groupby("strand"):
    P("   가닥 %-3s : array %3d · rho(array_pos, chrom_start) 중앙 %+.3f"
      % (s, len(g), g.rho.median()))
P("   ⇒ 부호가 가닥과 맞으면 units.tsv 가 이미 array 순서를 맞춰 둔 것이다")
P()

# ---------------------------------------------------------------- 메틸 site
def load_cpg(path):
    """modkit bedMethyl -> (chrom, pos, nmod, ncan)  'm' 행만"""
    rows = []
    op = gzip.open if path.endswith(".gz") else io.open
    with op(path, "rt") as f:
        for l in f:
            t = l.rstrip("\n").split("\t")
            code = t[3]
            if "allC" in path:                          # allC 만 문맥 표기
                if not code.startswith("m,"):
                    continue
                if code.split(",")[1] != "CG":
                    continue
            else:
                if code != "m":
                    continue
            rows.append((t[0], int(t[1]), int(t[11]), int(t[12])))
    return pd.DataFrame(rows, columns=["chrom", "pos", "nmod", "ncan"])


SAMP = []
for sp in sorted(os.listdir(R + "/05_bedmethyl")):
    d = R + "/05_bedmethyl/" + sp
    if not os.path.isdir(d):
        continue
    for f in sorted(os.listdir(d)):
        if f.endswith(".bed") or f.endswith(".bed.gz"):
            SAMP.append((sp, f.split(".")[0], d + "/" + f))
P("  bedmethyl %d 개" % len(SAMP))

rows = []
for sp, sm, path in SAMP:
    B = load_cpg(path)
    Us = U[U.species_key == sp]
    if not len(Us) or not len(B):
        P("   %-24s %-14s  (유닛 또는 site 없음)" % (sp, sm))
        continue
    tot = 0
    for ch, gb in B.groupby("chrom"):
        gu = Us[Us.chrom == ch]
        if not len(gu):
            continue
        pos = np.sort(gb.pos.values)
        order = np.argsort(gb.pos.values)
        nm = gb.nmod.values[order]
        nc = gb.ncan.values[order]
        cm = np.concatenate([[0], np.cumsum(nm)])
        cc = np.concatenate([[0], np.cumsum(nc)])
        cn = np.arange(len(pos) + 1)
        s0 = np.searchsorted(pos, gu.chrom_start.values - 1, "left")
        s1 = np.searchsorted(pos, gu.chrom_end.values, "right")
        for k, r in enumerate(gu.itertuples(index=False)):
            a, b = s0[k], s1[k]
            if b <= a:
                continue
            M = cm[b] - cm[a]
            Cn = cc[b] - cc[a]
            rows.append((sp, sm, r.array_id, r.array_pos, r.unit_len, r.strand,
                         int(cn[b] - cn[a]), int(M), int(M + Cn)))
            tot += 1
    P("   %-24s %-14s  site %9s -> 유닛 %s"
      % (sp, sm, format(len(B), ","), format(tot, ",")))

M = pd.DataFrame(rows, columns=["species_key", "sample", "array_id", "array_pos",
                                "unit_len", "strand", "n_site", "n_mod", "n_call"])
M["meth"] = M.n_mod / M.n_call.replace(0, np.nan)
P()
P("  유닛-시료 %s 행 · array %d" % (format(len(M), ","), M.array_id.nunique()))

# ---------------------------------------------------------------- 변이
cs = sorted(glob.glob(R + "/07_consensus_from_reads/*.unit_consensus.tsv"))
V = pd.concat([pd.read_csv(f, sep="\t", low_memory=False) for f in cs],
              ignore_index=True)
keep = ["species_key", "sample", "array_id", "array_pos", "refunit_len",
        "depth_median", "n_var_cons", "n_sub", "n_ins", "n_del", "n_var_strict"]
keep = [c for c in keep if c in V.columns]
V = V[keep]
for c in ("depth_median", "n_var_cons", "n_sub", "n_ins", "n_del", "n_var_strict"):
    if c in V.columns:
        V[c] = pd.to_numeric(V[c], errors="coerce")
P("  유닛 변이 %s 행 · array %d" % (format(len(V), ","), V.array_id.nunique()))

#  시료 이름이 판마다 다르다 (bedmethyl "Control1" vs 합본 "Control1.g12")
V["sample_key"] = V["sample"].astype(str).str.split(".").str[0]
M["sample_key"] = M["sample"].astype(str).str.split(".").str[0]
V2 = V.drop(columns=["sample"])
M = M.merge(V2, on=["species_key", "sample_key", "array_id", "array_pos"],
            how="left")
NM = {"zebrafish_GRCz12tu": "제브라피시", "zfinch_v03_dip": "금화조",
      "arabidopsis_ColCC": "애기장대", "bonobo_mPanPan1": "보노보",
      "chimp_mPanTro3": "침팬지", "gorilla_mGorGor1": "고릴라",
      "orang_bornean_mPonPyg2": "보르네오오랑",
      "orang_sumatran_mPonAbe1": "수마트라오랑", "siamang_mSymSyn1": "시아망"}
M["종"] = M.species_key.map(NM)

# 끝에서 센 거리 (양쪽) — 가닥/방향에 안 휘는 축
g = M.groupby(["array_id"])["array_pos"].transform("max")
M["N_unit"] = g
M["p_from_1"] = M.array_pos
M["p_from_N"] = M.N_unit - M.array_pos + 1
M["p_end"] = np.minimum(M.p_from_1, M.p_from_N)
M["rel"] = (M.array_pos - 1) / np.maximum(M.N_unit - 1, 1)
M.to_csv(OUT + "/tables/OTHERSP_UNIT_METH.tsv", sep="\t", index=False,
         float_format="%.5f")
P()
P("=" * 104)
P("종별 유닛 구조")
P("=" * 104)
P()
P("   종            array  유닛수 범위   유닛길이 중앙(범위)        유닛당 변이 중앙   깊이 중앙")
for sp, g2 in M.groupby("종", sort=False):
    au = g2.groupby("array_id").array_pos.max()
    P("   %-12s %4d   %4d-%-5d   %6.0f (%.0f-%.0f)    %s    %s"
      % (sp, g2.array_id.nunique(), au.min(), au.max(),
         g2.unit_len.median(), g2.unit_len.min(), g2.unit_len.max(),
         ("%5.1f" % g2.n_var_cons.median()) if g2.n_var_cons.notna().any() else "  -  ",
         ("%5.0f" % g2.depth_median.median()) if g2.depth_median.notna().any() else "  -  "))
P()
P("  저장 : 10_other_species/tables/OTHERSP_UNIT_METH.tsv")
LOG.close()
print("done")
