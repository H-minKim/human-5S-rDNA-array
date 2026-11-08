# -*- coding: utf-8 -*-
"""
v11_io - single loading gate for version **v11** (the v11 edition of g8_io.py)
Reference files (every md5 checked against measurement):
    map        03_<map>\\g8canonv10\\CONSENSUS_MAP_G8_full.tsv   91816b9ca9d3426272d956fd11fbcc90
    summary    03_<map>\\g8canonv10\\CONSENSUS_SUM_G8_full.tsv    096fa12e3cf61aa63ab6538401fc8d27
    canon list 02_<canon>\\hap_seq_node_v8.tsv                     bc11cdc3aec7d0245cce26f899f5725f
    merged     07_<merged>\\READMAP_V11_BYSOURCE.tsv                 bdf4504e57d18f6695e8e01714a39351
    methyl     08_<methyl>\\out\\READMAP_METH_V11.tsv                40143744a85c92a2217875ecb1d17666
The earlier release (*_G8V4* in DATA_LATEST) is not used: the map and merged table frames were misaligned.

The canonical unit sequence is **identical** to the old version (measured).
  canon_ctx.fa = [canonical unit last 200][canonical unit 2240][canonical unit first 200], so the middle is cut out.
  Therefore the segments (AluY 1-318, NTS1 319-1178, GArich 1179-1551, NTS2 1552-2121, coding 2122-2240),
  U84C=2205, and 2235=mature 5S position 114 remain valid as is.

Polarity: `array_pos 1` is the genomic **distal (DJ)** end.
  Genomic proximal index = N - array_pos + 1. Checked automatically with U84C and printed at every load.

Unit record contract (same as g8_io - so existing analysis code can be used unchanged)
    r[0] ap  = unit number from the genomic proximal (PJ) junction = N - array_pos + 1
    r[1] ad  = unit number from the genomic distal (DJ) junction = array_pos
    r[2] N, r[3] beta (nan if absent), r[4] variant set (v4 notation), r[5] n_obs(-1), r[6] meta

Analysis set filters
    germline : exclude cancer cohort (COLO829-BL, COLO829-T, HG008-T)
    stage    : exclude B_ (canonical invalid) and E_ (canonical mismatch) - the stage table is STAGE_ALL_V11.tsv,
               recomputed on v11 by `stage_v11.py` (rules unchanged from stage_all4.py). The union edition is primary, stage_rep also exists.
    no_break : exclude slots with n_break>0 in SUM (15 slots in v11)
    hg002_2 : excluded (canonical changed to the T2T paternal, but 0 read units)
    verify   : OK only (exclude the 13 v11 NO_FULL)
    n_unit   : exclude units with alt=N
"""
import io, os, re, hashlib, collections
import numpy as np

# Version v11-fix2 (2026-09-17) - canonical and G8 map are fix1, ARRAY_ONLY units, merged table and methylation are fix2 (section 5 of each version definition .md)
V = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11-fix2")            # merged table, methylation, read x unit methylation raw
V1 = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11-fix1")           # canonical, G8 map
V0 = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11")                 # canonical file list (hap_seq_node_v8 - file name unchanged in fix1 and fix2)
MAPF = os.path.join(V1, "03_maps", "g8canonv11f1", "CONSENSUS_MAP_G8_full.tsv")
MAP_MD5 = "bb68fe3a2a23945fe919efe8fd405a7b"
SUMF = os.path.join(V1, "03_maps", "g8canonv11f1", "CONSENSUS_SUM_G8_full.tsv")
SUM_MD5 = "4c8252462824bdfe70d734d95254f55a"
READMAPF = os.path.join(V, "07_merged", "READMAP_V11F2_BYSOURCE.tsv")
READMAP_MD5 = "a3f379c91b00f81b54c43a5158794b11"
METHF = os.path.join(V, "08_methylation", "out", "READMAP_METH_V11F2.tsv")
METH_MD5 = "eb8d77649a8f54baccdcfe8e3dbc22f1"
# methylation merged table moved to the fix4 edition (default F4, `V11_METH=F2` reverts)
#   Measured: same rows and columns (55,599, 12), symmetric difference of the 3-column key sets 0, keys with different values 237
#          all 237 lie outside the h1 cohort - sites with different values inside h1: 0 (largest Δ is HG008-T = cancer)
# Two-stage gate: h1 chain tables 269 / 269 unchanged (F2 and F4 give the same result)
#     ⇒ no risk in moving to the newer canonical. Cancer and trio presets may change, so state F2 explicitly for those.
V4 = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\for_paper\08_versions\v11-fix4")
METHF4 = os.path.join(V4, "08_methylation", "out", "READMAP_METH_V11F4.tsv")
METH4_MD5 = "fc7042fc9ef3b1ff03729bcece4aff12"
METHF2, METH2_MD5 = METHF, METH_MD5
if os.environ.get("V11_METH", "").upper() != "F2":
    METHF, METH_MD5 = METHF4, METH4_MD5
METHSITEF = os.path.join(V, "08_methylation", "out", "METH_SITE_V11F2.tsv.gz")
CTXF = os.path.join(V1, "03_maps", "g8canonv11f1", "canon_ctx.fa")
HSNF = os.path.join(V0, "02_canonical", "hap_seq_node_v8.tsv")
CANDIR = os.path.join(V1, "02_canonical", "cmp4d_v11_fix1")
STAGEF = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\5S_position_variant_work\2026-09-17_v11_fix2\STAGE_ALL_V11F2.tsv")
OUT = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\5S_position_variant_work\2026-09-17_v11_fix2")

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

# ---------- canonical unit ----------
_ctx = "".join(l.strip() for l in io.open(CTXF) if not l.startswith(">")).upper()
assert len(_ctx) == 2640, len(_ctx)
CAN = _ctx[200:2440]
assert len(CAN) == 2240
# Sequence-level check (on the sequence itself, not the file md5)
assert CAN[2121:2141] == "GTCTACGGCCATACCACCCT", "유전자부 시작이 다르다"
assert CAN[:20] == "TTTTCTTTGGCTTTTTGCTG", "유닛 앞머리가 다르다"
assert CAN[2204] == "T" and CAN[2234] == "T", "U84C/2235 자리 염기가 다르다"

# Repeat blocks - rule from poscat.py lines 17-29 (same as the old version)
_runs = []
for m in re.finditer(r"(.)\1{3,}", CAN): _runs.append((m.start() + 1, m.end()))
for k in (2, 3, 4, 5, 6):
    for m in re.finditer(r"((.{%d}))\1{2,}" % k, CAN): _runs.append((m.start() + 1, m.end()))
_runs.sort(); _mg = []
for a, b in _runs:
    if _mg and a <= _mg[-1][1] + 1: _mg[-1][1] = max(_mg[-1][1], b)
    else: _mg.append([a, b])
REPB = {}
for i, (a, b) in enumerate(_mg, 1):
    for p in range(a, b + 1): REPB[p] = i
BLOCKS = [tuple(x) for x in _mg]

G8TOK = re.compile(r"^(\d+)(?:([ACGTN])|D(\d+)|I(\d+))$")
def to_v4(t):
    """G8 token -> old v4 token (repeat-block indels are merged by block number). None if alt=N."""
    m = G8TOK.match(t)
    if not m: raise ValueError("토큰 형식 아님: %r" % t)
    p = int(m.group(1))
    if m.group(2):
        x = m.group(2)
        if x == "N": return None
        ref = CAN[p - 1]
        if ref == x: raise ValueError("alt==ref: %r" % t)
        return ("R%d%s%s" if p in REPB else "%d%s%s") % (p, ref, x)
    if m.group(3):
        k = int(m.group(3)); st = p; en = p + k - 1
        bid = REPB.get(st) or REPB.get(en)
        return ("DB%d:%d" % (bid, k)) if bid else ("D%d:%d" % (st, k))
    k = int(m.group(4)); qp = p - 1
    bid = REPB.get(qp) or REPB.get(qp + 1)
    return ("IB%d:%d" % (bid, k)) if bid else ("I%d:%d:?" % (qp, k))

def _rows(p):
    hd = None
    for L in io.open(p, encoding="utf8"):
        f = L.rstrip("\n").split("\t")
        if hd is None: hd = f; H = {c: i for i, c in enumerate(hd)}; continue
        yield f, H

def slot_tables():
    SUM = {}; STAGE = {}; COH = {}
    for f, H in _rows(SUMF): SUM[(f[0], f[1])] = f
    for f, H in _rows(STAGEF):
        STAGE[(f[0], f[1])] = (f[H["stage"]], f[H["stage_rep"]])
        COH[f[0]] = f[H["cohort"]]
    return SUM, STAGE, COH, H

PRESETS = {
    "primary": dict(germline=True, drop_stage=("B_", "E_"), no_break=True, drop_hg002_2=True),
    "ok_only": dict(germline=True, keep_stage=("OK",), no_break=True, drop_hg002_2=True),
    "all":     dict(germline=False, drop_stage=(), no_break=False, drop_hg002_2=False),
    # added - tiers of the common cohort rule (control ). Only intersects primary with the tier list.
    #   primary 564 ∩ H1 = 557 (= all of H1), primary 564 ∩ H3 = 490 (= all of H3, the same set as H2 within this base)
    #   The three existing presets have no "tier" key, so the tier branch of keep_slot never runs (primary output unchanged).
    "h1":      dict(germline=True, drop_stage=("B_", "E_"), no_break=True, drop_hg002_2=True, tier="H1"),
    "h3":      dict(germline=True, drop_stage=("B_", "E_"), no_break=True, drop_hg002_2=True, tier="H3"),
}

TIERF = (__import__("os").environ["RDNA_DATA_ROOT"] + r"\rDNA_registry\registry\TIERS_H0H3_SLOTS_2026-10-05.tsv")
TIER_MD5 = "fc34c6326111d7ad797e34ed0a954b4c"
_TIER_CACHE = {}

def tier_slots(tier, check_md5=True):
    """Return the (sample, slot) set from the control tier list. The original list is read only."""
    if tier in _TIER_CACHE: return _TIER_CACHE[tier]
    if check_md5:
        got = md5(TIERF)
        if got != TIER_MD5: raise RuntimeError("층 명단 md5 불일치 %s (기대 %s)" % (got, TIER_MD5))
    hd = None; out = set()
    for L in io.open(TIERF, encoding="utf8"):
        f = L.rstrip("\n").split("\t")
        if hd is None: hd = {c: i for i, c in enumerate(f)}; continue
        if f[hd["tier"]] == tier: out.add((f[hd["sample"]], f[hd["slot"]]))
    if not out: raise RuntimeError("층 %s 가 명단에 없다" % tier)
    _TIER_CACHE[tier] = out
    return out

def load(preset="primary", check_md5=True, verbose=True, **kw):
    """-> (UNITS, ASM, META)   UNITS[(sample,slot)] = [record, ...] in ascending ap"""
    opt = dict(PRESETS[preset]); opt.update(kw)
    if check_md5:
        for p, h in ((MAPF, MAP_MD5), (SUMF, SUM_MD5)):
            got = md5(p)
            if got != h: raise RuntimeError("md5 불일치 %s : %s (기대 %s)" % (p, got, h))
    SUM, STAGE, COH, SH = slot_tables()
    sc = opt.get("stage_col", "stage")
    si = 0 if sc == "stage" else 1
    TS = tier_slots(opt["tier"], check_md5=check_md5) if opt.get("tier") else None
    def keep_slot(s, sl):
        if TS is not None and (s, sl) not in TS: return "tier:" + opt["tier"]
        if opt.get("germline") and COH.get(s, "").startswith("암"): return "cancer"
        st = STAGE.get((s, sl), ("?", "?"))[si]
        if opt.get("keep_stage") and not st.startswith(tuple(opt["keep_stage"])): return "stage"
        if opt.get("drop_stage") and st.startswith(tuple(opt["drop_stage"])): return "stage"
        row = SUM.get((s, sl))
        if opt.get("no_break") and row is not None and int(row[7]) > 0: return "break"
        if opt.get("drop_hg002_2") and (s, sl) == ("hg002", "2"): return "hg002_2"
        return None
    U = collections.defaultdict(list)
    drop = collections.Counter(); dslot = collections.Counter(); seen = set()
    for f, H in _rows(MAPF):
        if f[H["seg_type"]] != "UNIT": continue
        s, sl = f[H["sample"]], f[H["slot"]]
        why = keep_slot(s, sl)
        if why:
            if (s, sl) not in seen: dslot[why] += 1; seen.add((s, sl))
            drop["slot:" + why] += 1; continue
        v = f[H["verify"]]
        if v != "OK": drop["verify:" + v] += 1; continue
        raw = [] if f[H["varlist"]] == "-" else f[H["varlist"]].split(",")
        if any(t.endswith("N") for t in raw): drop["alt=N 유닛"] += 1; continue
        vs = set(x for x in (to_v4(t) for t in raw) if x is not None)
        N = int(f[H["N"]]); apos = int(f[H["array_pos"]])
        pol = opt.get("polarity", "genomic")
        if pol == "genomic":  gp, gd = N - apos + 1, apos
        elif pol == "g8":     gp, gd = apos, N - apos + 1
        else: raise ValueError("polarity 는 genomic 또는 g8")
        stg = STAGE.get((s, sl), ("?", "?"))
        U[(s, sl)].append((gp, gd, N, np.nan, vs, -1,
                           dict(verify=v, raw=raw, stage=stg[0], stage_rep=stg[1],
                                cohort=COH.get(s, "?"), array_pos=apos)))
    for k in U: U[k].sort(key=lambda r: r[0])
    ASM = {k: (v[0][2], "v11") for k, v in U.items()}
    n = sum(len(v) for v in U.values())
    a = [("2205TC" in r[4]) for rs in U.values() for r in rs if r[0] <= 3]
    b = [("2205TC" in r[4]) for rs in U.values() for r in rs if r[1] <= 3]
    fa = float(np.mean(a)) if a else float("nan"); fb = float(np.mean(b)) if b else float("nan")
    META = dict(preset=preset, opt=opt, polarity=opt.get("polarity", "genomic"), map_version="v11",
                stage_col=sc, n_samples=len(set(k[0] for k in U)), n_slots=len(U), n_units=n,
                dropped_units=dict(drop), dropped_slots=dict(dslot),
                n_blocks=len(BLOCKS), block_bp=sum(b_ - a_ + 1 for a_, b_ in BLOCKS),
                polarity_check=dict(u84c_r0_le3=fa, u84c_r1_le3=fb))
    if verbose:
        print("[v11_io %s · 극성 %s · 단계열 %s] 개체 %d · 슬롯 %d · UNIT %d   (반복블록 %d 개 %d bp)" % (
            preset, META["polarity"], sc, META["n_samples"], META["n_slots"], n,
            META["n_blocks"], META["block_bp"]))
        print("   제외 슬롯:", dict(dslot), " · 제외 유닛:", dict(drop))
        print("   극성 점검: U84C 보유율 r[0](근위)<=3 %.3f · r[1](원위)<=3 %.3f  %s" % (
            fa, fb, ("OK — 게놈 극성과 일치" if (META["polarity"] == "genomic" and fa > fb)
                     else ("G8 이름표 그대로" if META["polarity"] == "g8" else "경고: 게놈 극성과 반대"))))
    return dict(U), ASM, META

# ---------- new: methylation ----------
DEDUP_DEFAULT = os.environ.get("V11_METH_DEDUP", "pool").lower()

def load_meth(check_md5=True, min_reads=1, verbose=True, dedup=None):
    """READMAP_METH_V11 -> M[(sample,slot,array_pos)] = dict(beta_call, beta_read_mean,
       beta_read_median, n_reads, n_calls, n_meth, data_source)
       The site key is array_pos, the same as the map (genomic distal based). Proximal index is N-array_pos+1.

       2026-10-06 - some sites have rows from two data sources (measured)
         The merged table `*_BYSOURCE` has two data sources for 19 (individual, slot) pairs.
         So keying on the 3 columns `(sample,slot,array_pos)` makes `M[k] = …` overwrite silently,
         and which row wins is set by file line order (not by a rule).
         Measured: 1,115 overlapping keys (all with different values), 776 sites inside the h1 cohort, 13 slots
                the source pair is always `1kGP_ONT + HPRC_ONT_submissions`
                largest |Δbeta_call| between first and last row 0.7903

       `dedup` picks the resolution rule (default = `pool`, also settable by env var `V11_METH_DEDUP`)
         "pool" default : β = Σn_meth / Σn_calls - hard call pooling the reads of both sources.
                 It uses the same scale as the main table (= hard call), so only this keeps the scale consistent.
                 `n_reads`, `n_calls` and `n_meth` are summed too, and `data_source` is joined with "+".
         "deep"        : the one source with larger `n_reads`
         "last"        : for reproducing the old incident - last row in the file (do not use in new analyses)
       Measured effect (h1, 49,886 sites with methylation): overlap is 1.56%
         overall mean β last 0.662407, deep 0.662443, pool 0.662462 → shift 0.000055
    """
    how = (dedup or DEDUP_DEFAULT).lower()
    if how not in ("pool", "deep", "last"):
        raise ValueError("dedup 은 pool · deep · last 중 하나")
    if check_md5:
        got = md5(METHF)
        if got != METH_MD5: raise RuntimeError("메틸 합본 md5 불일치 : %s" % got)
    G4 = collections.defaultdict(list)          # collect by the 4-column key first (measured to have 0 overlap)
    for f, H in _rows(METHF):
        n = int(f[H["n_reads"]])
        if n < min_reads: continue
        rec = dict(beta_call=float(f[H["beta_call"]]), beta_read_mean=float(f[H["beta_read_mean"]]),
                   beta_read_median=float(f[H["beta_read_median"]]), n_reads=n,
                   n_calls=int(f[H["n_calls"]]), n_meth=int(f[H["n_meth"]]),
                   data_source=f[H["data_source"]])
        G4[(f[H["sample"]], f[H["slot"]], int(f[H["array_pos"]]))].append(rec)
    M = {}; src = collections.Counter(); nr = 0; ndup = 0
    for k, L in G4.items():
        if len(L) == 1:
            rec = L[0]
        else:
            ndup += 1
            if how == "last":
                rec = L[-1]
            elif how == "deep":
                rec = max(L, key=lambda x: x["n_reads"])
            else:
                nc = sum(x["n_calls"] for x in L); nm = sum(x["n_meth"] for x in L)
                w = float(sum(x["n_reads"] for x in L)) or 1.0
                rec = dict(beta_call=(nm / nc if nc else float("nan")),
                           beta_read_mean=sum(x["beta_read_mean"] * x["n_reads"] for x in L) / w,
                           beta_read_median=sum(x["beta_read_median"] * x["n_reads"] for x in L) / w,
                           n_reads=sum(x["n_reads"] for x in L), n_calls=nc, n_meth=nm,
                           data_source="+".join(sorted(x["data_source"] for x in L)))
        M[k] = rec
        src[rec["data_source"]] += 1; nr += 1
    if verbose:
        print("[v11_io 메틸] 합본 %s · 겹친 자리 %d → dedup=%s" % (os.path.basename(METHF), ndup, how))
        b = np.array([v["beta_call"] for v in M.values()])
        print("[v11_io 메틸] 자리 %d · 개체 %d · 하플로타입 %d · read>=%d" % (
            nr, len(set(k[0] for k in M)), len(set((k[0], k[1]) for k in M)), min_reads))
        print("   beta_call 평균 %.4f · 중앙 %.4f · <0.3 %.1f%% · >0.7 %.1f%%" % (
            b.mean(), np.median(b), 100 * (b < 0.3).mean(), 100 * (b > 0.7).mean()))
    return M
