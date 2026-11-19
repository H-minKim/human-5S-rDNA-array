# -*- coding: utf-8 -*-
# Arabidopsis 5S: the three contexts CG, CHG and CHH per read x unit, for the whole unit and for the coding part, the 5S gene, separately.
#   (The animal BAMs use a 5mCG model and carry no CHG or CHH calls; only Arabidopsis was basecalled with dorado in all contexts.)
#   The BAM reading and the unit and gene coordinates use the functions of x11_coding_readunit.py unchanged, loaded with exec, with only the log lines altered.
#   The context rule matches units_matrix_g8.py: for a forward C, the next base G gives CG, otherwise the second next G gives CHG, otherwise CHH;
#                                         for a reverse position, a reference G, the previous C gives CG, otherwise the second previous C gives CHG, otherwise CHH.
#   The MM tag is C+m. so an omitted C counts as an unmethylated call, as in x11; without that the result does not match the existing table.
#   Check: the whole-unit (n_calls, n_5mC) per context must equal the existing read_unit.tsv.
#   The call follows the x10 rule (a second peak at least 0.30 away, a trough no higher than half that peak, at least 2% of molecules on its side) plus the overdispersion ratio against the binomial null.
#   Condition: the number of calls is at least half the median for that array, context and interval, and at least 3.
import io, os, sys, glob, struct, collections
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'Environment variables')" % k)
    return v


OUT = _need("RDNA_PREV_OUT_DIR")
src = io.open(OUT + "/code/x11_coding_readunit.py", encoding="utf8").read()
L11 = 'LOG = io.open(OUT + ("/X11_CODING_PILOT.txt" if PILOT else "/X11_CODING_BUILD.txt"), "a", encoding="utf8", newline="\\n")'
assert L11 in src
G = {"__name__": "x11hdr"}
exec(src.split("STAT = collections.Counter()")[0].replace(L11, "LOG = io.StringIO()"), G)
src10 = io.open(OUT + "/code/x10_readunit_per_array.py", encoding="utf8").read()
L10 = 'LOG = io.open(OUT + "/X10_READUNIT_PERARRAY.txt", "w", encoding="utf8", newline="\\n")'
assert L10 in src10
G10 = {}
exec(src10.split("P(\"=\" * 104)")[0].replace(L10, "LOG = io.StringIO()"), G10)
call = G10["call"]
rcParams["font.family"] = "Malgun Gothic"
rcParams["axes.unicode_minus"] = False
rcParams["figure.dpi"] = 110
LOG = io.open(OUT + "/X13_ATH_CONTEXT.txt", "w", encoding="utf8", newline="\n")
rng = np.random.RandomState(13)


def P(s=""):
    try:
        print(s)
    except UnicodeEncodeError:
        print(s.encode("ascii", "replace").decode())
    LOG.write(s + "\n")
    LOG.flush()


B = G["B"]
SP = "arabidopsis_ColCC"
A = G["load_species"](SP)
NIBCOMP = G["NIBCOMP"]
STAT = collections.Counter()
acc = collections.defaultdict(lambda: np.zeros(4, dtype=np.int64))   # (rid, aid, pos, ctx) -> n, m, n_gene, m_gene
bam = glob.glob(B + "/04_bam/%s/*.5S.bam" % SP)[0]
bychrom = collections.defaultdict(list)
for aid, a in A.items():
    bychrom[a["chrom"]].append(aid)
CTX = np.array(["CG", "CHG", "CHH"])
for refs, r in G["bam_records"](bam):
    refid, pos0, lrn, mapq, _bin, ncig, flag, lseq = struct.unpack("<iiBBHHHi", r[:20])
    STAT["rec"] += 1
    if refid < 0 or flag & 0x904 or refs[refid] not in bychrom:
        continue
    chrom = refs[refid]
    rid = r[32:32 + lrn - 1].decode()
    p = 32 + lrn
    cig = np.frombuffer(r[p:p + 4 * ncig], dtype="<u4"); p += 4 * ncig
    sq = r[p:p + (lseq + 1) // 2]; p += (lseq + 1) // 2 + lseq
    T = G["parse_tags"](r, p)
    if "MM" not in T or "ML" not in T:
        STAT["no_mm"] += 1
        continue
    ops, lens = cig & 15, (cig >> 4).astype(np.int64)
    if (ops == 5).any():
        STAT["hardclip"] += 1
        continue
    q2r = np.full(lseq, -1, dtype=np.int64)
    q, rr = 0, pos0
    for o, n in zip(ops.tolist(), lens.tolist()):
        if o in (0, 7, 8):
            q2r[q:q + n] = np.arange(rr, rr + n); q += n; rr += n
        elif o in (1, 4):
            q += n
        elif o in (2, 3):
            rr += n
    b = np.frombuffer(sq, dtype=np.uint8)
    nib = np.empty(2 * len(b), dtype=np.uint8); nib[0::2] = b >> 4; nib[1::2] = b & 15
    codes = nib[:lseq]
    rev = bool(flag & 16)
    ocodes = NIBCOMP[codes][::-1] if rev else codes
    cpos = np.flatnonzero(ocodes == 2)
    ml = T["ML"]; mi = 0
    m_idx = m_ml = None
    for ent in T["MM"].rstrip(";").split(";"):
        if not ent:
            continue
        parts = ent.split(",")
        head = parts[0]
        mode = "?"
        if head[-1] in "?.":
            mode = head[-1]; head = head[:-1]
        base, strand_, mc = head[0], head[1], head[2:]
        d = np.array([int(x) for x in parts[1:]], dtype=np.int64)
        nc = len(mc) if not mc.isdigit() else 1
        vals = ml[mi:mi + nc * len(d)].reshape(-1, nc) if len(d) else np.zeros((0, nc), np.uint8)
        mi += nc * len(d)
        if base != "C" or strand_ != "+" or "m" not in mc:
            continue
        STAT["mode" + mode] += 1
        k = np.cumsum(d + 1) - 1 if len(d) else np.zeros(0, dtype=np.int64)
        if len(k) and k[-1] >= len(cpos):
            STAT["mm_overflow"] += 1
            break
        m_idx, m_ml = cpos[k], vals[:, mc.index("m")]
        if mode == ".":
            full = np.zeros(len(cpos), dtype=np.uint8)
            full[k] = m_ml
            m_idx, m_ml = cpos, full
    if m_idx is None:
        continue
    sidx = (lseq - 1 - m_idx) if rev else m_idx
    ref0 = q2r[sidx]
    ok = ref0 >= 0
    ref0, mlv = ref0[ok], m_ml[ok]
    for aid in bychrom[chrom]:
        a = A[aid]
        s = a["seq"]
        k0 = ref0 - (a["lo1"] - 1)
        inr = (k0 >= 2) & (k0 < len(s) - 2)
        if not inr.any():
            continue
        kk, vv, rp = k0[inr], mlv[inr], ref0[inr] + 1
        if rev:
            base_ok = s[kk] == 71
            ctx = np.where(s[kk - 1] == 67, 0, np.where(s[kk - 2] == 67, 1, 2))
        else:
            base_ok = s[kk] == 67
            ctx = np.where(s[kk + 1] == 71, 0, np.where(s[kk + 2] == 71, 1, 2))
        kk, vv, rp, ctx = kk[base_ok], vv[base_ok], rp[base_ok], ctx[base_ok]
        ui = np.searchsorted(a["us"], rp, side="right") - 1
        good = (ui >= 0) & (rp <= a["ue"][np.clip(ui, 0, None)])
        ui, vv, rp, ctx = ui[good], vv[good], rp[good], ctx[good]
        if len(ui) == 0:
            continue
        ingene = np.zeros(len(rp), dtype=bool)
        for u in np.unique(ui):
            for g0, g1 in a["genes"][u]:
                ingene |= (ui == u) & (rp >= g0) & (rp <= g1)
        callm = vv >= 128
        for u in np.unique(ui):
            if not a["genes"][u]:
                continue
            for c in range(3):
                w = (ui == u) & (ctx == c)
                if not w.any():
                    continue
                wg = w & ingene
                acc[(rid, aid, int(a["pos"][u]), CTX[c])] += np.array([w.sum(), (callm & w).sum(), wg.sum(), (callm & wg).sum()])

D = pd.DataFrame([(k[0], k[1], k[2], k[3], *v) for k, v in acc.items()],
                 columns=["read_id", "array_id", "array_pos", "context", "n_calls", "n_5mC", "n_calls_gene", "n_5mC_gene"])
D.to_csv(OUT + "/tables/ATH_READUNIT_CONTEXT.tsv.gz", sep="\t", index=False, compression="gzip")
P("Arabidopsis BAM records %d, stats %s, rows %d" % (STAT["rec"], dict(STAT), len(D)))

# ---- check against the existing read_unit.tsv, per context
old = pd.read_csv(glob.glob(B + "/03_reads/%s/*.read_unit.tsv*" % SP)[0], sep="\t")
for c in ["CG", "CHG", "CHH"]:
    o = old[old.context == c]
    n = D[D.context == c]
    M = o.merge(n, on=["read_id", "array_id", "array_pos"], how="inner", suffixes=("_old", ""))
    same = ((M.n_calls_old == M.n_calls) & (M.n_5mC_old == M.n_5mC)).mean()
    near = ((M.n_calls_old - M.n_calls).abs() <= 1).mean()
    P("  check %-3s: reference %d, new %d, shared %d, exact match %.4f, call-count difference <=1 %.4f, reference mean β %.3f / new %.3f"
      % (c, len(o), len(n), len(M), same, near, (M.n_5mC_old / M.n_calls_old).mean(), (M.n_5mC / M.n_calls).mean()))
P()

# ---- the call per array, context and interval
ST = pd.read_csv(OUT + "/tables/OTHERSP_ARRAY_STRUCT.tsv", sep="\t")
UL = dict(zip(ST.array_id, ST.unit_len))
NU = dict(zip(ST.array_id, ST.n_unit))
D["n_sp"] = D.n_calls - D.n_calls_gene
ROWS, PANEL = [], {}
for (aid, c), g in D.groupby(["array_id", "context"]):
    for part, ncol_, mcol_ in [("whole", "n_calls", "n_5mC"), ("coding", "n_calls_gene", "n_5mC_gene")]:
        gg = g[g[ncol_] > 0]
        med = gg[ncol_].median()
        gg = gg[gg[ncol_] >= max(3, 0.5 * med)]
        if len(gg) == 0:
            continue
        x = (gg[mcol_] / gg[ncol_]).values
        nn = gg[ncol_].values.astype(int)
        r, _, _ = call(x, float(np.median(nn)))
        key = gg.array_pos.values
        pmean = pd.Series(x).groupby(key).transform("mean").values
        sim = rng.binomial(nn, np.clip(pmean, 0, 1)) / nn
        od = np.var(x - pmean) / max(np.var(sim - pmean), 1e-12)
        r.update(array_id=aid, context=c, part=part, unit_len=UL.get(aid), n_unit=NU.get(aid), n=len(x),
                 calls=float(np.median(nn)), mean=x.mean(), sd=x.std(), off=(x < 0.35).mean() if c == "CG" else np.nan,
                 q10=np.quantile(x, 0.1), q90=np.quantile(x, 0.9), od=od)
        if len(x) < 100:
            r["call"] = "표본 부족 (" + r["call"] + ")"
        ROWS.append(r)
        PANEL[(aid, c, part)] = (x, r)
T = pd.DataFrame(ROWS)
T.to_csv(OUT + "/tables/ATH_CONTEXT_PERARRAY.tsv", sep="\t", index=False, float_format="%.5f")
P("  array                          unit_x_n  context  interval     read-units  call_median  mean   sd    10%%-90%%     overdispersion  peak1  peak2  trough_ratio  minor_side  verdict")
order = T.groupby("array_id").n.max().sort_values(ascending=False).index.tolist()
for aid in order:
    for r in T[T.array_id == aid].sort_values(["context", "part"], ascending=[True, False]).itertuples(index=False):
        P("  %-30s %3.0fx%-5.0f %-4s %-7s %8s  %5.0f   %.3f  %.3f  %.2f-%.2f   %5.1f  %.2f %s %s %.3f  %s"
          % (aid.split("__")[-1], r.unit_len, r.n_unit, r.context, r.part, format(r.n, ","), r.calls, getattr(r, "mean"),
             r.sd, r.q10, r.q90, r.od, r.peak1, ("%.2f" % r.peak2) if r.peak2 == r.peak2 else "  - ",
             ("%.2f" % r.trough_ratio) if r.trough_ratio == r.trough_ratio else "  - ", r.minor_mass, r.call))
P()
P("  verdict tally")
for (c, part), gg in T.groupby(["context", "part"]):
    P("   %-4s %-7s : %s" % (c, part, " · ".join("%s %d" % (k, v) for k, v in gg.call.value_counts().items())))

# ---- figure: one page per interval, whole unit and coding, with arrays as rows and contexts as columns
COLC = {"CG": "#7e57c2", "CHG": "#26a69a", "CHH": "#ef6c00"}
for part in ["whole", "coding"]:
    nr = len(order)
    fig, axs = plt.subplots(nr, 3, figsize=(13, 2.35 * nr + 0.9), squeeze=False)
    for i, aid in enumerate(order):
        for j, c in enumerate(["CG", "CHG", "CHH"]):
            ax = axs[i][j]
            if (aid, c, part) not in PANEL:
                ax.axis("off"); continue
            x, r = PANEL[(aid, c, part)]
            nb = int(round(r["calls"]))
            edges = np.clip((np.arange(nb + 2) - 0.5) / nb, 0, 1) if nb < 40 else np.linspace(0, 1, 41)
            ax.hist(x, bins=edges, weights=np.ones(len(x)) / len(x), color=COLC[c], alpha=0.85)
            ax.axvline(r["mean"], color="k", lw=0.8, ls="--")
            ax.set_xlim(0, 1)
            acc_ = aid.split("__")[-1]
            ax.set_title("%s:%s  %.0fbp x %.0f · %s\nread-유닛 %s · 호출 %.0f · 평균 %.2f · %s"
                         % (acc_.split("_")[0], acc_.split("_")[1], r["unit_len"], r["n_unit"], c, format(r["n"], ","),
                            r["calls"], r["mean"], r["call"]), fontsize=7.6)
            ax.tick_params(labelsize=6.5)
            if i == nr - 1:
                ax.set_xlabel("read 안 유닛 β (%s)" % c, fontsize=8)
    fig.suptitle("애기장대 5S — CG / CHG / CHH 각 read 의 각 유닛 메틸 분포 · %s · 행 = array (점선 = 평균)"
                 % ("유닛 전체" if part == "whole" else "coding(5S 유전자)만"), fontsize=12)
    fig.tight_layout(rect=[0, 0.01, 1, 0.975])
    fig.savefig(OUT + "/figures/Fig_ath_context_%s.png" % part)
    plt.close(fig)
    P("  saved: figures/Fig_ath_context_%s.png" % part)
LOG.close()
