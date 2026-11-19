# -*- coding: utf-8 -*-
# Extract the complete-unit data per array.
#   The filter and the call are the same as in the v13 pipeline (x9A), but the dump also carries the array_id, the unit index and the start coordinate,
#   and the per-position profile is accumulated per array as well, so that topics A, B, C and E can be redrawn on the complete-unit basis.
#   (Elsewhere topics A, B, C and E use the looser filter n_calls >= 0.5 times the array median.)
#   Outputs: tables/UNITBETA_PERARRAY_v13.tsv.gz and tables/PROFILE_BP_PERARRAY_v13.tsv.gz
#          (The PROFILE_*_v13X tables are a control that must match the x9A output; x9A's output is the canonical one.)
#   The v13 pipeline, as in x9A:
# Methylation profile within the unit, v13: the call is made by a test rather than a threshold.
#   Thresholds on Ashman's D, Sarle's bimodality coefficient and the minor-component weight, with a null gate, move with the number of units
#   (the overdispersed-null pass rate for the Sumatran orangutan is 3.3% on the full set against 13.7% on a 3,000-unit subsample).
#   v13 fits a one-component against a two-component beta-binomial mixture per group and requires all three of:
#     1) a significant parametric bootstrap likelihood ratio (p <= 0.01; the null LR distribution is simulated from the fitted one-component model)
#     2) split-half reproducibility (the component means differ by at most a quarter of their separation, and the low fraction by at most 0.05)
#     3) at least 2 peaks in the molecule-level density implied by the fit
#   If only 1) holds, the call is recorded as undecided and the split is made as for a unimodal group, at the nearest cut. The threshold call is kept in its own column.
# A unimodal group is split with nearest_cut rather than with a hard cut on the ratio (ratio_split).
#              With a hard cut, ties at the cut value all go to HIGH, so the achieved LOW fraction is always below the target of 0.2191 (median 0.175, largest shortfall 0.154).
#   The pass rate against an overdispersed (beta-binomial) null is also recorded: it is the stronger null, asking whether a single broad distribution,
#              rather than two states, would already satisfy the three conditions. It is reported in the table only and does not enter the call.
#   The threshold call and its null gate:
#   Measured: a group with only 7-8 CpG sites satisfies the three bimodality conditions
#   99-100% of the time even when it is genuinely unimodal (table BIMODAL_NULL_CHECK.tsv).
#   So for each group a unimodal (binomial) simulation with the same n_cpg and the same mean beta is run to measure the pass rate,
#   and at a pass rate of 5% or more the group is not accepted as bimodal but treated as unimodal with the fraction split f. The pass rate is reported.
#   (For the apes, zebra finch and Arabidopsis, which have many CpG sites, the pass rate is 0.0%.)
#   1) Bimodality is decided per group (array): Ashman's D > 2 (two-component Gaussian EM), Sarle's bimodality coefficient > 0.555, minor component >= 5%.
# Definition check: human 5S (P10) gives D 3.002, BC 0.8311, w_minor 0.306, matching the reference values 3.00 / 0.83 / 0.306
#   2) If bimodal, units are split by the likelihood ratio of that group's own two binomial components (EM); p_low, p_high and the equivalent cut are reported.
#   3) If unimodal, the human 5S LOW fraction f is applied to that group's ranking (a complementary split, with no gap left in the middle).
#      f comes from the human pool itself (25,644 complete units in P10): p_low 0.1580, p_high 0.7867, cut 0.4610, f = 0.2191.
# (matches the reference values to the decimal; table HUMAN_RATIO_F.tsv)
#   4) In zebra finch and zebrafish the contrast between a highly and a weakly methylated array is a difference between arrays, not bimodality of a pooled distribution, so the call is made per group.
#   A complete unit is defined once: a read-unit with n_cpg >= MINC that observes at least 90% of that unit's reference CpG sites.
#   - For each group, x10's own density/call functions are run on that group's complete-unit beta values to obtain the uni/bimodal call and the two peaks.
#   - If bimodal, the split is by likelihood ratio: d = n_meth*ln(a1/a0) + (n_cpg-n_meth)*ln((1-a1)/(1-a0)); d > 0 gives hi.
# - if unimodal, the top and bottom 25% by rank are used; the grouping column of the table records which was used
#   - Hard call at ML >= 128, a complete unit observes 90% of the reference CpG sites, positions in bp, depth-weighted.
#   A value threshold (q1c/q4c) fails in groups where more than 75% of units have beta = 1.00 (zebrafish 177 bp, zebra finch chr2/chr4):
#   the lower and upper cuts both become 1.0, every unit falls into q1 and q4 is empty.
#   So the bottom and top 25% are cut by sorted rank instead, with ties broken in a random order from seed 0.
# A "90% coverage of common positions" criterion holds only in a canonical frame where positions are shared between units, as in the human case.
#     In the other species the unit sequences differ slightly, so the same CpG lands on different bp positions and the set of common positions grows to 2-8 times the CpG count per unit, which no unit can cover to 90%.
#   The same intent (excluding partial reads) is therefore expressed against that unit's own reference CpG count:
#     a complete unit is a read-unit observing at least 90% of the CpG sites in that unit's reference sequence, with n_cpg >= MINC.
#   Note: in the human version the coordinate c is in the canonical frame, so positions are shared between units, whereas
#     here the unit lengths differ (gorilla 2,248 / 2,246 / 2,244 / 2,183 ...) and one CpG spreads over several bp positions,
#     so the set of common positions grows to 2-8 times the CpG count per unit and no unit covers 90% of it.
#   Units in frame are therefore selected first from the g8 map: unit_len within +-2 of that array's mode, cov_start <= 2 and cov_end >= unit_len - 2.
#   1) Hard call: each CpG is 0 or 1 at ML >= 128; probabilities are never averaged.
#   2) The unit of aggregation is the bp position itself; nothing is pre-binned.
#   3) beta is depth-weighted: per position, the sum of methylated calls over the sum of calls. Never a mean of per-unit means.
#   4) A minimum-depth gate per position (at plotting time, calls within a +-HALF bp window >= MINCALL), applied in the plotting script x61.
#   5) The running mean is applied only when plotting (a depth-weighted window measured in bp), in x61.
#   The complete-unit filter; without it a spurious hypomethylated peak appears:
#       (a) n_cpg >= MINC, and (b) only units covering at least 90% of the common positions (those observed at least 0.10 x the maximum per-position count).
#       Note: CpG sites per unit range from 7 to 190 across species, so MINC is set to half that group's median CpG count, capped at 20, and reported.
#   Groups: all, q1 (bottom 25%) and q4 (top 25%). These are rank cuts, not two clusters found in the data.
#   Coordinate frame: the g8 unit frame with the gene at the end of the unit, flipped to the array strand, and the gene interval is written in the same frame.
#   Outputs: tables/PROFILE_BP_v13X.tsv.gz, holding the calls, the methylated calls and beta per group, class and bp position,
#          and tables/PROFILE_META_v13X.tsv, holding the complete-unit count, the number discarded, the common positions, the cut, Ashman's D and the gene interval.
import io, os, sys, glob, re, struct, collections, math
from scipy import stats
import numpy as np
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'Environment variables')" % k)
    return v


NEW = _need("RDNA_OUT_DIR")
B = _need("RDNA_DATA_ROOT")
OLD = _need("RDNA_PREV_OUT_DIR")
sys.path.insert(0, NEW + "/code")
import arch  # noqa
import bbmix  # Weighted-likelihood fit of a one- and two-component beta-binomial (values identical to xa1; table BBMIX_CHECK.tsv).
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
LOG = io.open(NEW + "/logs/X9D_PERARRAY.txt", "w", encoding="utf8", newline="\n")
REASON = "완전유닛 array 단위 덤프(A·B·C·E 완전유닛판 자료)"
MLTHR, COVMIN = 128, 0.9
FRATIO = 0.2191
NSIM, NCAP = 100, 3000   # The null simulation of the threshold call: replicates per group and the cap on sampled units (seed 0).
# ---- The v13 rule: the thresholds (D, BC) are recorded only, and the call requires all three of the following.
BOOT = 99          # bootstrap replicates, so the floor on p is 1/(BOOT+1) = 0.01
LRMIN = 30.0       # An LR below this is non-significant without bootstrapping; the largest null LR measured was 15.5.
PALPHA = 0.01      # significance threshold (= the floor on p)
WTOL = 0.05        # tolerance on the difference in the low fraction between split halves
SEPDIV = 4.0       # tolerance on the difference in component means between split halves = (component separation) / SEPDIV
# ---- Further gates; the thresholds are chosen from the measured distributions (table V13_RULE_V2.tsv).
NULLMAX_OK = 100.0  # If the largest null LR exceeds this the fit has diverged and the test cannot be used; values up to 6e8 were observed.
SEPRES = 2.0        # Require the component separation to be at least SEPRES x (1 / median CpG count), which excludes components created by quantisation.
MINW = 0.05 # lower bound on the minor-component weight, so that a 1.7% tail is not counted as a peak
VALLEY = 0.6        # An upper bound on the trough-to-peak ratio, so that a shoulder above 0.9 is not counted as two peaks.
NSPLIT, STARTS = 5, 5   # Five split-half replicates: with one, a local optimum in the fit can flip the call.


def em_two_binomial(kn, iters=200): # two-component binomial EM
    p0, p1, w = 0.20, 0.80, 0.40
    for _ in range(iters):
        a0 = b0 = a1 = b1 = sw = 0.0
        for k, n in kn:
            l0 = k * math.log(p0) + (n - k) * math.log(1 - p0) + math.log(max(w, 1e-12))
            l1 = k * math.log(p1) + (n - k) * math.log(1 - p1) + math.log(max(1 - w, 1e-12))
            m = max(l0, l1)
            e0, e1 = math.exp(l0 - m), math.exp(l1 - m)
            r0 = e0 / (e0 + e1)
            sw += r0; a0 += r0 * k; b0 += r0 * n; a1 += (1 - r0) * k; b1 += (1 - r0) * n
        w = sw / len(kn)
        p0 = a0 / b0 if b0 else p0
        p1 = a1 / b1 if b1 else p1
        if p0 > p1:
            p0, p1, w = p1, p0, 1 - w
    return p0, p1, w


def lr_split(kn, p0, p1):
    cut = (math.log((1 - p0) / (1 - p1)) / (math.log(p1 / p0) + math.log((1 - p0) / (1 - p1))))
    low = {i for i, (k, n) in enumerate(kn)
           if k * math.log(p1 / p0) + (n - k) * math.log((1 - p1) / (1 - p0)) <= 0}
    return low, cut


def ratio_split(betas, f):
    """Hard-cut split on the ratio. Kept for comparison; nearest_cut below is used instead."""
    b = sorted(betas)
    i = max(1, min(len(b) - 1, int(round(len(b) * f))))
    return b[i]


def nearest_cut(betas, f):
    """The cut whose attained fraction is closest to the target f.
       A hard cut on the ratio sends every tie at the cut value to HIGH,
       so the attained LOW fraction stays below the target.
       Every distinct β value is a candidate cut, and the c whose (β < c) fraction is closest to f is chosen.
       Ties are not split (units with the same value always fall in the same group). Between equal candidates the smaller cut is chosen, for reproducibility."""
    b = np.asarray(betas, dtype=float)
    uq = np.unique(b)
    fr = np.array([float((b < c).mean()) for c in uq])
    j = int(np.argmin(np.abs(fr - f)))
    return float(uq[j]), float(fr[j])


def betabinom_null_rate(kk, nn, nsim, rng, iters=60):
    """Overdispersed null: does a single broad (continuous) distribution, not 'two states', also pass the 3 conditions?
       A unimodal binomial null holds sampling error only, so its variance is too small (no biological variation). A beta-binomial is therefore
       fitted to the same (k,n) by the method of moments and simulated. Returns : (pass rate, α, β, usable)
       With α<1 or β<1 the beta itself is U-shaped (bimodal) and is not a 'unimodal null', so it is flagged as not usable."""
    p = kk / nn
    m = float(kk.sum() / nn.sum())
    s2 = float(np.var(p, ddof=1)) if len(p) > 1 else 0.0
    inv = float(np.mean(1.0 / nn))
    denom = m * (1 - m) * (1 - inv)
    rho = (s2 - m * (1 - m) * inv) / denom if denom > 1e-12 else 0.0
    rho = float(min(max(rho, 1e-6), 0.999))
    t = 1.0 / rho - 1.0
    a, bta = m * t, (1 - m) * t
    usable = bool((a >= 1.0) and (bta >= 1.0))
    hits = 0
    for _ in range(nsim):
        pp = rng.beta(a, bta, size=len(nn))
        sim = rng.binomial(nn.astype(int), pp) / nn
        if len(sim) > 4 and sim.std() > 1e-9:
            g1s = float(stats.skew(sim, bias=False)); kus = float(stats.kurtosis(sim, bias=False))
            ns = len(sim)
            BCs = (g1s ** 2 + 1) / (kus + 3 * ((ns - 1) ** 2) / ((ns - 2) * (ns - 3)))
            ms = np.array([np.quantile(sim, .2), np.quantile(sim, .8)]); sss = np.array([sim.std(), sim.std()]) + 1e-6
            wws = np.array([.5, .5])
            for _i in range(iters):
                pdfs = wws / (sss * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((sim[:, None] - ms) / sss) ** 2)
                rrs = pdfs / np.clip(pdfs.sum(1, keepdims=True), 1e-300, None)
                nks = rrs.sum(0) + 1e-12
                ms = (rrs * sim[:, None]).sum(0) / nks
                sss = np.sqrt((rrs * (sim[:, None] - ms) ** 2).sum(0) / nks) + 1e-6
                wws = nks / ns
            ADs = float(np.sqrt(2) * abs(ms[0] - ms[1]) / np.sqrt(sss[0] ** 2 + sss[1] ** 2))
            hits += (ADs > 2 and BCs > 0.555 and min(wws) >= 0.05)
    return hits / float(nsim), float(a), float(bta), usable


def P(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


src = io.open(OLD + "/code/x11_coding_readunit.py", encoding="utf8").read()
L11 = 'LOG = io.open(OUT + ("/X11_CODING_PILOT.txt" if PILOT else "/X11_CODING_BUILD.txt"), "a", encoding="utf8", newline="\\n")'
head = src.split("STAT = collections.Counter()")[0].replace(L11, "LOG = io.StringIO()")
G = {"__name__": "x11hdr"}
exec(compile(head, OLD + "/code/x11_coding_readunit.py", "exec"), G)
NIBCOMP, load_species, bam_records, parse_tags, read_fa = (
    G["NIBCOMP"], G["load_species"], G["bam_records"], G["parse_tags"], G["read_fa"])
# ---- x10's density and call functions are imported verbatim rather than re-implemented here.
_x10 = io.open(NEW + "/code/x10_readunit_per_array.py", encoding="utf8").read()
_hdr = _x10.split("def density(")[0]
_body = "def density(" + _x10.split("def density(")[1].split('P("=" * 104)')[0]
_ns = {"__name__": "x10fn"}
_pre = "import numpy as np\n" + "\n".join([l for l in _hdr.split("\n") if re.match(r"^SEP, TROUGH, MINOR *=", l)])
exec(compile(_pre + "\n" + _body, "x10_fn", "exec"), _ns)
x10_call = _ns["call"]

BC = pd.read_csv(NEW + "/tables/BONOBO_READ_BASECALLER.tsv", sep="\t")
B637 = set(BC.loc[BC.basecaller.astype(str).str.contains("6.3.7"), "read_id"])
MT3 = pd.read_csv(NEW + "/tables/PERARRAY_META_v3.tsv", sep="\t")      # The group definition is taken unchanged from x52.
GRP_OF = {}
for t in MT3.itertuples(index=False):
    for a in str(t.array목록).split("|"):
        GRP_OF[a] = t.묶음
GSPEC = {t.묶음: dict(종=t.종, 유닛길이=int(t.유닛길이)) for t in MT3.itertuples(index=False)}
SPK = [("고릴라", "gorilla_mGorGor1"), ("보노보", "bonobo_mPanPan1"), ("침팬지", "chimp_mPanTro3"),
       ("보르네오오랑우탄", "orang_bornean_mPonPyg2"), ("수마트라오랑우탄", "orang_sumatran_mPonAbe1"),
       ("시아망", "siamang_mSymSyn1"), ("애기장대", "arabidopsis_ColCC"),
       ("금화조", "zfinch_v03_dip"), ("제브라피시", "zebrafish_GRCz12tu")]


def ashman(x):
    """Ashman D = sqrt(2)|m1-m2| / sqrt(s1^2+s2^2) from a two-component normal mixture fitted by EM"""
    x = np.asarray(x, float)
    if len(x) < 50 or x.std() < 1e-6:
        return np.nan
    m = np.array([np.quantile(x, 0.2), np.quantile(x, 0.8)])
    s = np.array([x.std(), x.std()]) + 1e-6
    w = np.array([0.5, 0.5])
    for _ in range(200):
        pdf = w[None, :] / (s[None, :] * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x[:, None] - m[None, :]) / s[None, :]) ** 2)
        tot = pdf.sum(axis=1, keepdims=True)
        tot[tot == 0] = 1e-300
        r = pdf / tot
        nk = r.sum(axis=0) + 1e-12
        m = (r * x[:, None]).sum(axis=0) / nk
        s = np.sqrt((r * (x[:, None] - m[None, :]) ** 2).sum(axis=0) / nk) + 1e-6
        w = nk / len(x)
    return float(np.sqrt(2) * abs(m[0] - m[1]) / np.sqrt(s[0] ** 2 + s[1] ** 2))


ROWS, META, UROWS, AROWS = [], [], [], [] # UROWS = per unit (CpG count, methylated count, beta, group), for the "cut line on the beta distribution" panel
for ko, sp in SPK:
    Asp = load_species(sp)
    STR, FRAMEOK, FRAMEN = {}, {}, {}
    for f_ in sorted(glob.glob(B + "/02_consensus_map/%s/*.g8map.tsv" % sp)):
        if not re.search(r"_[0-9]+_[0-9]+[.]g8map[.]tsv$", f_.replace("\\", "/")):
            continue
        g_ = pd.read_csv(f_, sep="\t")
        aid_ = g_.array_id.iloc[0]
        STR[aid_] = g_.strand.iloc[0]
        # Units in frame only: length within +-2 of that array's mode, and map coverage spanning the whole unit.
        modal = int(g_.unit_len.mode().iloc[0])
        ok_ = g_[(g_.unit_len.sub(modal).abs() <= 2) & (g_.cov_start <= 2) & (g_.cov_end >= g_.unit_len - 2)]
        FRAMEOK[aid_] = set(int(x) for x in ok_.array_pos)
        FRAMEN[aid_] = (len(ok_), len(g_), modal)
    REFCPG = {}                      # (array, unit index) -> the CpG count of that unit's reference sequence
    for a in [x for x in Asp if x in GRP_OF]:
        v = Asp[a]; sq = v["seq"]
        cg = np.flatnonzero((sq[:-1] == 67) & (sq[1:] == 71)) + 1 + (v["lo1"] - 1)   # 1-based genome coordinate
        ui_ = np.searchsorted(v["us"], cg, side="right") - 1
        ok2 = (ui_ >= 0) & (cg <= v["ue"][np.clip(ui_, 0, None)])
        cnt_ = collections.Counter(ui_[ok2].tolist())
        REFCPG[a] = cnt_
    keep_aids = [a for a in Asp if a in GRP_OF]
    bychrom = collections.defaultdict(list)
    for a in keep_aids:
        bychrom[Asp[a]["chrom"]].append(a)
    # Collect read-units per group: (n_cpg, n_meth, position bits, methylation bits).
    UNITS = collections.defaultdict(list)
    for bam in sorted(glob.glob(B + "/04_bam/%s/*.5S.bam" % sp)):
        for refs, r in bam_records(bam):
            refid, pos0, lrn, mapq, _bin, ncig, flag, lseq = struct.unpack("<iiBBHHHi", r[:20])
            if refid < 0 or (flag & 0x904) or refs[refid] not in bychrom:
                continue
            rid = r[32:32 + lrn - 1].decode()
            if sp == "bonobo_mPanPan1" and rid in B637:
                continue
            p = 32 + lrn
            cig = np.frombuffer(r[p:p + 4 * ncig], dtype="<u4"); p += 4 * ncig
            sq = r[p:p + (lseq + 1) // 2]; p += (lseq + 1) // 2 + lseq
            T = parse_tags(r, p)
            if "MM" not in T or "ML" not in T:
                continue
            ops, lens = cig & 15, (cig >> 4).astype(np.int64)
            if (ops == 5).any():
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
            bb = np.frombuffer(sq, dtype=np.uint8)
            nib = np.empty(2 * len(bb), dtype=np.uint8); nib[0::2] = bb >> 4; nib[1::2] = bb & 15
            codes = nib[:lseq]
            rev = bool(flag & 16)
            ocodes = NIBCOMP[codes][::-1] if rev else codes
            cpos = np.flatnonzero(ocodes == 2)
            ml = T["ML"]; mi = 0
            m_idx = m_ml = None
            for ent in T["MM"].rstrip(";").split(";"):
                if not ent:
                    continue
                parts = ent.split(","); hd2 = parts[0]; mode = "?"
                if hd2[-1] in "?.":
                    mode = hd2[-1]; hd2 = hd2[:-1]
                base, st_, mc = hd2[0], hd2[1], hd2[2:]
                d = np.array([int(x) for x in parts[1:]], dtype=np.int64)
                nc = len(mc) if not mc.isdigit() else 1
                vals = ml[mi:mi + nc * len(d)].reshape(-1, nc) if len(d) else np.zeros((0, nc), np.uint8)
                mi += nc * len(d)
                if base != "C" or st_ != "+" or "m" not in mc or len(d) == 0:
                    continue
                k = np.cumsum(d + 1) - 1
                if k[-1] >= len(cpos):
                    break
                m_idx, m_ml = cpos[k], vals[:, mc.index("m")]
                if mode == ".":
                    full = np.zeros(len(cpos), dtype=np.uint8); full[k] = m_ml
                    m_idx, m_ml = cpos, full
            if m_idx is None:
                continue
            sidx = (lseq - 1 - m_idx) if rev else m_idx
            ref0 = q2r[sidx]
            ok = ref0 >= 0
            ref0, mlv = ref0[ok], m_ml[ok]
            for aid in bychrom[refs[refid]]:
                v, s_ = Asp[aid], Asp[aid]["seq"]
                key = GRP_OF[aid]
                L = GSPEC[key]["유닛길이"]
                k0 = ref0 - (v["lo1"] - 1)
                inr = (k0 >= 1) & (k0 < len(s_) - 1)
                if not inr.any():
                    continue
                kk, vv, rp = k0[inr], mlv[inr], ref0[inr] + 1
                cg = (s_[kk] == 71) & (s_[kk - 1] == 67) if rev else (s_[kk] == 67) & (s_[kk + 1] == 71)
                vv, rp = vv[cg], rp[cg]
                if not len(rp):
                    continue
                call = (vv >= MLTHR).astype(np.uint8)          # 1) hard call
                ui = np.searchsorted(v["us"], rp, side="right") - 1
                good = (ui >= 0) & (rp <= v["ue"][np.clip(ui, 0, None)])
                if not good.any():
                    continue
                ui2, call2, rp2 = ui[good], call[good], rp[good]
                inu = (rp2 - v["us"][ui2] + 1) if STR[aid] == "+" else (v["ue"][ui2] - rp2 + 1)
                for u in np.unique(ui2):
                    w = ui2 == u
                    pos = np.clip(inu[w].astype(np.int64), 1, L)
                    cl = call2[w]
                    bm_c = np.zeros(L + 2, dtype=bool); bm_m = np.zeros(L + 2, dtype=bool)
                    bm_c[pos] = True
                    bm_m[pos[cl.astype(bool)]] = True
                    UNITS[key].append((int(w.sum()), int(cl.sum()), np.packbits(bm_c), np.packbits(bm_m), int(REFCPG[aid].get(int(u), 0)), aid, int(u), int(v["us"][u])))
    # ---------------- per group: filter to complete units, then accumulate per position
    for key, lst in UNITS.items():
        L = GSPEC[key]["유닛길이"]
        nbits = L + 2
        TOT = np.zeros(nbits, dtype=np.int64)
        for nc, nm, bc, bmm, nref, _a, _u, _st in lst:
            TOT += np.unpackbits(bc)[:nbits]
        mx = TOT.max() if TOT.size else 0
        CPGB = np.flatnonzero(TOT >= 0.10 * mx)                 # a common position is one observed at least 0.10 x the maximum per-position count
        med_cpg = float(np.median([x[0] for x in lst])) if lst else 0.0
        med_ref = float(np.median([x[4] for x in lst])) if lst else 0.0
        MINC = int(min(20, max(3, round(0.5 * med_cpg))))       # adapted per species because the CpG count differs; reported in the table
        kept, keptinfo = [], []      # keptinfo = (array_id, unit index, start coordinate), for the complete-unit versions of A, B, C and E
        drop_nc = drop_cov = 0
        for nc, nm, bc, bmm, nref, aid_, u_, st_ in lst:
            if nc < MINC:
                drop_nc += 1
                continue
            if nref > 0 and nc / float(nref) < COVMIN:        # a complete unit observes at least 90% of the reference CpG sites
                drop_cov += 1
                continue
            bits = np.unpackbits(bc)[:nbits].astype(bool)
            kept.append((nc, nm, bits, np.unpackbits(bmm)[:nbits].astype(bool)))
            keptinfo.append((aid_, u_, st_))
        if not kept:
            P("  %s: complete units 0 (discarded nc %d, cov %d): skipped" % (key, drop_nc, drop_cov))
            continue
        bvals = np.array([nm / float(nc) for nc, nm, _b, _m in kept])
        bb = np.sort(bvals)
        q1c, q4c = float(bb[len(bb) // 4]), float(bb[(3 * len(bb)) // 4])
        # ---- Bimodality rule, three conditions (Gaussian EM Ashman D, Sarle BC, minor component)
        nunit = len(bvals)
        g1 = float(stats.skew(bvals, bias=False)); ku = float(stats.kurtosis(bvals, bias=False))
        BC = (g1 ** 2 + 1) / (ku + 3 * ((nunit - 1) ** 2) / ((nunit - 2) * (nunit - 3))) if nunit > 4 else np.nan
        mm = np.array([np.quantile(bvals, .2), np.quantile(bvals, .8)])
        ss = np.array([bvals.std(), bvals.std()]) + 1e-6
        ww = np.array([.5, .5])
        for _ in range(300):
            pdf = ww / (ss * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((bvals[:, None] - mm) / ss) ** 2)
            rr = pdf / np.clip(pdf.sum(1, keepdims=True), 1e-300, None)
            nk = rr.sum(0) + 1e-12
            mm = (rr * bvals[:, None]).sum(0) / nk
            ss = np.sqrt((rr * (bvals[:, None] - mm) ** 2).sum(0) / nk) + 1e-6
            ww = nk / nunit
        AD = float(np.sqrt(2) * abs(mm[0] - mm[1]) / np.sqrt(ss[0] ** 2 + ss[1] ** 2))
        wmin = float(min(ww))
        cond3 = (AD > 2) and (BC > 0.555) and (wmin >= 0.05)
        # ---- Null simulation: build a unimodal (binomial) set with the same n_cpg and mean, and measure the pass rate of the three conditions.
        kk_ = np.array([nm for _nc, nm, _b, _m in kept], dtype=float)
        nn_ = np.array([nc for nc, _nm, _b, _m in kept], dtype=float)
        rngs = np.random.default_rng(0)
        if len(nn_) > NCAP:
            sel_ = rngs.choice(len(nn_), NCAP, replace=False)
            kk_, nn_ = kk_[sel_], nn_[sel_]
        pbar = float(kk_.sum() / nn_.sum())
        hits = 0
        for _ in range(NSIM):
            sim = rngs.binomial(nn_.astype(int), pbar) / nn_
            if len(sim) > 4 and sim.std() > 1e-9:
                g1s = float(stats.skew(sim, bias=False)); kus = float(stats.kurtosis(sim, bias=False))
                ns = len(sim)
                BCs = (g1s ** 2 + 1) / (kus + 3 * ((ns - 1) ** 2) / ((ns - 2) * (ns - 3)))
                ms = np.array([np.quantile(sim, .2), np.quantile(sim, .8)]); sss = np.array([sim.std(), sim.std()]) + 1e-6
                wws = np.array([.5, .5])
                for _i in range(60):
                    pdfs = wws / (sss * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((sim[:, None] - ms) / sss) ** 2)
                    rrs = pdfs / np.clip(pdfs.sum(1, keepdims=True), 1e-300, None)
                    nks = rrs.sum(0) + 1e-12
                    ms = (rrs * sim[:, None]).sum(0) / nks
                    sss = np.sqrt((rrs * (sim[:, None] - ms) ** 2).sum(0) / nks) + 1e-6
                    wws = nks / ns
                ADs = float(np.sqrt(2) * abs(ms[0] - ms[1]) / np.sqrt(sss[0] ** 2 + sss[1] ** 2))
                hits += (ADs > 2 and BCs > 0.555 and min(wws) >= 0.05)
        nullrate = hits / float(NSIM)
        gate11 = cond3 and (nullrate < 0.05)       # The threshold calls, recorded for comparison only.
        # ---- The overdispersed null, reported only: would a single broad distribution also satisfy the three conditions?
        bbrate, bba, bbb, bbok = betabinom_null_rate(kk_, nn_, NSIM, np.random.default_rng(1))
        # ================= the v13 call: a test, not a threshold =================
        #   1) The parametric bootstrap likelihood ratio of a one- against a two-component beta-binomial mixture must be significant (p <= PALPHA).
        #   2) Split-half reproducibility: on a random half-split the two fitted components must agree to within a quarter of their separation, and the low fractions to within WTOL.
        #   3) The molecule-level density implied by the fit must have at least 2 peaks; two components is not the same as two peaks.
        kk_all = np.array([nm for _nc, nm, _b, _m in kept], dtype=float)
        nn_all = np.array([nc for nc, _nm, _b, _m in kept], dtype=float)
        ku, nu, cu = bbmix.pack(kk_all, nn_all)
        LRobs, _lg1, _lg2, ab1, prm = bbmix.lr_stat(ku, nu, cu)
        if LRobs >= LRMIN:                          # A small LR is non-significant without bootstrapping; the largest null LR measured was 15.5.
            rngb = np.random.default_rng(12345)
            nl = []
            for _ in range(BOOT):
                _ks, _ns, _cs = bbmix.boot_null(nn_all, ab1[0], ab1[1], rngb)
                nl.append(bbmix.lr_stat(_ks, _ns, _cs, starts=1)[0])
            nl = np.array(nl)
            pval = (np.sum(nl >= LRobs) + 1.0) / (BOOT + 1.0)
            nullmax = float(nl.max())
        else:
            pval, nullmax = 1.0, np.nan
        sep = float(prm["mu_high"] - prm["mu_low"])
        minw = float(min(prm["w_low"], 1 - prm["w_low"]))
        npk, vratio, pk0, pk1 = bbmix.mixture_peaks(prm)
        sepres = sep * med_cpg                       # how many times the resolution (1 / CpG count) the component separation is
        # ---- The gates are checked in order and each failure gets its own name, rather than one undifferentiated "undecided".
        sig = bool(pval <= PALPHA)
        dlo = dhi = dw = np.nan
        npass = -1
        if not (np.isfinite(nullmax) and nullmax < NULLMAX_OK) and LRobs >= LRMIN:
            verdict = "검정불가(귀무 적합 폭주)"
        elif not sig:
            verdict = "단봉(2성분 유의하지 않음)"
        elif sepres < SEPRES:
            verdict = "판정불가(분해능 부족 · 간격이 1/CpG 의 %.2f 배)" % sepres
        elif npk < 2 or (np.isfinite(vratio) and vratio > VALLEY):
            verdict = "단봉(봉우리 하나%s)" % ("" if npk < 2 else " · 어깨 골깊이비 %.2f" % vratio)
        elif minw < MINW:
            verdict = "단봉+소수꼬리(소수 성분 %.3f)" % minw
        else:
            dlos, dhis, dws, npass = [], [], [], 0
            for _s in range(NSPLIT):               # five replicates for reproducibility
                rh = np.random.default_rng(100 + _s)
                ix_ = rh.permutation(len(nn_all))
                hs = []
                for h_ in (ix_[:len(ix_) // 2], ix_[len(ix_) // 2:]):
                    _kh, _nh, _ch = bbmix.pack(kk_all[h_], nn_all[h_])
                    hs.append(bbmix.fit2(_kh, _nh, _ch, starts=STARTS)[1])
                a_ = abs(hs[0]["mu_low"] - hs[1]["mu_low"]); b_ = abs(hs[0]["mu_high"] - hs[1]["mu_high"])
                c_ = abs(hs[0]["w_low"] - hs[1]["w_low"])
                dlos.append(a_); dhis.append(b_); dws.append(c_)
                npass += int((a_ <= sep / SEPDIV) and (b_ <= sep / SEPDIV) and (c_ <= WTOL))
            dlo, dhi, dw = max(dlos), max(dhis), max(dws)
            verdict = "이봉" if npass == NSPLIT else "보류(재현 %d/%d)" % (npass, NSPLIT)
        is_bi = (verdict == "이봉")
        repro = (npass == NSPLIT)
        lab = np.full(nunit, "MID", dtype=object)
        if is_bi:
            kn = [(nm, nc) for nc, nm, _b, _m in kept]
            p0, p1, wlow = em_two_binomial(kn)
            lowset, cutv = lr_split(kn, p0, p1)
            for i in range(nunit):
                lab[i] = "q1" if i in lowset else "q4"
            method = ("우도비(이봉 · p %.4f/%.4f · 컷 %.4f · LRT p<=%.3f · LR %.0f)"
                      % (p0, p1, cutv, pval, LRobs))
        else:
            p0 = p1 = np.nan
            cutv, achv = nearest_cut(bvals, FRATIO)      # the cut whose achieved fraction comes closest to the target
            lab = np.where(bvals < cutv, "q1", "q4").astype(object)
            method = ("사람 5S 비율 목표 f=%.4f · 가장 가까운 컷 %.4f (달성 %.4f · %s)"
                      % (FRATIO, cutv, achv, verdict))
        acc = {g: np.zeros((2, nbits)) for g in ("all", "q1", "q4")}
        accA = collections.defaultdict(lambda: np.zeros((2, nbits)))   # per-array accumulation by position, for the complete-unit versions of B, C and E
        cnt = collections.Counter()
        cntA = collections.Counter()
        for (nc, nm, bits, mb), lb, (aid_, u_, st_) in zip(kept, lab, keptinfo):
            UROWS.append((GSPEC[key]["종"], key, aid_, int(u_), int(st_), int(nc), int(nm),
                          round(nm / float(nc), 6), lb))
            gs = ["all"] + ([lb] if lb in ("q1", "q4") else [])
            for g in gs:
                acc[g][0] += bits
                acc[g][1] += mb
                cnt[g] += 1
                accA[(aid_, g)][0] += bits
                accA[(aid_, g)][1] += mb
                cntA[(aid_, g)] += 1
        for g in ("all", "q1", "q4"):
            n, m = acc[g]
            for p_ in np.flatnonzero(n > 0):
                ROWS.append((GSPEC[key]["종"], key, g, L, int(p_), int(n[p_]), int(m[p_]), round(m[p_] / n[p_], 6)))
        for (aid_, g), arr in accA.items():                 # the per-position profile, per array
            n, m = arr
            for p_ in np.flatnonzero(n > 0):
                AROWS.append((GSPEC[key]["종"], key, aid_, g, L, int(p_), int(n[p_]), int(m[p_]),
                              round(m[p_] / n[p_], 6), int(cntA[(aid_, g)])))
        # the gene interval, in the same frame
        aid0 = [a for a in GRP_OF if GRP_OF[a] == key][0]
        v0 = Asp[aid0] if aid0 in Asp else None
        gw = []
        if v0 is not None:
            for i in range(len(v0["us"])):
                for (x, y) in v0["genes"][i]:
                    gw.append((x - v0["us"][i] + 1, y - v0["us"][i] + 1) if STR[aid0] == "+"
                              else (v0["ue"][i] - y + 1, v0["ue"][i] - x + 1))
        fr = [FRAMEN[a] for a in GRP_OF if GRP_OF[a] == key and a in FRAMEN]
        META.append(dict(종=GSPEC[key]["종"], 묶음=key, 유닛길이=L, 원유닛=len(lst), 완전유닛=len(kept),
                         틀맞춤유닛=sum(x[0] for x in fr), 지도유닛=sum(x[1] for x in fr),
                         버림_CpG적음=drop_nc, 버림_덮임부족=drop_cov, MINC=MINC, 흔한자리=int(len(CPGB)),
                         최대자리관측=int(mx), CpG중앙=med_cpg, 참조CpG중앙=med_ref, q1컷=round(q1c, 4), q4컷=round(q4c, 4), 무리방법=method,
                         판정=verdict.split("(")[0].strip(), 판정상세=verdict,   # call = bimodal / unimodal / unimodal with a minor tail / undecided / not callable / not testable
                         LRT_LR=round(float(LRobs), 2), LRT_귀무최대=(round(nullmax, 2) if np.isfinite(nullmax) else np.nan),
                         LRT_p=round(float(pval), 4), LRT_B=BOOT,
                         성분_저=round(prm["mu_low"], 4), 성분_고=round(prm["mu_high"], 4), 성분_저몫=round(prm["w_low"], 4),
                         소수성분몫=round(minw, 4), 간격배_분해능=round(float(sepres), 3),
                         반반_저차최대=(round(float(dlo), 4) if np.isfinite(dlo) else np.nan),
                         반반_고차최대=(round(float(dhi), 4) if np.isfinite(dhi) else np.nan),
                         반반_몫차최대=(round(float(dw), 4) if np.isfinite(dw) else np.nan),
                         반반허용=round(float(sep / SEPDIV), 4),
                         재현벌=("%d/%d" % (npass, NSPLIT) if npass >= 0 else "—"),
                         적합봉수=int(npk), 골깊이비=(round(float(vratio), 4) if np.isfinite(vratio) else np.nan),
                         한성분_알파=round(float(ab1[0]), 3), 한성분_베타=round(float(ab1[1]), 3),
                         옛판정_v12=("이봉" if gate11 else "단봉"), three조건=("통과" if cond3 else "탈락"), 귀무통과율=round(nullrate, 3),
                         과분산귀무통과율=round(bbrate, 3), 과분산귀무_알파=round(bba, 3), 과분산귀무_베타=round(bbb, 3),
                         과분산귀무_쓸수있나=("예" if bbok else "아니오(α·β<1 이라 귀무가 U자)"),
                         분할방식=("우도비" if is_bi else "가장 가까운 컷"),
                         달성LOW비율=round(float((lab == "q1").mean()), 4),
                         AshmanD=round(AD, 3), BC=round(BC, 4), 소수성분=round(wmin, 4),
                         컷=round(float(cutv), 4), p_low=(round(p0, 4) if np.isfinite(p0) else np.nan), p_high=(round(p1, 4) if np.isfinite(p1) else np.nan),
                         
                         n_all=cnt["all"], n_q1=cnt["q1"], n_q4=cnt["q4"],
                         유전자시작=(float(np.median([a for a, _b in gw])) if gw else np.nan),
                         유전자끝=(float(np.median([b for _a, b in gw])) if gw else np.nan)))
    P("%-10s groups %2d, complete units %s (discarded %s)"
      % (ko, sum(1 for m in META if m["종"] == ko),
         format(sum(m["완전유닛"] for m in META if m["종"] == ko), ","),
         format(sum(m["버림_CpG적음"] + m["버림_덮임부족"] for m in META if m["종"] == ko), ",")))

R = pd.DataFrame(ROWS, columns=["종", "묶음", "무리", "유닛길이", "위치", "호출수", "메틸수", "beta"])
M = pd.DataFrame(META)
TMP = os.environ.get("TEMP", ".") + "/x60"
os.makedirs(TMP, exist_ok=True)
R.to_csv(TMP + "/PROFILE_BP_v13X.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
M.to_csv(TMP + "/PROFILE_META_v13X.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
U = pd.DataFrame(UROWS, columns=["종", "묶음", "array_id", "유닛번호", "시작좌표", "n_cpg", "n_meth", "beta", "무리"])
AP = pd.DataFrame(AROWS, columns=["종", "묶음", "array_id", "무리", "유닛길이", "위치", "호출수", "메틸수", "beta", "유닛수"])
AP.to_csv(TMP + "/PROFILE_BP_PERARRAY_v13.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/PROFILE_BP_PERARRAY_v13.tsv.gz", "tables/PROFILE_BP_PERARRAY_v13.tsv.gz", REASON)
U.to_csv(TMP + "/UNITBETA_PERARRAY_v13.tsv.gz", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/PROFILE_BP_v13X.tsv.gz", "tables/PROFILE_BP_v13X.tsv.gz", REASON)
arch.put(TMP + "/PROFILE_META_v13X.tsv", "tables/PROFILE_META_v13X.tsv", REASON)
arch.put(TMP + "/UNITBETA_PERARRAY_v13.tsv.gz", "tables/UNITBETA_PERARRAY_v13.tsv.gz", REASON)
P("")
P("groups %d, position rows %s" % (len(M), format(len(R), ",")))
P("")
P("apes: gene region vs spacer (complete units, depth-weighted)")
if not len(M):
    P("  no complete unit at all: review the filter"); LOG.close(); sys.exit(0)
for t in M[M.종.isin(["고릴라", "보노보", "침팬지", "보르네오오랑우탄", "수마트라오랑우탄", "시아망"])].itertuples(index=False):
    d = R[R.묶음 == t.묶음]
    for g in ("all", "q1", "q4"):
        dd = d[d.무리 == g]
        ing = dd[(dd.위치 >= t.유전자시작) & (dd.위치 <= t.유전자끝)]
        out = dd[dd.위치 < t.유전자시작]
        gb = ing.메틸수.sum() / max(ing.호출수.sum(), 1)
        sb = out.메틸수.sum() / max(out.호출수.sum(), 1)
        P("  %-28s %-4s gene %.3f, spacer %.3f, difference %+.3f  (positions %d)"
          % (t.묶음, g, gb, sb, gb - sb, len(dd)))
LOG.close()
