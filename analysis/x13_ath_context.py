# -*- coding: utf-8 -*-
# 애기장대 5S : CG / CHG / CHH 세 문맥을 read x 유닛으로 · 유닛 전체 / coding(5S 유전자) 따로.
#   (동물 BAM 은 5mCG 모델이라 CHG/CHH 호출이 없다 — 애기장대만 dorado all-context)
#   BAM 읽기·유닛·유전자 좌표는 x11_coding_readunit.py 함수를 그대로 쓴다(exec · 로그 줄만 바꿈).
#   문맥 판정은 units_matrix_g8.py 와 같다 : 정방향 C → 다음 G=CG · 아니고 다다음 G=CHG · 그 외 CHH
#                                         역방향(참조 G) → 앞 C=CG · 아니고 앞앞 C=CHG · 그 외 CHH
#   MM 이 `C+m.` 이므로 생략된 C 는 비메틸 호출로 센다(x11 과 같음 · 안 세면 기존표와 안 맞음).
#   검산 : 문맥별 유닛 전체 (n_calls, n_5mC) 가 기존 read_unit.tsv 와 같아야 한다.
#   판정 : x10 규칙(두 번째 봉우리 >= 0.30 떨어짐 · 골 <= 그 봉우리 절반 · 그쪽 분자 >= 2%) + 이항 귀무 과분산비
#   조건 : 호출 수 >= 그 array·문맥·구간 중앙값의 절반 (최소 3)
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
                       "see README.md, section 'How to run it')" % k)
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
P("애기장대 BAM 기록 %d · 통계 %s · 행 %d" % (STAT["rec"], dict(STAT), len(D)))

# ---- 검산 : 기존 read_unit.tsv 문맥별
old = pd.read_csv(glob.glob(B + "/03_reads/%s/*.read_unit.tsv*" % SP)[0], sep="\t")
for c in ["CG", "CHG", "CHH"]:
    o = old[old.context == c]
    n = D[D.context == c]
    M = o.merge(n, on=["read_id", "array_id", "array_pos"], how="inner", suffixes=("_old", ""))
    same = ((M.n_calls_old == M.n_calls) & (M.n_5mC_old == M.n_5mC)).mean()
    near = ((M.n_calls_old - M.n_calls).abs() <= 1).mean()
    P("  검산 %-3s : 옛 %d · 새 %d · 공통 %d · 완전일치 %.4f · 호출수 차 <=1 %.4f · 옛 평균 β %.3f / 새 %.3f"
      % (c, len(o), len(n), len(M), same, near, (M.n_5mC_old / M.n_calls_old).mean(), (M.n_5mC / M.n_calls).mean()))
P()

# ---- array x 문맥 x 구간 판정
ST = pd.read_csv(OUT + "/tables/OTHERSP_ARRAY_STRUCT.tsv", sep="\t")
UL = dict(zip(ST.array_id, ST.unit_len))
NU = dict(zip(ST.array_id, ST.n_unit))
D["n_sp"] = D.n_calls - D.n_calls_gene
ROWS, PANEL = [], {}
for (aid, c), g in D.groupby(["array_id", "context"]):
    for part, ncol_, mcol_ in [("전체", "n_calls", "n_5mC"), ("coding", "n_calls_gene", "n_5mC_gene")]:
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
P("  array                          유닛x개수  문맥  구간     read-유닛  호출중앙  평균   sd    10%%-90%%     과분산  봉1  봉2  골비  소수쪽  판정")
order = T.groupby("array_id").n.max().sort_values(ascending=False).index.tolist()
for aid in order:
    for r in T[T.array_id == aid].sort_values(["context", "part"], ascending=[True, False]).itertuples(index=False):
        P("  %-30s %3.0fx%-5.0f %-4s %-7s %8s  %5.0f   %.3f  %.3f  %.2f-%.2f   %5.1f  %.2f %s %s %.3f  %s"
          % (aid.split("__")[-1], r.unit_len, r.n_unit, r.context, r.part, format(r.n, ","), r.calls, getattr(r, "mean"),
             r.sd, r.q10, r.q90, r.od, r.peak1, ("%.2f" % r.peak2) if r.peak2 == r.peak2 else "  - ",
             ("%.2f" % r.trough_ratio) if r.trough_ratio == r.trough_ratio else "  - ", r.minor_mass, r.call))
P()
P("  판정 집계")
for (c, part), gg in T.groupby(["context", "part"]):
    P("   %-4s %-7s : %s" % (c, part, " · ".join("%s %d" % (k, v) for k, v in gg.call.value_counts().items())))

# ---- 그림 : 구간(전체/coding)마다 한 장 · 행 = array · 열 = 문맥
COLC = {"CG": "#7e57c2", "CHG": "#26a69a", "CHH": "#ef6c00"}
for part in ["전체", "coding"]:
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
                 % ("유닛 전체" if part == "전체" else "coding(5S 유전자)만"), fontsize=12)
    fig.tight_layout(rect=[0, 0.01, 1, 0.975])
    fig.savefig(OUT + "/figures/Fig_ath_context_%s.png" % part)
    plt.close(fig)
    P("  저장 : figures/Fig_ath_context_%s.png" % part)
LOG.close()
