# -*- coding: utf-8 -*-
# — 보노보를 guppy 6.3.8 read 만으로 다시 계산 (+ 6.3.7 만 · 혼합)
#   read 판본 대리 표지 = BAM MM 태그 머리 : "C+m?..." 만 = 6.3.8 · "C+h?...;C+m?..." = 6.3.7
# counted again here
# Metric, with the original definition
#     A1 array − flank : AVF 정의(array_vs_flank.py) — MM 의 모든 m 호출 · 정렬 일치 자리 · ML >= 128 ·
#                        구간 = AVF.tsv 의 reg_start..reg_end (1기반 양끝 포함). 혼합판이 AVF n_c · n_5mC 를 재현해야 함.
#     read_unit 층 (03_reads 보노보 read_unit · CG · x10 거름 : n_calls >= array 중앙값 절반)
#       A7 평균 β · OFF(<0.35) · x10 이봉 판정 · A8 과분산비(x8 정의)
#       A4 끝 섞임비율(§5-C 인라인 계산을 코드로 옮김) : 끝거리 = min(pos, N−pos+1) · 끝 1–5 vs 안쪽 21+ OFF 비율 · ON(β>0.65) 평균
#       A3 끝 − 안쪽(read_unit 판) : 유닛별 호출가중 β → 끝 1–5 평균 − 안쪽 20+ 평균
#          (원 A3 값 −0.474 는 자리별 bedMethyl 합산판 x5/x6 — 판본별로는 bedMethyl 을 못 나눠 read_unit 판으로 비교)
#     A9 coding : OTHERSP_READUNIT_CODING.tsv.gz (x11) · x12 거름 · coding OFF · 판정
#   부수 : §5-C 끝 섞임비율 표(유인원 6 종 · 혼합)를 코드 파일로 재현
#   출력 : tables/BONOBO_BY_BASECALLER.tsv · tables/APE_END_MIXING.tsv · X17_BONOBO_BY_BASECALLER.txt
import io, struct, collections
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


OUT = _need("RDNA_PREV_OUT_DIR")
B = _need("RDNA_DATA_ROOT")
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
call10 = G10["call"]
LOG = io.open(OUT + "/X17_BONOBO_BY_BASECALLER.txt", "w", encoding="utf8", newline="\n")
rng = np.random.RandomState(20260922)


def P(s=""):
    print(s)
    LOG.write(s + "\n")
    LOG.flush()


SP = "bonobo_mPanPan1"
BAM = B + "/04_bam/%s/%s.5S.bam" % (SP, SP)
AVF = pd.read_csv(B + "/09_array_vs_flank/AVF.tsv", sep="\t")
AVF = AVF[AVF.species_key == SP]
REG = [(r.region, r.chrom, int(r.reg_start), int(r.reg_end), int(r.n_c), int(r.n_5mC)) for r in AVF.itertuples()]
NIBCOMP = G["NIBCOMP"]

# ---------------------------------------------------------------- BAM : read 판본 + AVF 구간 호출
sig = {}
cnt = collections.defaultdict(lambda: np.zeros(2, dtype=np.int64))   # (grp, region) -> n, m
for refs, r in G["bam_records"](BAM):
    refid, pos0, lrn, mapq, _bin, ncig, flag, lseq = struct.unpack("<iiBBHHHi", r[:20])
    if refid < 0 or flag & 0x904:
        continue
    rid = r[32:32 + lrn - 1].decode()
    p = 32 + lrn
    cig = np.frombuffer(r[p:p + 4 * ncig], dtype="<u4"); p += 4 * ncig
    sq = r[p:p + (lseq + 1) // 2]; p += (lseq + 1) // 2 + lseq
    T = G["parse_tags"](r, p)
    if "MM" not in T or "ML" not in T:
        continue
    head = T["MM"].split(",")[0]
    grp = "6.3.7" if T["MM"].startswith("C+h") else ("6.3.8" if T["MM"].startswith("C+m") else "기타:" + head)
    sig[rid] = grp
    ops, lens = cig & 15, (cig >> 4).astype(np.int64)
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
    best = {}
    for ent in T["MM"].rstrip(";").split(";"):
        if not ent:
            continue
        parts = ent.split(",")
        hd = parts[0]
        if hd[-1] in "?.":
            hd = hd[:-1]
        base, st_, mc = hd[0], hd[1], hd[2:]
        d = np.array([int(x) for x in parts[1:]], dtype=np.int64)
        nc = len(mc)
        vals = ml[mi:mi + nc * len(d)].reshape(-1, nc) if len(d) else np.zeros((0, nc), np.uint8)
        mi += nc * len(d)
        if base != "C" or st_ != "+" or not len(d):
            continue
        k = np.cumsum(d + 1) - 1
        idx = cpos[k]
        for j, c in enumerate(mc):
            if c == "m":
                for q_, v in zip(idx.tolist(), vals[:, j].tolist()):
                    best[q_] = max(best.get(q_, 0), v)
            elif c == "h":
                for q_ in idx.tolist():
                    best.setdefault(q_, 0)          # AVF 는 h 만 있는 자리도 호출 1 개로 센다(pm=0)
    if not best:
        continue
    qi = np.fromiter(best.keys(), dtype=np.int64); mv = np.fromiter(best.values(), dtype=np.int64)
    sidx = (lseq - 1 - qi) if rev else qi
    rp1 = q2r[sidx] + 1
    ok = q2r[sidx] >= 0
    chrom = refs[refid]
    for reg, c, lo, hi, _, _ in REG:
        if c != chrom:
            continue
        w = ok & (rp1 >= lo) & (rp1 <= hi)
        if w.any():
            for gname in (grp, "혼합"):
                cnt[(gname, reg)] += np.array([w.sum(), (mv[w] >= 128).sum()])

P("보노보 BAM read 판본 : " + " · ".join("%s %d" % kv for kv in collections.Counter(sig.values()).items()))
P()
P("A1 검산 — 혼합판이 AVF.tsv 를 재현하나")
for reg, c, lo, hi, n_c, n_m in REG:
    n, m = cnt[("혼합", reg)]
    P("   %-18s AVF n_c %9d n_5mC %9d | BAM %9d %9d | %s" % (reg, n_c, n_m, n, m, "같음" if (n == n_c and m == n_m) else "다름"))
ROWS = []
for gname in ("혼합", "6.3.8", "6.3.7"):
    a = cnt[(gname, "array")]; fl = sum((cnt[(gname, r)] for r in ("flank_left_100kb", "flank_right_100kb")), np.zeros(2, dtype=np.int64))
    ROWS.append(dict(grp=gname, metric="A1 array β (AVF 정의)", value=a[1] / a[0], n=int(a[0])))
    ROWS.append(dict(grp=gname, metric="A1 flank β (L+R)", value=fl[1] / fl[0], n=int(fl[0])))
    ROWS.append(dict(grp=gname, metric="A1 array − flank", value=a[1] / a[0] - fl[1] / fl[0], n=int(a[0])))

# ---------------------------------------------------------------- read_unit 층
RU = pd.read_csv(B + "/03_reads/%s/%s__%s.read_unit.tsv" % (SP, SP, SP), sep="\t")
RU = RU[RU.context == "CG"].drop_duplicates(["read_id", "array_id", "array_pos"]).copy()
RU["grp"] = RU.read_id.map(sig).fillna("BAM에없음")
P()
P("read_unit read 판본 : " + " · ".join("%s %d read" % (k, v) for k, v in RU.groupby("grp").read_id.nunique().items()))
allcall = RU.copy()
med = RU.n_calls.median()
RU = RU[RU.n_calls >= 0.5 * med].copy()
RU["b"] = RU.n_5mC / RU.n_calls
N = RU.array_pos.max()
RU["p_end"] = np.minimum(RU.array_pos, N - RU.array_pos + 1)
CD = pd.read_csv(OUT + "/tables/OTHERSP_READUNIT_CODING.tsv.gz", sep="\t")
CD = CD[(CD.species_key == SP) & (CD.n_gene_hits > 0)].drop_duplicates(["read_id", "array_id", "array_pos"]).copy()
CD["grp"] = CD.read_id.map(sig).fillna("BAM에없음")
gmed = CD.n_calls_gene.median()
CD = CD[CD.n_calls_gene >= max(3, 0.5 * gmed)].copy()
CD["bg"] = CD.n_5mC_gene / CD.n_calls_gene


def od_of(g):
    p = g.groupby("array_pos").b.transform("mean").values
    n = g.n_calls.values.astype(int)
    sim = rng.binomial(n, np.clip(p, 0, 1)) / n
    return np.var(g.b.values - p) / max(np.var(sim - p), 1e-12)


for gname in ("혼합", "6.3.8", "6.3.7"):
    g = RU if gname == "혼합" else RU[RU.grp == gname]
    a = allcall if gname == "혼합" else allcall[allcall.grp == gname]
    ROWS.append(dict(grp=gname, metric="array 호출가중 β (read_unit · 거름 없음 · CMD01 정의)", value=a.n_5mC.sum() / a.n_calls.sum(), n=int(a.read_id.nunique())))
    r, _, _ = call10(g.b.values, float(g.n_calls.median()))
    ROWS.append(dict(grp=gname, metric="A7 read-유닛 평균 β", value=g.b.mean(), n=len(g)))
    ROWS.append(dict(grp=gname, metric="A7 OFF(<0.35)", value=(g.b < 0.35).mean(), n=len(g)))
    ROWS.append(dict(grp=gname, metric="A7 x10 판정 · 봉우리", value=np.nan, n=len(g), note="%s · %.2f/%s · 골비 %s" % (r["call"], r["peak1"], ("%.2f" % r["peak2"]) if r["peak2"] == r["peak2"] else "-", ("%.2f" % r["trough_ratio"]) if r["trough_ratio"] == r["trough_ratio"] else "-")))
    ROWS.append(dict(grp=gname, metric="A8 과분산비(x8 정의)", value=od_of(g), n=len(g)))
    e, i = g[g.p_end <= 5], g[g.p_end >= 21]
    ROWS.append(dict(grp=gname, metric="A4 끝1–5 OFF", value=(e.b < 0.35).mean(), n=len(e)))
    ROWS.append(dict(grp=gname, metric="A4 안쪽21+ OFF", value=(i.b < 0.35).mean(), n=len(i)))
    ROWS.append(dict(grp=gname, metric="A4 ON(β>0.65) 평균 끝 / 안쪽", value=np.nan, n=len(e), note="%.3f / %.3f" % (e.b[e.b > 0.65].mean(), i.b[i.b > 0.65].mean())))
    u = g.groupby("array_pos").agg(m=("n_5mC", "sum"), n=("n_calls", "sum"), pe=("p_end", "first"))
    u["b"] = u.m / u.n
    ROWS.append(dict(grp=gname, metric="A3 끝1–5 − 안쪽20+ (read_unit 판 · 유닛 호출가중)", value=u[u.pe <= 5].b.mean() - u[u.pe >= 20].b.mean(), n=len(u)))
    c = CD if gname == "혼합" else CD[CD.grp == gname]
    rc, _, _ = call10(c.bg.values, float(c.n_calls_gene.median()))
    ROWS.append(dict(grp=gname, metric="A9 coding OFF(<0.35)", value=(c.bg < 0.35).mean(), n=len(c), note=rc["call"]))

T = pd.DataFrame(ROWS)
T.to_csv(OUT + "/tables/BONOBO_BY_BASECALLER.tsv", sep="\t", index=False, float_format="%.5f")
P()
P("%-52s %10s %10s %10s" % ("지표", "혼합", "6.3.8 만", "6.3.7 만"))
for m in T.metric.unique():
    v = []
    for gname in ("혼합", "6.3.8", "6.3.7"):
        rr = T[(T.metric == m) & (T.grp == gname)].iloc[0]
        v.append(rr.note if isinstance(rr.get("note"), str) and rr.value != rr.value else ("%+.4f" % rr.value if "−" in m or "- " in m else "%.4f" % rr.value))
    P("%-52s %10s %10s %10s" % (m, *v))

# ---------------------------------------------------------------- §5-C 끝 섞임비율 (유인원 6 종 · 혼합) 코드판
NM = {"bonobo_mPanPan1": "보노보", "chimp_mPanTro3": "침팬지", "gorilla_mGorGor1": "고릴라",
      "orang_bornean_mPonPyg2": "보르네오오랑", "orang_sumatran_mPonAbe1": "수마트라오랑", "siamang_mSymSyn1": "시아망"}
D = pd.concat([pd.read_csv(B + "/03_reads/%s/%s__%s.read_unit.tsv" % (k, k, k), sep="\t") for k in NM])
D = D[D.context == "CG"]
m_ = D.groupby("array_id").n_calls.transform("median"); D = D[D.n_calls >= 0.5 * m_].copy()
D["b"] = D.n_5mC / D.n_calls
Nn = D.groupby("array_id").array_pos.transform("max")
D["p_end"] = np.minimum(D.array_pos, Nn - D.array_pos + 1)
E = []
P()
P("§5-C 끝 섞임비율 (코드판 · 원 인라인 계산과 같은 정의)")
for k, nm in NM.items():
    g = D[D.species_key == k]
    e, i = g[g.p_end <= 5], g[g.p_end >= 21]
    lo_e, lo_i = (e.b < 0.35).mean(), (i.b < 0.35).mean()
    E.append(dict(species=nm, off_end1_5=lo_e, off_in21=lo_i, fold=lo_e / max(lo_i, 1e-9),
                  on_beta_end=e.b[e.b > 0.65].mean(), on_beta_in=i.b[i.b > 0.65].mean()))
    P("   %-8s 끝1-5 %.3f vs 안쪽21+ %.3f (%.1f배) · 고메틸 β 끝 %.3f vs 안쪽 %.3f" % (nm, lo_e, lo_i, lo_e / max(lo_i, 1e-9), E[-1]["on_beta_end"], E[-1]["on_beta_in"]))
pd.DataFrame(E).to_csv(OUT + "/tables/APE_END_MIXING.tsv", sep="\t", index=False, float_format="%.5f")
LOG.close()
