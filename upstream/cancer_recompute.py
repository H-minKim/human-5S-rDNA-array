# -*- coding: utf-8 -*-
# 암 시료 4 가지 직접 재계산 ( · 이전 표를 보지 않고 원자료에서)
# Input : v11-add1_tissue_cancer\08_methylation\node3\FINAL\RU_{HG008,CASTLE,CASTLE1,COLOr}.tsv.gz (read x unit, calls column)
#                v11-add1_tissue_cancer\03_maps\AOG8_inputs_merged\node3\AOG8_*\MAP_G8.tsv               (varlist · N)
#                정본 유닛 서열 canon_v2.fa
#   HG008 는 추가1(GIAB_HG008_ONT · 애초에 CPGONLY)을 쓴다 — v11 코호트의 GIAB_SOMATIC_ONT(모든-C)가 아니다
#   방법은 v100_build_ru_f4.py 와 같은 코드 논리 (구획 정의·ML 문턱·varlist 좌표 옮김)
#   산출 4 가지
#     ① 전체 메틸 정도        CR_OVERALL.tsv
#     ② 구획별(스페이서·유전자부 등) CR_SEG.tsv
#     ③ 근위/원위 경계         CR_BOUNDARY.tsv
#     ④ 하플로 안 read 분포     CR_READ.tsv  (+ 요약 CR_READDIST.tsv)
import io, os, re, sys, gzip, glob, collections, time
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
A1 = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11-add1_tissue_cancer/"
CANON = __import__("os").environ["CANON_FA"]
LOG = io.open(os.path.join(HERE, "CANCER_RECOMPUTE.txt"), "w", encoding="utf8", newline="\n")
def P(s=""):
    print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
MLCUT = 128
THR = 0.508 # unit on/off threshold, the canonical trough

# ---------- 정본 CpG · 구획 ----------
canon = "".join(x.strip() for x in io.open(CANON, encoding="utf8") if x[0] != ">")
CPG = sorted(i + 1 for i in range(len(canon) - 1) if canon[i:i + 2] == "CG")
CPGSET = set(CPG)
SEG = [("alu", 1, 318), ("nts1", 319, 1178), ("ga", 1179, 1551), ("nts2", 1552, 2121), ("gene", 2122, 2240)]
UPS = (2032, 2121)
SEGOF = {p: nm for nm, a, b in SEG for p in CPG if a <= p <= b}
UPSET = {p for p in CPG if UPS[0] <= p <= UPS[1]}
P("정본 길이 %d · CpG %d · 구획 %s · 프로모터 CpG %d"
  % (len(canon), len(CPG), " ".join("%s %d" % (n, sum(1 for p in CPG if a <= p <= b)) for n, a, b in SEG), len(UPSET)))

# ---------- 대상 시료 ----------
TUM = ["HG008-T", "COLO829-Tr", "COLO829-T", "H1437-T", "H2009-T", "HCC1395-T", "HCC1937-T", "HCC1954-T", "Hs578T-T"]
NOR = ["HG008-N-P", "HG008-N-D", "COLO829-BLr", "COLO829-BL", "H1437-N", "H2009-N",
       "HCC1395-N", "HCC1937-N", "HCC1954-N", "Hs578T-N"]
WANT = set(TUM) | set(NOR)

# ---------- 지도 ----------
IN = re.compile(r"^(\d+)I(\d+)$"); DE = re.compile(r"^(\d+)D(\d+)$")
UNITN, UMAP = {}, {}
t0 = time.time(); nu = nbad = 0
for mp in sorted(glob.glob(A1 + "03_maps/AOG8_inputs_merged/node3/AOG8_*/MAP_G8.tsv")):
    with io.open(mp, encoding="utf8") as h:
        hd = h.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hd)}
        for line in h:
            p = line.rstrip("\n").split("\t")
            if p[ix["seg_type"]] != "UNIT": continue
            smp = p[ix["sample"]]
            if smp not in WANT: continue
            sl = p[ix["slot"]]; ap = int(p[ix["array_pos"]]); ul = int(p[ix["unit_len"]])
            if (smp, sl, ap) in UMAP: continue          # 먼저 읽은 그룹 것을 쓴다
            UNITN[(smp, sl)] = int(p[ix["N"]])
            vl = p[ix["varlist"]]; ins = {}; dele = {}
            if vl and vl != "-":
                for v in vl.split(","):
                    m = IN.match(v)
                    if m: ins[int(m.group(1))] = int(m.group(2)); continue
                    m = DE.match(v)
                    if m: dele[int(m.group(1))] = int(m.group(2))
            arr = np.zeros(ul + 2, dtype=np.int16)
            sp = 1; cp = 1; skip = 0
            while cp <= 2240 and sp <= ul:
                if skip > 0: skip -= 1; cp += 1; continue
                if cp in dele: skip = dele[cp] - 1; cp += 1; continue
                arr[sp] = cp; sp += 1
                if cp in ins:
                    for _ in range(ins[cp]):
                        if sp <= ul: arr[sp] = 0; sp += 1
                cp += 1
            if sp - 1 != ul: nbad += 1
            UMAP[(smp, sl, ap)] = arr; nu += 1
P("지도 UNIT %s · 길이 불일치 %d · 슬롯 %d · %.0f 초"
  % (format(nu, ","), nbad, len(UNITN), time.time() - t0))
P("  슬롯별 N : %s" % ", ".join("%s|%s=%d" % (k[0], k[1], v) for k, v in sorted(UNITN.items())))

# ---------- read x 유닛 훑기 ----------
SEGK = ["alu", "nts1", "ga", "nts2", "gene", "ups", "sp"]
rows = []
nline = nmap0 = nskip = 0
for f in sorted(glob.glob(A1 + "08_methylation/node3/FINAL/RU_*.tsv.gz")):
    t1 = time.time(); nf = 0
    with gzip.open(f, "rt", encoding="utf8") as h:
        hd = h.readline().rstrip("\n").split("\t"); ix = {c: i for i, c in enumerate(hd)}
        for line in h:
            p = line.rstrip("\n").split("\t")
            smp = p[ix["sample"]]
            if smp not in WANT: continue
            sl = p[ix["slot"]]; ap = int(p[ix["array_pos"]]); rid = p[ix["read_id"]]
            cs = p[ix["calls"]]
            if not cs or cs == ".": nskip += 1; continue
            arr = UMAP.get((smp, sl, ap))
            if arr is None: nmap0 += 1; continue
            L = len(arr) - 1
            tot = met = 0; sumv = 0
            segn = collections.Counter(); segm = collections.Counter()
            for it in cs.split(","):
                a, v = it.split(":"); a = int(a); v = int(v)
                if a < 1 or a > L: continue
                cp = int(arr[a])
                if cp == 0 or cp not in CPGSET: continue
                m = 1 if v >= MLCUT else 0
                tot += 1; met += m; sumv += v
                sg = SEGOF[cp]; segn[sg] += 1; segm[sg] += m
                if cp in UPSET: segn["ups"] += 1; segm["ups"] += m
                if sg != "gene": segn["sp"] += 1; segm["sp"] += m
            if tot == 0: continue
            nline += 1; nf += 1
            r = dict(sample=smp, src=p[ix["data_source"]], slot=sl, ap=ap, read_id=rid,
                     N=UNITN.get((smp, sl), 0), n_cpg=tot, n_meth=met,
                     b_call=met / tot, b_prob=sumv / (255.0 * tot))
            for k in SEGK:
                r["n_" + k] = segn[k]
                r["b_" + k] = (segm[k] / segn[k]) if segn[k] else np.nan
            rows.append(r)
    P("  %-18s 대상 행 %s · %.0f 초" % (os.path.basename(f), format(nf, ","), time.time() - t1))
D = pd.DataFrame(rows)
P()
P("read x 유닛 행 %s · 지도없음 %d · calls 빈 행 %d" % (format(len(D), ","), nmap0, nskip))
D["role"] = np.where(D["sample"].isin(TUM), "종양", "정상")
D["hap"] = D["sample"] + "|" + D.slot.astype(str)
D["ad"] = D.N - D.ap + 1                      # 원위 끝에서 센 번호(지도틀 array_pos = ap 는 원위 index 가 아니라 지도 번호)
D.to_parquet(os.path.join(HERE, "CR_READUNIT.parquet"), index=False)

# ---------- ① 전체 메틸 정도 ----------
O = D.groupby(["sample", "role", "slot", "src"]).apply(lambda x: pd.Series({
    "N": int(x.N.iloc[0]), "유닛": x.ap.nunique(), "read": x.read_id.nunique(),
    "호출": int(x.n_cpg.sum()),
    "β_call": x.n_meth.sum() / x.n_cpg.sum(),
    "β_prob": np.average(x.b_prob, weights=x.n_cpg),
    "유닛β중앙": x.groupby("ap").apply(lambda y: y.n_meth.sum() / y.n_cpg.sum()).median(),
    "켜진유닛비율": (x.groupby("ap").apply(lambda y: y.n_meth.sum() / y.n_cpg.sum()) >= THR).mean(),
})).reset_index()
P()
P("[①] 전체 메틸 정도 (시료·슬롯)")
P(O.to_string(index=False))
O.to_csv(os.path.join(HERE, "CR_OVERALL.tsv"), sep="\t", index=False, float_format="%.4f")

# ---------- ② 구획별 ----------
srows = []
for (smp, role, sl), g in D.groupby(["sample", "role", "slot"]):
    r = dict(sample=smp, role=role, slot=sl)
    for k in SEGK:
        n = g["n_" + k].sum()
        r[k] = (g["b_" + k] * g["n_" + k]).sum() / n if n else np.nan
        r["n_" + k] = int(n)
    srows.append(r)
S = pd.DataFrame(srows)
S["gene_minus_sp"] = S["gene"] - S["sp"]
P()
P("[②] 구획별 β (호출 가중) — sp = 유전자부 밖 전부(스페이서) · ups = 프로모터(NTS2 끝 90bp)")
P(S[["sample", "role", "slot", "alu", "nts1", "ga", "nts2", "gene", "ups", "sp", "gene_minus_sp"]]
  .round(4).to_string(index=False))
S.to_csv(os.path.join(HERE, "CR_SEG.tsv"), sep="\t", index=False, float_format="%.4f")

# ---------- ③ 근위/원위 경계 ----------
def boundary(u, side):
    """u : ap -> 유닛 β (오름차순 ap). side='ad' 면 원위(ap 작은 쪽)에서, 'ap' 면 근위(ap 큰 쪽)에서
       경계 = 끝에서부터 세어 처음으로 β >= THR 가 되는 자리 번호 (1 부터)."""
    s = u.sort_index()
    seq = s.values if side == "ad" else s.values[::-1]
    for i, v in enumerate(seq, 1):
        if v >= THR: return i
    return np.nan
brows = []
for (smp, role, sl), g in D.groupby(["sample", "role", "slot"]):
    u = g.groupby("ap").apply(lambda y: y.n_meth.sum() / y.n_cpg.sum())
    cov = g.groupby("ap").read_id.nunique()
    u5 = u[cov >= 3]                                   # read 3 개 이상인 자리만
    brows.append(dict(sample=smp, role=role, slot=sl, N=int(g.N.iloc[0]), 자리=len(u), 자리_read3=len(u5),
                      원위경계=boundary(u5, "ad"), 근위경계=boundary(u5, "ap"),
                      원위끝5유닛β=u5.head(5).mean(), 근위끝5유닛β=u5.tail(5).mean(),
                      안쪽β=u5.iloc[10:-10].mean() if len(u5) > 25 else np.nan))
B = pd.DataFrame(brows)
P()
P("[③] 경계 (유닛 β 가 %.3f 를 처음 넘는 자리 · read 3 개 이상 자리만)" % THR)
P(B.round(4).to_string(index=False))
B.to_csv(os.path.join(HERE, "CR_BOUNDARY.tsv"), sep="\t", index=False, float_format="%.4f")

# ---------- ④ 하플로 안 read 분포 ----------
R = D.groupby(["sample", "role", "slot", "read_id"]).apply(lambda x: pd.Series({
    "유닛": x.ap.nunique(), "호출": int(x.n_cpg.sum()),
    "β": x.n_meth.sum() / x.n_cpg.sum()})).reset_index()
R = R[R.호출 >= 50]
RD = R.groupby(["sample", "role", "slot"]).apply(lambda x: pd.Series({
    "read": len(x), "β중앙": x.β.median(), "β_sd": x.β.std(),
    "10%": x.β.quantile(.1), "90%": x.β.quantile(.9),
    "저(β<0.2)": (x.β < .2).mean(), "고(β>0.8)": (x.β > .8).mean(),
    "가운데(0.2-0.8)": ((x.β >= .2) & (x.β <= .8)).mean()})).reset_index()
P()
P("[④] 하플로 안 read 분포 (read 당 호출 50 개 이상)")
P(RD.round(4).to_string(index=False))
R.to_csv(os.path.join(HERE, "CR_READ.tsv"), sep="\t", index=False, float_format="%.4f")
RD.to_csv(os.path.join(HERE, "CR_READDIST.tsv"), sep="\t", index=False, float_format="%.4f")
P()
P("[저장] CR_READUNIT.parquet · CR_OVERALL · CR_SEG · CR_BOUNDARY · CR_READ · CR_READDIST")
LOG.close()
