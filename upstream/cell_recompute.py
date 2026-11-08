# -*- coding: utf-8 -*-
# 세포주 read×유닛 만들기 — cancer_recompute.py 와 같은 코드 논리 (구획 정의·ML 문턱·varlist 좌표 옮김 동일)
#   입력 : RU_{GM12878FS,GM12878,K562r}.tsv.gz  ($P11_DIR; as-ran value in ENV_RDNA.tsv)
# MAP_G8_CELLLINES.tsv (the cell-line rows of ~/p11_v11/AOG8_*/MAP_G8.tsv)
#   자체검정 : 내가 계산한 n_cpg / beta 가 RU 파일이 이미 적어 둔 n_cpg / beta 와 같은가
#                 (RU 는 파이프라인이 같은 지도 로 만든 것이므로, 같으면 틀이 맞은 것이다)
import io, os, re, sys, gzip, glob, collections, hashlib
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
P11 = __import__("os").environ["P11_DIR"] + "/"
CANON = __import__("os").environ["CANON_FA"]
LOG = io.open(os.path.join(HERE, "CELL_RECOMPUTE.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
MLCUT = 128
THR = 0.508
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()

canon = "".join(x.strip() for x in io.open(CANON, encoding="utf8") if x[0] != ">")
CPG = sorted(i + 1 for i in range(len(canon) - 1) if canon[i:i + 2] == "CG")
CPGSET = set(CPG)
SEG = [("alu", 1, 318), ("nts1", 319, 1178), ("ga", 1179, 1551), ("nts2", 1552, 2121), ("gene", 2122, 2240)]
UPS = (2032, 2121)
SEGOF = {p: nm for nm, a, b in SEG for p in CPG if a <= p <= b}
UPSET = set(p for p in CPG if UPS[0] <= p <= UPS[1])
SEGK = ["alu", "nts1", "ga", "nts2", "gene", "ups", "sp"]
P("정본 %d bp · CpG %d · 구획 %s · 프로모터 CpG %d"
  % (len(canon), len(CPG), " ".join("%s %d" % (n, sum(1 for p in CPG if a <= p <= b)) for n, a, b in SEG), len(UPSET)))

RUS = {"GM12878FS": P11 + "RU_GM12878FS.tsv.gz",
       "GM12878": P11 + "RU_GM12878.tsv.gz",
       "K562r": P11 + "RU_K562r.tsv.gz"}
MAPF = os.path.join(HERE, "MAP_G8_CELLLINES.tsv")
P()
P("[입력 md5]")
for k, v in list(RUS.items()) + [("MAP_G8_CELLLINES", MAPF)]:
    P("   %-18s %s  %s" % (k, md5(v), os.path.basename(v)))

# ---------- 지도 ----------
IN = re.compile(r"^(\d+)I(\d+)$"); DE = re.compile(r"^(\d+)D(\d+)$")
UNITN, UMAP = {}, {}
nbad = 0
mp = pd.read_csv(MAPF, sep="\t", dtype={"slot": str, "varlist": str})
for _, p in mp[mp.seg_type == "UNIT"].iterrows():
    smp = p["sample"]; sl = str(p["slot"]); ap = int(p["array_pos"]); ul = int(p["unit_len"])
    UNITN[(smp, sl)] = int(p["N"])
    vl = p["varlist"]; ins = {}; dele = {}
    if isinstance(vl, str) and vl and vl != "-":
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
    UMAP[(smp, sl, ap)] = arr
P()
P("[지도] UNIT %s · 길이 불일치 %d · 슬롯 %d" % (format(len(UMAP), ","), nbad, len(UNITN)))
P("   슬롯별 N : %s" % ", ".join("%s|%s=%d" % (k[0], k[1], v) for k, v in sorted(UNITN.items())))

# ---------- read x 유닛 ----------
rows = []; chk = []
for smp_file, f in RUS.items():
    nmap0 = nrow = 0
    with gzip.open(f, "rt", encoding="utf8") as h:
        ix = {c: i for i, c in enumerate(h.readline().rstrip("\n").split("\t"))}
        for line in h:
            p = line.rstrip("\n").split("\t")
            smp = p[ix["sample"]]; sl = p[ix["slot"]]; ap = int(p[ix["array_pos"]])
            cs = p[ix["calls"]]
            if not cs or cs == ".": continue
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
            nrow += 1
            N = UNITN.get((smp, sl), 0)
            r = dict(sample=smp, src=p[ix["data_source"]], slot=sl, ap=ap, read_id=p[ix["read_id"]],
                     N=N, n_cpg=tot, n_meth=met, b_call=met / tot, b_prob=sumv / (255.0 * tot),
                     ad=ap, prox=N - ap + 1)
            for k in SEGK:
                r["n_" + k] = segn[k]
                r["b_" + k] = (segm[k] / segn[k]) if segn[k] else np.nan
            rows.append(r)
            chk.append((smp, tot, int(p[ix["n_cpg"]]), met, int(p[ix["n_meth"]])))
    P("   %-22s 행 %s · 지도에 없어 버림 %d" % (os.path.basename(f), format(nrow, ","), nmap0))
D = pd.DataFrame(rows)

# ---------- 자체검정 ----------
P()
P("[자체검정] 내가 센 CpG·메틸 수 ↔ RU 파일이 적어 둔 값")
C = pd.DataFrame(chk, columns=["sample", "n_mine", "n_file", "m_mine", "m_file"])
for s, g in C.groupby("sample"):
    P("   %-12s 유닛 %5d · CpG 수 완전일치 %.2f%% · 메틸 수 완전일치 %.2f%% · |ΔCpG| 평균 %.3f"
      % (s, len(g), 100 * (g.n_mine == g.n_file).mean(), 100 * (g.m_mine == g.m_file).mean(),
         (g.n_mine - g.n_file).abs().mean()))
P("   이 값이 100%% 에 가까우면 RU 와 지도가 같은 틀 이라는 뜻이다.")

D.to_parquet(os.path.join(HERE, "CELL_READUNIT.parquet"), index=False)
P()
P("[요약] 세포주 read×유닛 %s · 시료 %s" % (format(len(D), ","), sorted(D["sample"].unique())))
D["state"] = np.where(D.b_call >= THR, "고메틸 유닛", "저메틸 유닛")
for s, g in D.groupby("sample"):
    P("   %-12s 유닛 %5d · read %4d · 전체 β %.4f · 저메틸유닛몫 %.3f · gene β %.4f · sp β %.4f · gene−sp %+.4f"
      % (s, len(g), g.read_id.nunique(), g.n_meth.sum() / g.n_cpg.sum(),
         (g.state == "저메틸 유닛").mean(),
         (g.b_gene * g.n_gene).sum() / g.n_gene.sum(),
         (g.b_sp * g.n_sp).sum() / g.n_sp.sum(),
         (g.b_gene * g.n_gene).sum() / g.n_gene.sum() - (g.b_sp * g.n_sp).sum() / g.n_sp.sum()))

# ---------- ML 척도 ----------
P()
P("[ML 척도] 호출 확률 분포 (0–255) — 판본·플랫폼이 다르면 β 가 생물학 없이 밀린다")
ML = {}
for f in RUS.values():
    with gzip.open(f, "rt", encoding="utf8") as h:
        ix = {c: i for i, c in enumerate(h.readline().rstrip("\n").split("\t"))}
        for line in h:
            p = line.rstrip("\n").split("\t")
            cs = p[ix["calls"]]
            if not cs or cs == ".": continue
            a = ML.setdefault(p[ix["sample"]], [0, 0, 0, 0])
            for it in cs.split(","):
                v = int(it.split(":")[1]); a[0] += 1; a[1] += v; a[2] += (v > 191); a[3] += (v < 64)
M = pd.DataFrame([dict(sample=k, 호출=v[0], ML평균=v[1] / v[0], 확신고=v[2] / v[0], 확신저=v[3] / v[0])
                  for k, v in ML.items()])
P(M.round(4).sort_values("sample").to_string(index=False))
M.to_csv(os.path.join(HERE, "CELL_ML.tsv"), sep="\t", index=False, float_format="%.5f")
P()
P("[저장] CELL_READUNIT.parquet · CELL_ML.tsv")
LOG.close()
