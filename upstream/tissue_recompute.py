# -*- coding: utf-8 -*-
# 조직 비교 자료 (HPRC Year7 같은 사람 10 명의 LCL ↔ PBMC) 를 같은 코드 논리로 read×유닛 으로
# Input : v11-add1_tissue_cancer 08_methylation/node3/FINAL/RU_{LCL,PBMC}.tsv.gz
#          같은 판본의 AOG8_LCL / AOG8_PBMC MAP_G8.tsv
import io, os, re, sys, gzip, glob, collections, time
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
A1 = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11-add1_tissue_cancer/"
CANON = __import__("os").environ["CANON_FA"]
LOG = io.open(os.path.join(HERE, "TISSUE_RECOMPUTE.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
MLCUT = 128

canon = "".join(x.strip() for x in io.open(CANON, encoding="utf8") if x[0] != ">")
CPG = sorted(i + 1 for i in range(len(canon) - 1) if canon[i:i + 2] == "CG")
CPGSET = set(CPG)
SEG = [("alu", 1, 318), ("nts1", 319, 1178), ("ga", 1179, 1551), ("nts2", 1552, 2121), ("gene", 2122, 2240)]
UPS = (2032, 2121)
SEGOF = {p: nm for nm, a, b in SEG for p in CPG if a <= p <= b}
UPSET = set(p for p in CPG if UPS[0] <= p <= UPS[1])
SEGK = ["alu", "nts1", "ga", "nts2", "gene", "ups", "sp"]

# ---------- 대상 시료 : RU 안에 있는 것 전부 ----------
WANT = set()
SRCOF = {}
for f in (A1 + "08_methylation/node3/FINAL/RU_LCL.tsv.gz", A1 + "08_methylation/node3/FINAL/RU_PBMC.tsv.gz"):
    with gzip.open(f, "rt", encoding="utf8") as h:
        ix = {c: i for i, c in enumerate(h.readline().rstrip("\n").split("\t"))}
        for line in h:
            p = line.rstrip("\n").split("\t")
            WANT.add(p[ix["sample"]]); SRCOF[p[ix["sample"]]] = p[ix["data_source"]]
P("대상 시료 %d : %s" % (len(WANT), sorted(WANT)))

# ---------- 지도 ----------
IN = re.compile(r"^(\d+)I(\d+)$"); DE = re.compile(r"^(\d+)D(\d+)$")
UNITN, UMAP = {}, {}
for mp in (A1 + "03_maps/AOG8_inputs_merged/node3/AOG8_LCL/MAP_G8.tsv",
           A1 + "03_maps/AOG8_inputs_merged/node3/AOG8_PBMC/MAP_G8.tsv"):
    with io.open(mp, encoding="utf8") as h:
        ix = {c: i for i, c in enumerate(h.readline().rstrip("\n").split("\t"))}
        for line in h:
            p = line.rstrip("\n").split("\t")
            if p[ix["seg_type"]] != "UNIT": continue
            smp = p[ix["sample"]]
            if smp not in WANT: continue
            sl = p[ix["slot"]]; ap = int(p[ix["array_pos"]]); ul = int(p[ix["unit_len"]])
            if (smp, sl, ap) in UMAP: continue
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
            UMAP[(smp, sl, ap)] = arr
P("지도 UNIT %s · 슬롯 %d" % (format(len(UMAP), ","), len(UNITN)))

rows = []
nhit = ntot = 0
for f in (A1 + "08_methylation/node3/FINAL/RU_LCL.tsv.gz", A1 + "08_methylation/node3/FINAL/RU_PBMC.tsv.gz"):
    t0 = time.time(); nf = nmap0 = 0
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
                ntot += 1
                if a < 1 or a > L: continue
                cp = int(arr[a])
                if cp == 0 or cp not in CPGSET: continue
                nhit += 1
                m = 1 if v >= MLCUT else 0
                tot += 1; met += m; sumv += v
                sg = SEGOF[cp]; segn[sg] += 1; segm[sg] += m
                if cp in UPSET: segn["ups"] += 1; segm["ups"] += m
                if sg != "gene": segn["sp"] += 1; segm["sp"] += m
            if tot == 0: continue
            nf += 1
            N = UNITN.get((smp, sl), 0)
            r = dict(sample=smp + ("_LCL" if "LCL" in p[ix["data_source"]] else "_PBMC"),
                     person=smp, tissue=("LCL" if "LCL" in p[ix["data_source"]] else "PBMC"),
                     src=p[ix["data_source"]], slot=sl, ap=ap, read_id=p[ix["read_id"]],
                     N=N, n_cpg=tot, n_meth=met, b_call=met / tot, b_prob=sumv / (255.0 * tot))   # 자 대조 생산자 — 두 열을 함께 낸다(복원)
            for k in SEGK:
                r["n_" + k] = segn[k]
                r["b_" + k] = (segm[k] / segn[k]) if segn[k] else np.nan
            rows.append(r)
    P("   %-18s 행 %s · 지도없음 %d · %.0f 초" % (os.path.basename(f), format(nf, ","), nmap0, time.time() - t0))
D = pd.DataFrame(rows)
P()
P("[틀 검정] 호출이 정본 CpG 에 붙는 비율 %.4f  (기준선 0.967–0.979)" % (nhit / ntot))
D.to_parquet(os.path.join(HERE, "TISSUE_READUNIT.parquet"), index=False)
P("[저장] TISSUE_READUNIT.parquet · read×유닛 %s" % format(len(D), ","))
P()
P("[요약] 사람마다 두 조직")
for pr, g in D.groupby("person"):
    a = g[g.tissue == "LCL"]; b = g[g.tissue == "PBMC"]
    def wb(x, k):
        n = x["n_" + k].sum()
        return ((x["b_" + k] * x["n_" + k]).sum() / n) if n else np.nan
    P("   %-9s LCL 유닛 %5d β %.4f gene-sp %+.4f | PBMC 유닛 %5d β %.4f gene-sp %+.4f"
      % (pr, len(a), a.n_meth.sum() / a.n_cpg.sum(), wb(a, "gene") - wb(a, "sp"),
         len(b), b.n_meth.sum() / b.n_cpg.sum(), wb(b, "gene") - wb(b, "sp")))
LOG.close()
