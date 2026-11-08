# -*- coding: utf-8 -*-
# v143 (2026-09-25) junction 창의 서열 주석 — v142 와 같은 bp 창 에 CpG 밀도와 Alu 유사 구간을 센다
#   참조 = rc(flank_final2 fasta) if rcf=Y else fasta  (v116b 와 같은 규칙) · a0 = ref_front · a1 = a0 + canon_len
#   Alu 탐지 : 5S 정본 유닛의 alu 구간(1–318 · v100 SEG 정의)을 질의로 16-mer 씨앗 일치(양 가닥)
#              16-mer 무작위 기대치는 50 kb 당 0.003 개라 사실상 위양성 없음 · 씨앗이 있는 자리를 Alu 유사로 센다
#   출력 02_tables/JUNCTION_SEQ_ANNOT.tsv · 로그 04_logs/V143_SEQ.txt
import io, os, tarfile
import numpy as np, pandas as pd
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"; T = OUT + "/02_tables"
TAR = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11/12_flank/final/flank_final2.tar"
OFF = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11/13_junction50k/mapping_table/FLANK_OFFSET.tsv"
CANON = __import__("os").environ["CANON_FA"]
LOG = io.open(OUT + "/04_logs/V143_SEQ.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
CO = str.maketrans("ACGTacgtN", "TGCAtgcaN")
rc = lambda s: s.translate(CO)[::-1]
OUTMAX, INMAX, KS = 50000, 45000, 16

def binof(d):
    a = abs(d)
    if a < 2000: w = 50
    elif a < 10000: w = 200
    elif a <= 50000: w = 1000
    else: return None
    lo = (a // w) * w
    return lo if d > 0 else -lo

canon = "".join(x.strip() for x in io.open(CANON, encoding="utf8") if x[0] != ">").upper()
alu = canon[0:318]                                   # v100 SEG : ("alu", 1, 318)
KM = set()
for i in range(len(alu) - KS + 1):
    KM.add(alu[i:i + KS]); KM.add(rc(alu[i:i + KS]))
P("[0] 정본 유닛 %d bp · alu 구간 %d bp · 16-mer 씨앗 %d 개" % (len(canon), len(alu), len(KM)))

O = pd.read_csv(OFF, sep="\t", dtype={"slot": str})
key = {(r.sample, str(r.slot)): r for r in O.itertuples(index=False)}
def hd(txt):
    h = txt.split("\n", 1)[0]
    m = {}
    for kv in h[1:].split("|"):
        if "=" in kv: k, v = kv.split("=", 1); m[k] = v
    return m, "".join(txt.split("\n")[1:])

rows = []; n_ok = 0
with tarfile.open(TAR) as tf:
    for mi, m in enumerate(tf.getmembers()):
        if not m.name.endswith(".fa"): continue
        txt = tf.extractfile(m).read().decode()
        hdr, seq = hd(txt)
        base = os.path.basename(m.name)[:-3]
        tail = base.split("__")[-1]
        smp, slot = tail.rsplit("_", 1)
        smp = smp.replace("hiarr_", "").replace("arr_", "")
        r = key.get((smp, slot))
        if r is None: continue
        ref = (rc(seq) if r.rcf == "Y" else seq).upper()
        a0, a1 = int(r.ref_front), int(r.ref_front) + int(r.canon_len)
        if a1 > len(ref): continue
        n_ok += 1
        # Alu 씨앗 자리 표시
        hit = np.zeros(len(ref), dtype=bool)
        for i in range(len(ref) - KS + 1):
            if ref[i:i + KS] in KM: hit[i:i + KS] = True
        acc = {}
        for side, lo_out, hi_out, inside in (("원위", 0, a0, (a0, min(a0 + INMAX, a1))),
                                             ("근위", a1, len(ref), (max(a1 - INMAX, a0), a1))):
            for p in range(max(0, lo_out - 0), hi_out):
                pass
            # 바깥(flank)
            if side == "원위":
                rng = range(max(0, a0 - OUTMAX), a0)
                dfun = lambda p: a0 - p
            else:
                rng = range(a1, min(len(ref), a1 + OUTMAX))
                dfun = lambda p: p - a1 + 1
            for p in rng:
                b = binof(dfun(p))
                if b is None: continue
                a = acc.setdefault((side, b), [0, 0, 0])
                a[0] += 1
                if ref[p] == "C" and p + 1 < len(ref) and ref[p + 1] == "G": a[1] += 1
                if hit[p]: a[2] += 1
            # 안쪽(array)
            i0, i1 = inside
            for p in range(i0, i1):
                d = -((p - a0 + 1) if side == "원위" else (a1 - p))
                if abs(d) > INMAX: continue
                b = binof(d)
                if b is None: continue
                a = acc.setdefault((side, b), [0, 0, 0])
                a[0] += 1
                if ref[p] == "C" and p + 1 < len(ref) and ref[p + 1] == "G": a[1] += 1
                if hit[p]: a[2] += 1
        for (side, b), (bp, cg, al) in acc.items():
            rows.append(dict(sample=smp, slot=slot, side=side, b=b, bp=bp, n_cpg=cg, n_alu=al))
        if n_ok % 50 == 0: P("   %d 하플로 ... %s|%s" % (n_ok, smp, slot))
F = pd.DataFrame(rows)
F.to_csv(T + "/JUNCTION_SEQ_ANNOT.tsv", sep="\t", index=False)
P("[1] 저장 02_tables/JUNCTION_SEQ_ANNOT.tsv · 하플로 %d · 행 %d" % (n_ok, len(F)))
g = F.groupby(["side", "b"]).agg(bp=("bp", "sum"), cg=("n_cpg", "sum"), al=("n_alu", "sum"))
g["cpg_per_kb"] = 1000 * g.cg / g.bp; g["alu_frac"] = g.al / g.bp
for side in ("근위", "원위"):
    s = g.loc[side].sort_index()
    sel = [x for x in (-4000, -2000, -1000, 0, 500, 1000, 2000, 4000, 6000, 10000, 20000, 30000, 40000) if x in s.index]
    P("   %s CpG/kb : %s" % (side, " ".join("%d:%.1f" % (k, s.loc[k, "cpg_per_kb"]) for k in sel)))
    P("   %s Alu 비율: %s" % (side, " ".join("%d:%.2f" % (k, s.loc[k, "alu_frac"]) for k in sel)))
LOG.close()
