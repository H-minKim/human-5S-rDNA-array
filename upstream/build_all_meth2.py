# -*- coding: utf-8 -*-
# 완전 합본 v2 — Year7(조직) + 세포주 3 시료까지 넣는다
#   바탕 = v11-수정4 `READMAP_METH_V11F4+추가1+추가2.tsv` (59,486 행 · 글자 그대로 보존)
#   더함 = RU 에서 유도한 자리표 (유도 코드는 공식 merge_meth.py 와 같음이 검산됨
#          — HG008 700 · CASTLE 1,433 · COLOr 266 자리를 글자까지 재현)
#   Y7 은 공식 methmerge 산출 대신 RU 유도판을 쓴다 (이유는 문서 참고 · 둘 다 남긴다)
import io, os, sys, hashlib, collections
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
HERE = os.path.dirname(os.path.abspath(__file__))
V = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/"
BASE = V + "v11-fix4/08_methylation/READMAP_METH_V11F4+add1+add2.tsv"
OUTDIR = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar/02_tables/"
OUT = OUTDIR + "READMAP_METH_ALL_v2_2026-09-27.tsv"
LOG = io.open(os.path.join(HERE, "BUILD_ALL_METH2.txt"), "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()

base = io.open(BASE, encoding="utf8").read().rstrip("\n").split("\n")
hd = base[0].split("\t"); body = base[1:]
P("[바탕] %s · md5 %s · 행 %s" % (os.path.basename(BASE), md5(BASE), format(len(body), ",")))
keys5 = set()
for ln in body:
    p = ln.split("\t"); keys5.add((p[0], p[1], p[2], p[3], p[-1]))
P("   5 열 키(+dataset) 중복 : %d" % (len(body) - len(keys5)))

ADD = [("RM_Y7_LCL.tsv", "HPRC_Y7_LCL"), ("RM_Y7_PBMC.tsv", "HPRC_Y7_PBMC"),
       ("RM_CELL_K562.tsv", "CELLLINE_K562"),
       ("RM_CELL_GM12878_ONT.tsv", "CELLLINE_GM12878"),
       ("RM_CELL_GM12878_FS.tsv", "CELLLINE_GM12878FS")]
new = []
for fn, ds in ADD:
    L = io.open(os.path.join(HERE, fn), encoding="utf8").read().rstrip("\n").split("\n")
    assert L[0].split("\t") == hd[:-1], "열이 다르다"
    hit = 0
    for ln in L[1:]:
        p = ln.split("\t")
        k = (p[0], p[1], p[2], p[3], ds)
        if k in keys5: hit += 1; continue
        keys5.add(k); new.append(ln + "\t" + ds)
    P("[더함] %-26s dataset=%-20s 행 %s · 키 겹쳐 뺌 %d"
      % (fn, ds, format(len(L) - 1 - hit, ","), hit))

allrows = body + new
io.open(OUT, "w", encoding="utf8", newline="\n").write(
    "\t".join(hd) + "\n" + "\n".join(allrows) + "\n")
P()
P("[저장] %s" % OUT)
P("   행 %s (바탕 %s + 새로 %s) · md5 %s"
  % (format(len(allrows), ","), format(len(body), ","), format(len(new), ","), md5(OUT)))
chk = io.open(OUT, encoding="utf8").read().rstrip("\n").split("\n")
P("   바탕 행이 글자 그대로 남아 있나 : %s" % ("예" if chk[1:1 + len(body)] == body else "아니오"))
ds = collections.Counter(); smp = set(); src = collections.Counter()
for ln in chk[1:]:
    p = ln.split("\t"); ds[p[-1]] += 1; smp.add(p[0]); src[p[1]] += 1
P("   dataset : %s" % dict(ds))
P("   data_source : %s" % dict(src))
P("   시료 %d" % len(smp))
k4 = collections.Counter(); k5 = collections.Counter()
for ln in chk[1:]:
    p = ln.split("\t")
    k4[(p[0], p[1], p[2], p[3])] += 1; k5[(p[0], p[1], p[2], p[3], p[-1])] += 1
P("   키 중복 4 열 %d 자리 · 5 열(+dataset) %d 자리"
  % (sum(1 for v in k4.values() if v > 1), sum(1 for v in k5.values() if v > 1)))

# Y7 공식판도 같이 남긴다 (견주어 볼 수 있게)
OFF = OUTDIR + "READMAP_METH_Y7_공식methmerge_2026-09-27.tsv"
rows = ["\t".join(hd)]
for grp, ds_ in (("LCL", "HPRC_Y7_LCL_officialmerge"), ("PBMC", "HPRC_Y7_PBMC_officialmerge")):
    f = V + "v11-add1_tissue_cancer/08_methylation/node3/JX_%s/methmerge/out/READMAP_METH_G8V4.tsv" % grp
    L = io.open(f, encoding="utf8").read().rstrip("\n").split("\n")
    for ln in L[1:]: rows.append(ln + "\t" + ds_)
io.open(OFF, "w", encoding="utf8", newline="\n").write("\n".join(rows) + "\n")
P()
P("[같이 저장] %s (공식 methmerge 산출 · 견줌용) · md5 %s" % (os.path.basename(OFF), md5(OFF)))
P()
P("※ 정본 판본 폴더에는 쓰지 않았다 — 편입은 사항.")
LOG.close()
