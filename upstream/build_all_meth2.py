# -*- coding: utf-8 -*-
# The full merge, version 2: the Year 7 tissues and three cell-line samples are added.
#   The basis is the v11 revision 4 table READMAP_METH_V11F4 with the two additions, 59,486 rows, preserved character for character.
#   Added: the per-site table derived from RU. The derivation was checked to match the official merge_meth.py
#          by reproducing 700 HG008, 1,433 CASTLE and 266 COLO sites character for character.
#   For Year 7 the RU-derived version is used in place of the official methmerge output; both are kept.
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
P("[base] %s, md5 %s, rows %s" % (os.path.basename(BASE), md5(BASE), format(len(body), ",")))
keys5 = set()
for ln in body:
    p = ln.split("\t"); keys5.add((p[0], p[1], p[2], p[3], p[-1]))
P("   duplicates of the 5-column key (+dataset): %d" % (len(body) - len(keys5)))

ADD = [("RM_Y7_LCL.tsv", "HPRC_Y7_LCL"), ("RM_Y7_PBMC.tsv", "HPRC_Y7_PBMC"),
       ("RM_CELL_K562.tsv", "CELLLINE_K562"),
       ("RM_CELL_GM12878_ONT.tsv", "CELLLINE_GM12878"),
       ("RM_CELL_GM12878_FS.tsv", "CELLLINE_GM12878FS")]
new = []
for fn, ds in ADD:
    L = io.open(os.path.join(HERE, fn), encoding="utf8").read().rstrip("\n").split("\n")
    assert L[0].split("\t") == hd[:-1], "columns differ"
    hit = 0
    for ln in L[1:]:
        p = ln.split("\t")
        k = (p[0], p[1], p[2], p[3], ds)
        if k in keys5: hit += 1; continue
        keys5.add(k); new.append(ln + "\t" + ds)
    P("[added] %-26s dataset=%-20s rows %s, dropped for key overlap %d"
      % (fn, ds, format(len(L) - 1 - hit, ","), hit))

allrows = body + new
io.open(OUT, "w", encoding="utf8", newline="\n").write(
    "\t".join(hd) + "\n" + "\n".join(allrows) + "\n")
P()
P("[saved] %s" % OUT)
P("   rows %s (base %s + new %s), md5 %s"
  % (format(len(allrows), ","), format(len(body), ","), format(len(new), ","), md5(OUT)))
chk = io.open(OUT, encoding="utf8").read().rstrip("\n").split("\n")
P("   base rows kept verbatim: %s" % ("yes" if chk[1:1 + len(body)] == body else "아니오"))
ds = collections.Counter(); smp = set(); src = collections.Counter()
for ln in chk[1:]:
    p = ln.split("\t"); ds[p[-1]] += 1; smp.add(p[0]); src[p[1]] += 1
P("   dataset : %s" % dict(ds))
P("   data_source : %s" % dict(src))
P("   samples %d" % len(smp))
k4 = collections.Counter(); k5 = collections.Counter()
for ln in chk[1:]:
    p = ln.split("\t")
    k4[(p[0], p[1], p[2], p[3])] += 1; k5[(p[0], p[1], p[2], p[3], p[-1])] += 1
P("   key duplicates: 4 columns %d positions, 5 columns (+dataset) %d positions"
  % (sum(1 for v in k4.values() if v > 1), sum(1 for v in k5.values() if v > 1)))

# the official Year 7 version is kept alongside, for comparison
OFF = OUTDIR + "READMAP_METH_Y7_공식methmerge_2026-09-27.tsv"
rows = ["\t".join(hd)]
for grp, ds_ in (("LCL", "HPRC_Y7_LCL_officialmerge"), ("PBMC", "HPRC_Y7_PBMC_officialmerge")):
    f = V + "v11-add1_tissue_cancer/08_methylation/node3/JX_%s/methmerge/out/READMAP_METH_G8V4.tsv" % grp
    L = io.open(f, encoding="utf8").read().rstrip("\n").split("\n")
    for ln in L[1:]: rows.append(ln + "\t" + ds_)
io.open(OFF, "w", encoding="utf8", newline="\n").write("\n".join(rows) + "\n")
P()
P("[also saved] %s (official methmerge output, for comparison), md5 %s" % (os.path.basename(OFF), md5(OFF)))
P()
P("note: not written to the canonical version folder.")
LOG.close()
