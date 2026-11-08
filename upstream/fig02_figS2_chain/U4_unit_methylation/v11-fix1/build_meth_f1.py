# -*- coding: utf-8 -*-
# 12 v11-fix1 methylation raw files: from v11 raw, only the rows of two individuals are new (inserted in place of the old rows, headers confirmed identical). The rest are copied
import io, os, gzip, shutil, collections
NL = chr(10); TAB = chr(9)
V = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11/08_methylation/raw"
G = __import__("os").environ["RDNA_WORK_DIR"] + "/cmp/g8place/fix1"
M = G + "/meth_new"; OUT = G + "/meth_raw_f1"
os.makedirs(OUT, exist_ok=True)
SMP = {"ip-172-31-32-12": "HG00133", "ip-172-31-58-58": "HG03270"}
log = []
for fn in sorted(os.listdir(V)):
    src = os.path.join(V, fn); dst = os.path.join(OUT, fn)
    base = fn.replace(".tsv.gz", "")          # METH_READUNIT_<host>_<k>
    kind, host, k = ("READUNIT" if "READUNIT" in fn else "SITE"), base.split("_")[2], base.split("_")[-1]
    if host not in SMP:
        shutil.copyfile(src, dst); log.append("복사 " + fn); continue
    smp = SMP[host]
    newp = os.path.join(M, "NEW_%s_%s_%s.tsv.gz" % (kind, host, k))
    old = gzip.open(src, "rt", encoding="utf8").read().split(NL)
    if old[-1] == "": old = old[:-1]
    new = gzip.open(newp, "rt", encoding="utf8").read().split(NL)
    if new[-1] == "": new = new[:-1]
    assert old[0] == new[0], "머리줄 다름 " + fn
    body_new = new[1:]
    assert all(L.split(TAB, 1)[0] == smp for L in body_new), "다른 개체 행 " + fn
    idx = [i for i, L in enumerate(old) if i > 0 and L.split(TAB, 1)[0] == smp]
    if idx:
        assert idx[-1] - idx[0] + 1 == len(idx), "옛 행 연속 아님 " + fn
        out = old[:idx[0]] + body_new + old[idx[-1] + 1:]
    else:
        out = old + body_new
    with gzip.open(dst, "wt", encoding="utf8", newline=NL) as g:
        g.write(NL.join(out) + NL)
    log.append("교체 %s : %s 옛 %d -> 새 %d 행 · 전체 %d -> %d" % (fn, smp, len(idx), len(body_new), len(old), len(out)))
io.open(OUT + "/_build_log.txt", "w", encoding="utf8").write(NL.join(log) + NL)
print(NL.join(log))
