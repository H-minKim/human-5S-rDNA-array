# -*- coding: utf-8 -*-
# v11-fix1 read unit tables (9 files): from the v11 final unit folder (path V below), replace only the rows of two individuals with new rows (inserted at the old rows' position). Other files are copied unchanged
import io, os, sys, shutil, hashlib, collections
NL = chr(10)
V = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11/06_units_final"
G = __import__("os").environ["RDNA_WORK_DIR"] + "/cmp/g8place/fix1"
OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
REPL = {
    "AO_G8_UNITS_ip-172-31-32-12.mid.tsv": ("HG00133", G + "/units_b/AO_G8_UNITS_ip-172-31-32-12.pad.tsv"),
    "AO_G8_UNITS_ip-172-31-58-58.mid.tsv": ("HG03270", G + "/units_b/AO_G8_UNITS_ip-172-31-58-58.pad.tsv"),
    "JX_G8_UNITS_main_ip-172-31-32-12.mid.tsv": ("HG00133", G + "/units/JX_G8_UNITS_main_ip-172-31-32-12.pad.tsv"),
    "JX_G8_UNITS_main_ip-172-31-58-58.mid.tsv": ("HG03270", G + "/units/JX_G8_UNITS_main_ip-172-31-58-58.pad.tsv"),
}
log = []
for fn in sorted(os.listdir(V)):
    src = os.path.join(V, fn); dst = os.path.join(OUT, fn)
    if fn not in REPL:
        shutil.copyfile(src, dst); log.append("복사 %s" % fn); continue
    smp, newp = REPL[fn]
    old = io.open(src, encoding="utf8").read().split(NL)
    if old and old[-1] == "": old = old[:-1]
    new = [L for L in io.open(newp, encoding="utf8").read().split(NL) if L]
    assert all(L.split("\t")[0] == smp for L in new), "새 표에 다른 개체"
    idx = [i for i, L in enumerate(old) if L.split("\t", 1)[0] == smp]
    assert idx and idx[-1] - idx[0] + 1 == len(idx), "옛 행이 연속이 아님 %s" % fn
    out = old[:idx[0]] + new + old[idx[-1] + 1:]
    io.open(dst, "w", encoding="utf8", newline=NL).write(NL.join(out) + NL)
    log.append("교체 %s : %s 옛 %d 행 -> 새 %d 행 · 전체 %d -> %d" % (fn, smp, len(idx), len(new), len(old), len(out)))
io.open(os.path.join(OUT, "_build_log.txt"), "w", encoding="utf8").write(NL.join(log) + NL)
print(NL.join(log))
