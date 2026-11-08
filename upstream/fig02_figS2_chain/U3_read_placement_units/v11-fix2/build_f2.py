# -*- coding: utf-8 -*-
# v11-fix2: in the fix1 ARRAY_ONLY unit table (folder 06, final units) and ARRAY_ONLY methylation raw data (folder 08, methylation, raw *_ao),
#            replace only the rows of 9 samples with the new fix2 rows (inserted at the old row block, with an assert that each sample's old rows are contiguous).
#            junction files and rows of other samples are unchanged.
#   Arguments: kind (units|meth), new file folder (per-node subfolder <host>/), output folder
import io, os, sys, gzip, shutil, collections
NL = chr(10); TAB = chr(9)
T1 = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11-fix1"
KIND, NEWD, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)
HOSTS = {"ip-172-31-32-12": ["HG01081", "NA20282", "HG00133"],
         "ip-172-31-58-58": ["HG02559", "HG03270"],
         "ip-172-31-57-143": ["HG01106", "HG02624", "NA20850", "HG00350"]}
def rlines(p):
    op = gzip.open(p, "rt", encoding="utf8") if p.endswith(".gz") else io.open(p, encoding="utf8")
    x = op.read().split(NL)
    if x and x[-1] == "": x = x[:-1]
    return x
def wlines(p, lines):
    op = gzip.open(p, "wt", encoding="utf8", newline=NL) if p.endswith(".gz") else io.open(p, "w", encoding="utf8", newline=NL)
    op.write(NL.join(lines) + NL); op.close()
if KIND == "units":
    SRC = T1 + "/06_units_final"; files = sorted(os.listdir(SRC))
    def target(fn):
        for h in HOSTS:
            if fn == "AO_G8_UNITS_%s.mid.tsv" % h: return h, os.path.join(NEWD, h, "AO_G8_UNITS_%s.pad.tsv" % h), False
        return None
else:
    SRC = T1 + "/08_methylation/raw"; files = sorted(os.listdir(SRC))
    def target(fn):
        for h in HOSTS:
            if fn == "METH_READUNIT_%s_ao.tsv.gz" % h: return h, os.path.join(NEWD, h, "METH_READUNIT_%s.tsv.gz" % h), True
            if fn == "METH_SITE_%s_ao.tsv.gz" % h: return h, os.path.join(NEWD, h, "METH_SITE_%s.tsv.gz" % h), True
        return None
log = []
for fn in files:
    src = os.path.join(SRC, fn); dst = os.path.join(OUT, fn)
    t = target(fn)
    if t is None:
        shutil.copyfile(src, dst); log.append("복사 " + fn); continue
    h, newp, has_head = t
    old = rlines(src); new = rlines(newp)
    if has_head:
        assert old[0] == new[0], "머리줄 다름 " + fn
        head, old, new = [old[0]], old[1:], new[1:]
    else:
        head = []
    smps = HOSTS[h]
    byn = collections.defaultdict(list)
    for L in new:
        s = L.split(TAB, 1)[0]
        assert s in smps, "새 파일에 다른 개체 %s %s" % (fn, s)
        byn[s].append(L)
    out = list(old)
    for s in smps:
        idx = [i for i, L in enumerate(out) if L.split(TAB, 1)[0] == s]
        if idx:
            assert idx[-1] - idx[0] + 1 == len(idx), "옛 행 연속 아님 %s %s" % (fn, s)
            out = out[:idx[0]] + byn[s] + out[idx[-1] + 1:]
        else:
            out = out + byn[s]
        log.append("교체 %s : %s 옛 %d -> 새 %d 행" % (fn, s, len(idx), len(byn[s])))
    wlines(dst, head + out)
    log.append("   %s 전체 %d -> %d 행" % (fn, len(old), len(out)))
io.open(os.path.join(OUT, "_build_log.txt"), "w", encoding="utf8").write(NL.join(log) + NL)
print(NL.join(log))
