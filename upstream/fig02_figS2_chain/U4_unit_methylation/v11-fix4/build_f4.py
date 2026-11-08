# -*- coding: utf-8 -*-
# Build v11-fix4 (control order , 2026-09-22,)
#   Only the HG008-T, GIAB_SOMATIC_ONT rows are switched to CpG-only (CPGONLY=1). No other row changes by a single character.
#   1) raw : 12 files in the v11-fix2 methylation raw folder (R2 below) → in the 4 READUNIT and SITE files of 58-58 ao/jx only the HG008-T GIAB block is replaced, the other 8 files are copied unchanged
#   2) merged table : run merge_meth.py(4bc9cf68) on a folder holding only the HG008-T GIAB rows to get just that block, then splice it into the V11F2 merged table
#      (the merge_meth key is (sample, source, slot, array_pos), so rows of other samples are unaffected. Pilot : the CPGONLY=0 version must equal the V11F2 block)
#   3) combined tables : replace the dataset=v11 block of the extra1, extra2 and extra1+extra2 combined tables with V11F4
import io, os, sys, gzip, shutil, hashlib, subprocess, collections

NL = "\n"; TAB = "\t"
F = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions"
N = F + "/v11-fix4"
R2 = F + "/v11-fix2/08_methylation/raw"
C = N + "/01_code/local"
Q = N + "/09_tests"
KEY = ("HG008-T", "GIAB_SOMATIC_ONT")
LOG = []


def log(s):
    print(s); LOG.append(s)


def md5b(b):
    return hashlib.md5(b).hexdigest()


def md5L(L):
    return hashlib.md5("".join(x + NL for x in L).encode("utf8")).hexdigest()


def gzlines(p):
    with gzip.open(p, "rt", encoding="utf8", newline="") as f:
        return f.read().split(NL)[:-1]


def lines(p):
    return io.open(p, encoding="utf8", newline="").read().split(NL)[:-1]


def must_new(p):
    assert not os.path.exists(p), "이미 있음 " + p
    return p


def is_key(L):
    r = L.split(TAB, 2)
    return (r[0], r[1]) == KEY


# ---------------- 1. raw ----------------
RAWN = N + "/08_methylation/raw"
os.makedirs(RAWN, exist_ok=True)
for fn in sorted(os.listdir(R2)):
    src = R2 + "/" + fn
    dst = RAWN + "/" + fn
    if not fn.endswith(".tsv.gz"):
        continue
    EXIST = os.path.exists(dst) # raw already written, not rewritten; only its content is re-checked
    if "58-58" in fn:
        kind = "ao" if fn.endswith("_ao.tsv.gz") else "jx"
        typ = "METH_READUNIT" if fn.startswith("METH_READUNIT") else "METH_SITE"
        newf = Q + "/CPGONLY1/out_%s/%s_REPRO_ip-172-31-58-58.tsv.gz" % (kind, typ)
        oldL = gzlines(src)
        nw = gzlines(newf)
        assert nw[0] == oldL[0], "머리줄 다름 " + fn
        nw = nw[1:]
        assert all(is_key(L) for L in nw)
        idx = [i for i, L in enumerate(oldL) if i > 0 and is_key(L)]
        assert idx == list(range(idx[0], idx[-1] + 1))
        outL = oldL[:idx[0]] + nw + oldL[idx[-1] + 1:]
        if not EXIST:
            with gzip.open(dst, "wt", encoding="utf8", newline=NL) as g:
                g.write("".join(x + NL for x in outL))
        assert gzlines(dst) == outL
        back = gzlines(dst)
        o_rest = [L for L in oldL if not is_key(L)]
        n_rest = [L for L in back if not is_key(L)]
        log("[raw] %s : 교체 HG008-T GIAB %d → %d 행 · 전체 %d → %d · 다른 행 %d / %d md5 %s / %s 같음=%s"
            % (fn, len(idx), len(nw), len(oldL) - 1, len(back) - 1, len(o_rest), len(n_rest),
               md5L(o_rest), md5L(n_rest), o_rest == n_rest))
        assert o_rest == n_rest
        assert [L for L in back if is_key(L)] == nw
    else:
        if not EXIST:
            shutil.copyfile(src, dst)
        a, b = md5b(open(src, "rb").read()), md5b(open(dst, "rb").read())
        log("[raw] %s : copied as is, md5 %s same=%s" % (fn, b, a == b))
        assert a == b

# ---------------- 2. merged table ----------------
MERGE = C + "/merge_meth.py"
assert md5b(open(MERGE, "rb").read()) == "0ad651bd8314335432c779105c5133a7"


def merge_key_only(tag, ao_ru, jx_ru):
    """Build an input folder holding only the HG008-T GIAB READUNIT rows and run merge_meth. v11 rule : file names _a_jx and _b_ao (the file read first wins DUAL, i.e. junction first). The first run read _ao first and the pilot did not match"""
    d_in = Q + "/merged_%s/in" % tag
    d_out = Q + "/merged_%s/out" % tag
    os.makedirs(d_in, exist_ok=False)
    for nm, rows, hd in (("METH_READUNIT_ip-172-31-58-58_b_ao.tsv.gz",) + ao_ru, ("METH_READUNIT_ip-172-31-58-58_a_jx.tsv.gz",) + jx_ru):
        with gzip.open(d_in + "/" + nm, "wt", encoding="utf8", newline=NL) as g:
            g.write("".join(x + NL for x in [hd] + rows))
    r = subprocess.run([sys.executable, MERGE, d_in, d_out], capture_output=True, text=True, encoding="utf8",
                       env=dict(os.environ, PYTHONIOENCODING="utf-8"))
    io.open(Q + "/merged_%s/merge.log" % tag, "w", encoding="utf8", newline=NL).write(r.stdout + r.stderr)
    assert r.returncode == 0, r.stderr
    return lines(d_out + "/READMAP_METH_G8V4.tsv"), gzlines(d_out + "/METH_SITE_G8V4.tsv.gz")


def ru_rows(p):
    L = gzlines(p)
    return (L[1:], L[0])


raw2 = {k: gzlines(R2 + "/METH_READUNIT_ip-172-31-58-58_%s.tsv.gz" % k) for k in ("ao", "jx")}
old_ru = {k: ([L for L in raw2[k][1:] if is_key(L)], raw2[k][0]) for k in raw2}
new_ru = {k: ru_rows(Q + "/CPGONLY1/out_%s/METH_READUNIT_REPRO_ip-172-31-58-58.tsv.gz" % k) for k in ("ao", "jx")}

V2 = lines(F + "/v11-fix2/08_methylation/out/READMAP_METH_V11F2.tsv")
S2 = gzlines(F + "/v11-fix2/08_methylation/out/METH_SITE_V11F2.tsv.gz")
v2blk = [L for L in V2[1:] if is_key(L)]
s2blk = [L for L in S2[1:] if is_key(L)]

# Pilot : merge the block from the old rows (same as CPGONLY=0) → must equal the V11F2 block
pm, ps = merge_key_only("파일럿_옛행", old_ru["ao"], old_ru["jx"])
log("[파일럿] 옛 HG008-T 행 덩어리 합본 : READMAP %d 행 vs V11F2 덩어리 %d 행 같음=%s · SITE %d vs %d 같음=%s"
    % (len(pm) - 1, len(v2blk), pm[1:] == v2blk, len(ps) - 1, len(s2blk), ps[1:] == s2blk))
assert pm[0] == V2[0] and pm[1:] == v2blk and ps[1:] == s2blk

nm_, ns_ = merge_key_only("CPGONLY1", new_ru["ao"], new_ru["jx"])
log("[합본] 새 HG008-T 덩어리 : READMAP %d 행 · SITE %d 행" % (len(nm_) - 1, len(ns_) - 1))


def splice(old, newblk, name):
    idx = [i for i, L in enumerate(old) if i > 0 and is_key(L)]
    assert idx == list(range(idx[0], idx[-1] + 1))
    out = old[:idx[0]] + newblk + old[idx[-1] + 1:]
    o_rest = [L for L in old[1:] if not is_key(L)]
    n_rest = [L for L in out[1:] if not is_key(L)]
    log("[%s] 덩어리 %d → %d 행 · 전체 %d → %d · 다른 행 같음=%s (md5 %s)"
        % (name, len(idx), len(newblk), len(old) - 1, len(out) - 1, o_rest == n_rest, md5L(n_rest)))
    assert o_rest == n_rest
    # Check that sort order is kept (same as the sorted key of merge_meth)
    def k(L):
        r = L.split(TAB)
        return (r[0], r[1], r[2], int(r[3])) + ((int(r[4]),) if len(r) == 8 else ())
    ks = [k(L) for L in out[1:]]
    assert ks == sorted(ks), name + " 정렬 어긋남"
    return out


OUTD = N + "/08_methylation/out"
os.makedirs(OUTD, exist_ok=True)
V4 = splice(V2, nm_[1:], "READMAP_METH_V11F4")
S4 = splice(S2, ns_[1:], "METH_SITE_V11F4")
p = must_new(OUTD + "/READMAP_METH_V11F4.tsv")
io.open(p, "w", encoding="utf8", newline=NL).write("".join(x + NL for x in V4))
p = must_new(OUTD + "/METH_SITE_V11F4.tsv.gz")
with gzip.open(p, "wt", encoding="utf8", newline=NL) as g:
    g.write("".join(x + NL for x in S4))

# ---------------- 3. combined tables (with a dataset column) ----------------
for src, dstn in ((F + "/v11-add1_tissue_cancer/08_methylation/READMAP_METH_V11F2+add1_cancer4groups_noY7.tsv", "READMAP_METH_V11F4+추가1_암4군_Y7없음.tsv"),
                  (F + "/v11-add2_trio/08_methylation/READMAP_METH_V11F2+extra2_trios.tsv", "READMAP_METH_V11F4+extra2_trios.tsv"),
                  (F + "/v11-add2_trio/08_methylation/READMAP_METH_V11F2+add1+add2.tsv", "READMAP_METH_V11F4+추가1+추가2.tsv")):
    A = lines(src)
    assert A[0] == V2[0] + TAB + "dataset"
    vpart = [L.rsplit(TAB, 1)[0] for L in A[1:] if L.rsplit(TAB, 1)[1] == "v11"]
    assert vpart == V2[1:], "합본 속 v11 덩어리가 V11F2 와 다름 " + src
    vi = [i for i, L in enumerate(A) if i > 0 and L.rsplit(TAB, 1)[1] == "v11"]
    assert vi == list(range(1, len(V2)))
    out = [A[0]] + [L + TAB + "v11" for L in V4[1:]] + A[len(V2):]
    o_rest = A[len(V2):]
    log("[합본표] %s : v11 덩어리 %d → %d 행 · 다른 dataset %d 행 그대로 · 전체 %d → %d"
        % (dstn, len(V2) - 1, len(V4) - 1, len(o_rest), len(A) - 1, len(out) - 1))
    p = must_new(N + "/08_methylation/" + dstn)
    io.open(p, "w", encoding="utf8", newline=NL).write("".join(x + NL for x in out))

io.open(must_new(Q + "/build_f4.log"), "w", encoding="utf8", newline=NL).write(NL.join(LOG) + NL)
print("완료")
