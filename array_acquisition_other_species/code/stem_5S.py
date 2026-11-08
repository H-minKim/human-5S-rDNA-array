# -*- coding: utf-8 -*-
# GSE175951 Stage III(GV) 난모세포 STEM-seq 가공 CpG 표(danRer7) 로 5S 구역 메틸을 잰다.
#   구역은 danRer7 서열에서 직접 찾은 유전자 히트(hits_dr7.tsv)로 만든다.
#   chr4 : 모계형 유전자가 3 kb 안으로 이어진 무리(≥10 개)를 array 로 본다.
#   코딩부 = 유전자 · 스페이서 = array 안 나머지 · flank = 양쪽 100 kb(경계 5 kb 띄움, 다른 array 제외)
#   구간 : CpG 자리 부트스트랩(피복이 자리당 ~1 이라 자리 ≈ read 1 개)
import bisect
import collections
import gzip
import io
import random

random.seed(0)
TAB = "\t"
hits = collections.defaultdict(list)
for L in io.open("hits_dr7.tsv"):
    c, a, b, st, typ, dm, ds = L.rstrip("\n").split(TAB)
    if (c == "chr4" and typ == "maternal") or (c == "chr18" and typ == "somatic"):
        hits[c].append((int(a), int(b)))
regions = collections.defaultdict(list)   # (locus, part) -> [(c,s,e)]
for c, g in hits.items():
    g.sort()
    groups = [[g[0]]]
    for x in g[1:]:
        # chr18 유닛은 1,980 bp 이고 danRer7 에서 빠진 유닛이 있어 간격이 6 kb 까지 벌어진다
        if x[0] - groups[-1][-1][1] <= (3000 if c == "chr4" else 12000):
            groups[-1].append(x)
        else:
            groups.append([x])
    groups = [gr for gr in groups if len(gr) >= (10 if c == "chr4" else 5)]
    spans = [(gr[0][0], gr[-1][1]) for gr in groups]
    loc = "chr4_maternal" if c == "chr4" else "chr18_somatic"
    for gr in groups:
        s0, e0 = gr[0][0], gr[-1][1]
        p = s0
        for a, b in gr:
            regions[(loc, "코딩부")].append((c, a, b))
            if a > p:
                regions[(loc, "스페이서")].append((c, p, a - 1))
            p = b + 1
        for fs, fe, side in ((s0 - 105000, s0 - 5001, "flankL"), (e0 + 5001, e0 + 105000, "flankR")):
            # 다른 array 와 겹치는 flank 부분은 뺀다
            pieces = [(max(1, fs), fe)]
            for (xs, xe) in spans:
                new = []
                for (ps, pe) in pieces:
                    if xe < ps or xs > pe:
                        new.append((ps, pe))
                    else:
                        if xs > ps:
                            new.append((ps, xs - 1))
                        if xe < pe:
                            new.append((xe + 1, pe))
                pieces = new
            for ps, pe in pieces:
                regions[(loc, "바깥_" + side)].append((c, ps, pe))
    print(loc, "array", len(groups), "· 유전자", sum(len(gr) for gr in groups),
          "· span 합 %d bp" % sum(e - s + 1 for s, e in spans))

# 빠른 판정용 색인
idx = {}
for key, iv in regions.items():
    by = collections.defaultdict(list)
    for c, s, e in iv:
        by[c].append((s, e))
    for c in by:
        by[c].sort()
    idx[key] = {c: ([s for s, e in v], [e for s, e in v]) for c, v in by.items()}


def inside(key, c, p):
    t = idx[key].get(c)
    if not t:
        return False
    st, en = t
    i = bisect.bisect_right(st, p) - 1
    return i >= 0 and p <= en[i]


out = io.open("STEM_5S.tsv", "w", encoding="utf8", newline="\n")
out.write(TAB.join(["sample", "locus", "part", "n_site", "n_call", "n_meth", "frac", "lo95", "hi95"]) + "\n")
gw = {}
for samp, fn in (("oocyte_StageIII_rep1", "GSM5351778.mCG.txt.gz"), ("oocyte_StageIII_rep2", "GSM5351780.mCG.txt.gz")):
    sites = collections.defaultdict(list)
    tm = tu = 0
    for L in gzip.open(fn, "rt"):
        f = L.split()
        c, p, m, u = f[0], int(f[1]), int(f[2]), int(f[3])
        tm += m
        tu += u
        if c not in ("chr4", "chr18"):
            continue
        for key in regions:
            if inside(key, c, p):
                sites[key].append((m, u))
                break
    gw[samp] = tm / float(tm + tu)
    print("%s 전 게놈 CpG 메틸 %.3f" % (samp, gw[samp]))
    for key in sorted(regions):
        v = sites.get(key, [])
        m = sum(x[0] for x in v)
        n = sum(x[0] + x[1] for x in v)
        if n == 0:
            out.write(TAB.join([samp, key[0], key[1], "0", "0", "0", "-", "-", "-"]) + "\n")
            continue
        bs = []
        for _ in range(2000):
            s = [v[random.randrange(len(v))] for _ in range(len(v))]
            nn = sum(x[0] + x[1] for x in s)
            bs.append(sum(x[0] for x in s) / float(nn))
        bs.sort()
        out.write(TAB.join([samp, key[0], key[1], str(len(v)), str(n), str(m), "%.4f" % (m / float(n)),
                            "%.4f" % bs[50], "%.4f" % bs[1949]]) + "\n")
out.close()
print("완료")

# ── 코딩부/스페이서 − flank(L+R) 차이 : 자리 부트스트랩, 독립 표본(정렬 전 차이) ──
import numpy as np
rng = np.random.default_rng(1)


def bsamp(v, nb=4000):
    a = np.array(v, dtype=float)
    if len(a) == 0:
        return None
    ix = rng.integers(0, len(a), size=(nb, len(a)))
    return a[ix, 0].sum(1) / (a[ix, 0] + a[ix, 1]).sum(1)


dout = io.open("STEM_5S_DIFF.tsv", "w", encoding="utf8", newline="\n")
dout.write(TAB.join(["sample", "locus", "contrast", "diff", "lo95", "hi95"]) + "\n")
for samp, fn in (("oocyte_StageIII_rep1", "GSM5351778.mCG.txt.gz"), ("oocyte_StageIII_rep2", "GSM5351780.mCG.txt.gz")):
    sites = collections.defaultdict(list)
    for L in gzip.open(fn, "rt"):
        f = L.split()
        c = f[0]
        if c not in ("chr4", "chr18"):
            continue
        p = int(f[1])
        for key in regions:
            if inside(key, c, p):
                sites[key].append((int(f[2]), int(f[3])))
                break
    for loc in ("chr4_maternal", "chr18_somatic"):
        fl = sites.get((loc, "바깥_flankL"), []) + sites.get((loc, "바깥_flankR"), [])
        bf = bsamp(fl)
        for part in ("코딩부", "스페이서"):
            v = sites.get((loc, part), [])
            bv = bsamp(v)
            if bv is None or bf is None:
                continue
            pt = sum(x[0] for x in v) / float(sum(x[0] + x[1] for x in v)) - sum(x[0] for x in fl) / float(sum(x[0] + x[1] for x in fl))
            lo, hi = np.percentile(bv - bf, [2.5, 97.5])
            dout.write(TAB.join([samp, loc, part + "-flank", "%+.4f" % pt, "%+.4f" % lo, "%+.4f" % hi]) + "\n")
dout.close()
print("차이 완료")
