# -*- coding: utf-8 -*-
"""FiberHMM BED12 -> 5S 유닛 좌표 · 전량판 (2026-10-11)

ex_<sample>__<pacbio|nanopore>/<tag>.calls_{tf,nucleosome,msp}.bed 를 읽어
MAP_V11 유닛 경계로 자르고 TSS 기준 rel 로 옮긴다.
재정렬 없음 — 좌표는 배치기(ao_g8cut) 가 준 것 그대로.
TSS = 유닛 서열에서 GTCTACGGC 의 마지막 출현 (a_lib 과 같은 규칙)
"""
import io
import os
import glob
import collections

NL = chr(10)
TB = chr(9)
H = os.path.expanduser
W = H("~/fhmm")
MAPF = H("~/MAP_V11.tsv")
ANCH = "GTCTACGGC"
HOST = os.uname()[1]
REFD = [H("~/pbmeth_v12d/bam/ref"), H("~/fsmeth_v12c/bam/ref")]


def rd_fa(p):
    d, nm, buf = {}, None, []
    for l in io.open(p, encoding="utf8"):
        if l.startswith(">"):
            if nm:
                d[nm] = "".join(buf).upper()
            nm, buf = l[1:].strip().split()[0], []
        else:
            buf.append(l.strip())
    if nm:
        d[nm] = "".join(buf).upper()
    return d


U = collections.defaultdict(list)
f = io.open(MAPF, encoding="utf8")
ix = dict((n, i) for i, n in enumerate(f.readline().rstrip(NL).split(TB)))
for L in f:
    r = L.rstrip(NL).split(TB)
    if r[ix["seg_type"]] == "UNIT":
        U[(r[0], r[1])].append((int(r[ix["unit_start"]]), int(r[ix["unit_end"]]),
                                r[ix["array_pos"]]))

REFP = {}
for d in REFD:
    for p in glob.glob(os.path.join(d, "*.fa")):
        REFP.setdefault(os.path.basename(p).split(".")[0], []).append(p)

out = io.open(os.path.join(W, "FHMM_UNITS_%s.tsv" % HOST), "w", encoding="utf8", newline=NL)
out.write(TB.join(["sample", "plat", "slot", "array_pos", "read_id", "cls",
                   "rel_start", "rel_end", "size"]) + NL)
cov = io.open(os.path.join(W, "FHMM_READCOV_%s.tsv" % HOST), "w", encoding="utf8", newline=NL)
cov.write(TB.join(["sample", "plat", "slot", "array_pos", "read_id",
                   "rel_lo", "rel_hi"]) + NL)

TOT = collections.Counter()
SEEN = set()
for exd in sorted(glob.glob(os.path.join(W, "out", "ex_*"))):
    tag = os.path.basename(exd)[3:]
    if "__" not in tag:
        continue
    smp, plat = tag.split("__", 1)
    paths = REFP.get(smp)
    if not paths:
        continue
    want = [p for p in paths if (("pbmeth" in p) == (plat == "pacbio"))]
    if not want:
        continue
    SEQ = rd_fa(want[0])
    UT = {}
    for rn, s in SEQ.items():
        sm_, sl = rn.rsplit("_h", 1)
        for (us, ue, ap) in U.get((sm_, sl), []):
            u = s[us - 1:ue]
            p = u.rfind(ANCH)
            if p < 0:
                p = len(u) - 119
            UT.setdefault(rn, []).append((us, ue, ap, us + p))
    for rn in UT:
        UT[rn].sort()

    for cls, suf in (("tf", "_tf.bed"), ("nuc", "_nucleosome.bed"), ("msp", "_msp.bed")):
        bp = os.path.join(exd, tag + ".calls" + suf)
        if not os.path.exists(bp):
            continue
        for L in io.open(bp, encoding="utf8"):
            x = L.rstrip(NL).split(TB)
            if len(x) < 12:
                continue
            rn, cs, rid = x[0], int(x[1]), x[3]
            if rn not in UT:
                continue
            sm_, sl = rn.rsplit("_h", 1)
            bl = [int(v) for v in x[10].rstrip(",").split(",") if v != ""]
            bs = [int(v) for v in x[11].rstrip(",").split(",") if v != ""]
            # read 가 덮는 범위는 어느 class 에서든 한 번만 적는다
            #   (tf BED 에만 적으면 TF 가 없는 read 가 분모에서 빠진다 - 2026-10-11 수정)
            r0, r1 = cs + 1, int(x[2])
            for (us, ue, ap, tss) in UT[rn]:
                if min(r1, ue) - max(r0, us) < 200:        # 이 유닛을 충분히 덮는 read 만
                    continue
                if (sm_, plat, sl, ap, rid) in SEEN:
                    continue
                SEEN.add((sm_, plat, sl, ap, rid))
                # 덮개를 유닛 안으로 자르지 않는다 — 탠덤이라 flank 는 이웃 유닛의 실제 서열이다.
                #   자르면 유닛 끝(rel ~ +118)에서 분모가 끊겨 가짜 봉우리가 선다 (2026-10-11 수정)
                a, b = max(r0, tss - 1200), min(r1, tss + 1200)
                cov.write(TB.join([sm_, plat, sl, ap, rid,
                                   "%d" % (a - tss), "%d" % (b - tss)]) + NL)
            for sz, st in zip(bl, bs):
                g0 = cs + st + 1
                g1 = g0 + sz - 1
                for (us, ue, ap, tss) in UT[rn]:
                    if g1 < tss - 1200 or g0 > tss + 1200:
                        continue
                    out.write(TB.join([sm_, plat, sl, ap, rid, cls,
                                       "%d" % (g0 - tss), "%d" % (g1 - tss),
                                       "%d" % sz]) + NL)
                    TOT[(plat, cls)] += 1
out.close()
cov.close()
for k in sorted(TOT):
    print("  %-9s %-4s %9d" % (k[0], k[1], TOT[k]))
print("[%s] -> FHMM_UNITS_%s.tsv" % (HOST, HOST))
