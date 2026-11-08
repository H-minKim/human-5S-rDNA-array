# -*- coding: utf-8 -*-
# 같은 대상, 서열과 좌표 전부
#   단계 : ① array 잘라내기
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 과제 1-A: 조립본이 있는 모든 하플로타입에서 5S array 서열을 뽑아 보관한다
#   (asm_units3.py 와 같은 방법 · 다만 서열과 좌표를 전부 저장한다)
#
#   이것이 조립본 좌표계 consensus map 의 원자료가 된다.
#   좌표를 남겨야 나중에 bigwig 메틸을 같은 자리에 붙일 수 있다
#   (asm_units3.py 는 a2-a1 만 남기고 좌표를 버려서 다시 구해야 했다).
#
# URL 은 반드시 색인 CSV 에서 읽는다. 지어내지 않는다.
# 유닛 계수는 Python sorted(). awk asorti 는 문자열 정렬이라 좌표를 뒤섞는다.
import io
import os
import sys
import subprocess
import collections
import time

W = os.path.expanduser("~/asmarr")
os.makedirs(W, exist_ok=True)
GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
IDX = os.path.expanduser("~/telomere_pilot/data/assemblies_release2_v1.0.index.csv")
PP = os.path.expanduser("~/asmunits3/probe.fa")
TH = os.environ.get("TH", "6")
TAG = os.environ.get("TAG", "x")
NL = "\n"

gs = "".join(x.strip() for x in io.open(GENE, encoding="utf8") if x[0] != ">")
comp = str.maketrans("ACGT", "TGCA")
MOT = []
for off in (5, 30, 55, 80):
    m = gs[off:off + 25]
    if len(m) == 25:
        MOT.append(m)
        MOT.append(m.translate(comp)[::-1])

URL = collections.defaultdict(list)
for i, L in enumerate(io.open(IDX, encoding="utf8")):
    f = L.rstrip(NL).split(",")
    if i == 0:
        col = dict((c, k) for k, c in enumerate(f))
        continue
    if len(f) <= col["assembly"]:
        continue
    u = f[col["assembly"]].replace(
        "s3://human-pangenomics/",
        "https://s3-us-west-2.amazonaws.com/human-pangenomics/")
    URL[f[col["sample_id"]]].append((f[col["haplotype"]], f[col["assembly_name"]], u))


def units_of(seq, tag):
    """유닛 시작좌표 사슬을 돌려준다 (array 안 0-based)"""
    fa = os.path.join(W, "a_%s.fa" % tag)
    io.open(fa, "w", newline=NL).write(">arr" + NL + seq + NL)
    db = os.path.join(W, "d_%s" % tag)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   check=True, capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                        "-evalue", "1e-5", "-max_target_seqs", "100000000",
                        "-num_threads", TH], capture_output=True, text=True, check=True)
    st = set()
    pl = mi = 0
    for L in r.stdout.splitlines():
        f = L.split("\t")
        if len(f) >= 2:
            a, b = int(f[0]), int(f[1])
            st.add(min(a, b))
            if a < b:
                pl += 1
            else:
                mi += 1
    K = sorted(x - 1 for x in st)
    ok = sum(1 for i in range(len(K) - 1) if 2100 <= K[i + 1] - K[i] <= 2400)
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    return ok, K, (1 if pl >= mi else -1)


OUT = io.open(os.path.join(W, "asmarr_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
for smp in sys.argv[1:]:
    if smp not in URL:
        sys.stderr.write("  색인 없음 %s\n" % smp)
        continue
    for hap, anm, u in URL[smp]:
        tag = "%s_%s" % (smp, hap)
        if os.path.exists(os.path.join(W, "arr_%s.fa" % tag)):
            continue
        t0 = time.time()
        pat = "|".join(MOT)
        cmd = "curl -sfL " + repr(u) + " | zcat | grep -E " + repr("^>|" + pat)
        pr = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True,
                              stderr=subprocess.DEVNULL)
        cur = None
        cnt = collections.Counter()
        for L in pr.stdout:
            if L[0] == ">":
                cur = L[1:].split()[0]
            elif cur:
                cnt[cur] += 1
        pr.wait()
        if not cnt:
            OUT.write("\t".join([smp, hap, anm, "-", "-", "-", "-", "-",
                                 "NO_MOTIF"]) + NL)
            OUT.flush()
            continue
        fai = os.path.join(W, tag + ".fai")
        gzi = os.path.join(W, tag + ".gzi")
        try:
            subprocess.run(["curl", "-sfL", "-o", fai, u + ".fai"],
                           check=True, timeout=300)
            subprocess.run(["curl", "-sfL", "-o", gzi, u + ".gzi"],
                           check=True, timeout=300)
        except Exception:
            OUT.write("\t".join([smp, hap, anm, "-", "-", "-", "-", "-",
                                 "NO_INDEX"]) + NL)
            OUT.flush()
            continue
        SZ = {}
        for L in io.open(fai, encoding="utf8"):
            f = L.rstrip(NL).split("\t")
            if len(f) >= 2:
                SZ[f[0]] = int(f[1])
        found = None
        for c, n in cnt.most_common(4):
            if c not in SZ:
                continue
            rr = subprocess.run(["samtools", "faidx", "--fai-idx", fai,
                                 "--gzi-idx", gzi, u, "%s:1-%d" % (c, SZ[c])],
                                capture_output=True, text=True)
            if rr.returncode != 0 or len(rr.stdout) < 10000:
                continue
            reg = os.path.join(W, tag + ".reg.fa")
            io.open(reg, "w", newline=NL).write(rr.stdout)
            paf = subprocess.run(["minimap2", "-x", "asm20", "-t", TH,
                                  "--secondary=no", PP, reg],
                                 capture_output=True, text=True)
            d = {}
            for LL in paf.stdout.splitlines():
                f = LL.split("\t")
                if len(f) < 12 or int(f[9]) < 20000:
                    continue
                if f[5] not in d or int(f[9]) > d[f[5]][2]:
                    d[f[5]] = (int(f[2]), int(f[3]), int(f[9]))
            if "PROX50k" in d and "DIST50k" in d:
                a = sorted([d["PROX50k"][0], d["PROX50k"][1],
                            d["DIST50k"][0], d["DIST50k"][1]])
                seq = "".join(io.open(reg, encoding="utf8").read().split(NL)[1:])
                # 근위 탐침이 array 앞쪽인지 뒤쪽인지 기록한다 (참조 창 방향 함정)
                pmid = (d["PROX50k"][0] + d["PROX50k"][1]) // 2
                side = "PROX_FIRST" if pmid < a[1] else "DIST_FIRST"
                found = (c, a[1], a[2], seq[a[1]:a[2]], side)
                os.remove(reg)
                break
            os.remove(reg)
        for f_ in (fai, gzi):
            try:
                os.remove(f_)
            except OSError:
                pass
        if not found:
            top = ",".join("%s:%d" % (c, n) for c, n in cnt.most_common(2))
            OUT.write("\t".join([smp, hap, anm, top, "-", "-", "-", "-",
                                 "NO_BOTH_PROBE"]) + NL)
            OUT.flush()
            sys.stderr.write("  %s h%s 탐침 양쪽 없음 (%.0fs)\n"
                             % (smp, hap, time.time() - t0))
            sys.stderr.flush()
            continue
        c, a1, a2, seq, side = found
        ok, K, sg = units_of(seq, tag)
        io.open(os.path.join(W, "arr_%s.fa" % tag), "w", newline=NL).write(
            ">%s|%s|%d|%d|%s|%s%s%s%s" % (tag, c, a1, a2, side,
                                          "+" if sg > 0 else "-", NL, seq, NL))
        io.open(os.path.join(W, "units_%s.txt" % tag), "w", newline=NL).write(
            NL.join(str(x) for x in K) + NL)
        OUT.write("\t".join([smp, hap, anm, c, str(a1), str(a2),
                             str(a2 - a1), str(ok), "OK"]) + NL)
        OUT.flush()
        sys.stderr.write("  %s h%s  %s · %d-%d (%d bp) · 유닛 %d · %s %s (%.0fs)\n"
                         % (smp, hap, c.split("#")[-1][:18], a1, a2, a2 - a1, ok,
                            side, "+" if sg > 0 else "-", time.time() - t0))
        sys.stderr.flush()
OUT.close()
print("done")
