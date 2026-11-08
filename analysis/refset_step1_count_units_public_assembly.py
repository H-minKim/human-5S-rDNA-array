# -*- coding: utf-8 -*-
# HPRC 공개 조립본에서 유닛 수 (노드 ~/asm_units3.py 에서 가져옴)
#   단계 : ① array 잘라내기
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# HPRC release2 조립본에서 5S array 만 뽑아 flank-flank 길이와 유닛 수를 센다
#
# 왜 v3 인가
#  v1: 조립본 3 GB 전체를 minimap2 질의로 -> 하플로타입당 200 초. 느리다.
#  v2: '가장 큰 contig 끝 45 Mb' 로 좁힘 -> 6 초. 그러나 조립본마다 array 가 있는
#      contig 가 다르다 (HG01943 은 11 번째 contig JBHIIK010000077.1).
#      contig 를 추측하면 안 된다.
#  v3: contig 를 서열로 찾는다. 스트림을 grep 으로 훑어 5S 모티프가 있는 contig 를
#      찾고, 그 contig 만 remote faidx 로 받아 탐침 정렬한다.
#
# 유닛 계수는 Python sorted(). awk asorti 는 인덱스를 문자열 정렬해 좌표를 뒤섞는다
#   (2026-08-27 정정문 — 옛 "조립본 14 유닛 축약"의 원인).
# URL 은 반드시 인덱스 CSV 에서 읽는다. 지어내지 않는다.
import io, os, sys, subprocess, collections, time

W = os.path.expanduser("~/asmunits3")
os.makedirs(W, exist_ok=True)
REF = os.path.expanduser("~/ectopic_sp/chm13v2.0.fa")
GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
IDX = os.path.expanduser("~/telomere_pilot/data/assemblies_release2_v1.0.index.csv")
TH = os.environ.get("TH", "8")
TAG = os.environ.get("TAG", "x")
NL = "\n"
PJ, DJ = 227744884, 228025900

# ---- 탐침: junction 바깥 50 kb (array 는 넣지 않는다 — 탠덤이라 아무데나 붙는다) ----
pp = os.path.join(W, "probe.fa")
if not os.path.exists(pp):
    with io.open(pp, "w", newline=NL) as g:
        for nm, reg in (("PROX50k", "chr1:%d-%d" % (PJ - 50000, PJ - 1)),
                        ("DIST50k", "chr1:%d-%d" % (DJ + 1, DJ + 50000))):
            s = subprocess.run(["samtools", "faidx", REF, reg],
                               capture_output=True, text=True, check=True).stdout
            g.write(">" + nm + NL + "".join(s.split(NL)[1:]) + NL)

# ---- 5S 모티프: 유전자 서열에서 25-mer 몇 개 (양가닥) ----
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
        col = dict((c, k) for k, c in enumerate(f)); continue
    if len(f) <= col["assembly"]:
        continue
    u = f[col["assembly"]].replace("s3://human-pangenomics/",
                                   "https://s3-us-west-2.amazonaws.com/human-pangenomics/")
    URL[f[col["sample_id"]]].append((f[col["haplotype"]], f[col["assembly_name"]], u))


def find_contig(u, tag):
    """스트림을 grep 으로 훑어 5S 모티프가 가장 많은 contig 를 돌려준다."""
    pat = "|".join(MOT)
    cmd = ("curl -sfL '%s' | zcat | grep -E '^>|%s'" % (u, pat))
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
    return cnt


def units_of(seq, tag):
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
    for L in r.stdout.splitlines():
        f = L.split("\t")
        if len(f) >= 2:
            st.add(min(int(f[0]), int(f[1])))
    K = sorted(st)
    ok = sum(1 for i in range(len(K) - 1) if 2100 <= K[i + 1] - K[i] <= 2400)
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    return ok, len(K) - 1, len(K)


OUT = io.open(os.path.join(W, "asm_units_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
for smp in sys.argv[1:]:
    if smp not in URL:
        sys.stderr.write("  색인 없음 %s\n" % smp)
        continue
    for hap, anm, u in URL[smp]:
        t0 = time.time()
        tag = "%s_%s" % (smp, hap)
        cnt = find_contig(u, tag)
        t1 = time.time()
        if not cnt:
            OUT.write("\t".join([smp, hap, anm, "-", "-", "-", "-", "-", "NO_MOTIF"]) + NL)
            OUT.flush()
            sys.stderr.write("  %s hap%s 모티프 없음 (%.0fs)\n" % (smp, hap, t1 - t0))
            continue
        fai = os.path.join(W, tag + ".fai")
        gzi = os.path.join(W, tag + ".gzi")
        try:
            subprocess.run(["curl", "-sfL", "-o", fai, u + ".fai"], check=True, timeout=300)
            subprocess.run(["curl", "-sfL", "-o", gzi, u + ".gzi"], check=True, timeout=300)
        except Exception:
            OUT.write("\t".join([smp, hap, anm, "-", "-", "-", "-", "-", "NO_INDEX"]) + NL)
            OUT.flush(); continue
        SZ = {}
        for L in io.open(fai, encoding="utf8"):
            f = L.rstrip(NL).split("\t")
            if len(f) >= 2:
                SZ[f[0]] = int(f[1])
        found = None
        for c, n in cnt.most_common(4):
            if c not in SZ:
                continue
            rr = subprocess.run(["samtools", "faidx", "--fai-idx", fai, "--gzi-idx", gzi,
                                 u, "%s:1-%d" % (c, SZ[c])], capture_output=True, text=True)
            if rr.returncode != 0 or len(rr.stdout) < 10000:
                continue
            reg = os.path.join(W, tag + ".reg.fa")
            io.open(reg, "w", newline=NL).write(rr.stdout)
            paf = subprocess.run(["minimap2", "-x", "asm20", "-t", TH, "--secondary=no",
                                  pp, reg], capture_output=True, text=True)
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
                found = (c, a[1], a[2], seq[a[1]:a[2]], n)
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
            sys.stderr.write("  %s hap%s 탐침 양쪽 없음 [%s] (%.0fs)\n"
                             % (smp, hap, top, time.time() - t0))
            continue
        c, a1, a2, seq, nmot = found
        ok, gaps, nh = units_of(seq, tag)
        OUT.write("\t".join([smp, hap, anm, c, str(a2 - a1), str(nh), str(gaps),
                             str(ok), "OK"]) + NL)
        OUT.flush()
        sys.stderr.write("  %s hap%s  contig %s (모티프 %d) · flank-flank %d bp · 유닛 %d "
                         "(grep %.0fs 총 %.0fs)\n"
                         % (smp, hap, c.split("#")[-1], nmot, a2 - a1, ok,
                            t1 - t0, time.time() - t0))
        sys.stderr.flush()
OUT.close()
print("done")
