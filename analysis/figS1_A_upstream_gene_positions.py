# -*- coding: utf-8 -*-
# makeblastdb + blastn → genepos.tsv
#   단계 : S1 A 상류
#   Method document : analysis/METHOD_CONSENSUS_MAP_G8_REBUILD.md

# -*- coding: utf-8 -*-
# 유전자 시작 위치를 그대로 기록한다 — 유닛 수 차이가 진짜인지 경계 못 그은 것인지 가른다
#   왜 : fasta 전체 길이 비교는 flank 를 얼마나 물었나가 섞여 들어간다.
#        array 안쪽만 보려면 첫 유전자 ~ 마지막 유전자 거리를 써야 한다.
#   내는 것 : 파일, 유전자시작 개수, K[0], K[-1], 안쪽거리(K[-1]-K[0]), 전체길이,
#             앞 flank(K[0]), 뒤 flank(len-K[-1]),
#             간격 분포 : 정상(2100-2400) · 큰간격(>2400, 놓친 경계) · 작은간격(<2100)
#             그리고 큰 간격들을 유닛 환산으로 적는다 (4480 ≈ 2 유닛 = 하나 놓침)
import io, os, sys, glob, subprocess
W = os.path.expanduser("~/psw3"); os.makedirs(W, exist_ok=True)
D = os.path.expanduser("~/cmp4d")
R = os.path.expanduser("~/s5s_tuc/ref")
GENE = os.path.join(R, "5S_gene_query.fa")
NL = "\n"; TH = os.environ.get("TH", "3")
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(CMP)[::-1]


def starts(p, tag):
    sq = "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])
    for _ in range(2):
        fa = os.path.join(W, "g_%s.fa" % tag); db = os.path.join(W, "gd_%s" % tag)
        io.open(fa, "w", newline=NL).write(">x" + NL + sq + NL)
        subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                       check=True, capture_output=True)
        r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                            "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                            "-evalue", "1e-5", "-max_target_seqs", "100000000",
                            "-num_threads", TH], capture_output=True, text=True)
        plus = minus = 0; st = set()
        for L in r.stdout.splitlines():
            g = L.split("\t")
            if len(g) < 2:
                continue
            a, b = int(g[0]), int(g[1])
            if a <= b:
                plus += 1
            else:
                minus += 1
            st.add(min(a, b) - 1)
        for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
            try:
                os.remove(db + e)
            except OSError:
                pass
        os.remove(fa)
        if minus > plus:
            sq = rc(sq); continue
        break
    return sorted(st), len(sq)


FL = [os.path.join(D, x.strip()) for x in
      io.open(os.path.expanduser("~/genepos_list.txt")) if x.strip()]
out = io.open(os.path.expanduser("~/genepos.tsv"), "w", encoding="utf8", newline=NL)
out.write("file\tn_gene\tk0\tklast\tinner_bp\ttotal_bp\tflank5\tflank3"
          "\tgap_ok\tgap_big\tgap_small\tbig_gaps\n")
for i, p in enumerate(FL):
    b = os.path.basename(p)
    try:
        K, L = starts(p, "gp%d" % i)
    except Exception:
        out.write("%s\t-1\t-1\t-1\t-1\t-1\t-1\t-1\t0\t0\t0\t-\n" % b); out.flush(); continue
    if not K:
        out.write("%s\t0\t-1\t-1\t-1\t%d\t-1\t-1\t0\t0\t0\t-\n" % (b, L)); out.flush(); continue
    g = [K[j + 1] - K[j] for j in range(len(K) - 1)]
    ok = sum(1 for x in g if 2100 <= x <= 2400)
    big = [x for x in g if x > 2400]
    small = sum(1 for x in g if x < 2100)
    out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n"
              % (b, len(K), K[0], K[-1], K[-1] - K[0], L, K[0], L - K[-1],
                 ok, len(big), small,
                 ",".join("%.2f" % (x / 2240.0) for x in sorted(big, reverse=True)[:8]) or "-"))
    out.flush()
    sys.stderr.write("%d/%d %s\n" % (i + 1, len(FL), b)); sys.stderr.flush()
print("done")
