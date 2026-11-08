# -*- coding: utf-8 -*-
# 전 게놈 read 스트림에서 5S read 만 건진다 — BBDuk k=25 hdist=2 관문 + 유닛/유전자 질의. 이 글이 만든 `*.hits.fastq.gz` 가 hifiasm 의 입력이다
#   단계 : ① read 추출 관문

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ONT 추출 v5 — BBDuk k-mer 단일 관문 (tRNA 제거 · 스트림 blastn 제거).

바뀐 이유 (전부 실측, 같은 383MB 입력)
  · 정답 = 현행 파이프라인이 ARRAY 로 판정해 보관한 read 2,262 개.
      BBDuk k31 hdist1 (판별세그먼트 NTS)  188 MB/s · 재현율 1.0000 · 위양성 8
      blastn ws=16                          8 MB/s · 재현율 1.0000 · 위양성 728
      megablast ws=28                      14 MB/s · 재현율 0.9881 (놓침 27)
    → 스트리밍 관문을 BBDuk 로 바꾸면 재현율 손실 없이 23 배 빨라진다.
  · the earlier note "k-mer selection unsuitable" refers to **dispersed pseudogenes** diverged by 76-92%.
    지금 보관 대상인 array read 는 정본과 거의 동일해 k-mer 가 완벽하게 듣는다.
  · blastn 은 통과한 read(전체의 ~0.3%)에만 돌려 5s.tsv 필드를 채운다 — 비용 무시 가능.
  · tRNA is excluded. The MM and ML tags are preserved.
"""
import sys, os, subprocess, time, collections, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--nts", required=True, help="5S_unit_discseg.fa (판별세그먼트)")
ap.add_argument("--gene", required=True, help="5S 유전자 fasta (양성대조 겸 후처리 질의)")
ap.add_argument("--gene-db", required=True)
ap.add_argument("--bbduk", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--chunk", type=int, default=16000)
ap.add_argument("--threads", type=int, default=6)
ap.add_argument("--k", type=int, default=25)      # ONT 실측: k25 hdist2 에서 span>500 재현율 1.0000
ap.add_argument("--hdist", type=int, default=2)   # (k31 h1 은 0.9158 — 발산한 짧은 NTS 를 놓친다)
ap.add_argument("--tmp", default="/tmp")
a = ap.parse_args()

for f in (a.nts, a.gene, a.bbduk, a.gene_db + ".nin"):
    if not os.path.exists(f):
        sys.stderr.write("FATAL: missing %s\n" % f); raise SystemExit(9)
if subprocess.run(["bash", "-lc", "command -v blastn"], capture_output=True).returncode != 0:
    sys.stderr.write("FATAL: blastn not on PATH\n"); raise SystemExit(9)

NTSSEQ = "".join(l.strip() for l in open(a.nts) if not l.startswith(">"))
POSCTL = "POSCTL_NTS"

fq = open(a.out + ".hits.fastq", "w")
t0 = time.time()
n_read = n_bp = nhit = nchunk = 0
chunk = []


def flush(chunk):
    global nhit, nchunk
    if not chunk: return
    nchunk += 1
    tag = "%d_%d" % (os.getpid(), nchunk % 2)
    fa = os.path.join(a.tmp, "v5_%s.fa" % tag)
    with open(fa, "w") as fh:
        for rid, seq, qual, cmt in chunk:
            fh.write(">%s\n%s\n" % (rid, seq))
        fh.write(">%s\n%s\n" % (POSCTL, NTSSEQ))
    tout = os.path.join(a.tmp, "v5m_%s.fa" % tag)
    p = subprocess.Popen(["bash", a.bbduk, "-Xmx3g", "in=" + fa, "ref=" + a.nts,
                          "k=%d" % a.k, "hdist=%d" % a.hdist, "mkh=2", "rcomp=t", "outm=" + tout,
                          "overwrite=t", "threads=%d" % a.threads],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _, err = p.communicate()
    if p.returncode != 0:
        sys.stderr.write("FATAL bbduk: %s\n" % (err or b"")[-400:].decode("utf8", "replace"))
        raise SystemExit(9)
    keep = set()
    if os.path.exists(tout):
        for L in open(tout):
            if L.startswith(">"): keep.add(L[1:].strip().split()[0])
    # 양성대조: 판별세그먼트 자신이 안 잡히면 참조/설정이 깨진 것이다
    if POSCTL not in keep:
        sys.stderr.write("FATAL: 양성대조 미검출 (chunk %d)\n" % nchunk); raise SystemExit(9)
    keep.discard(POSCTL)
    for rid, seq, qual, cmt in chunk:
        if rid in keep:
            nhit += 1
            fq.write("@%s%s\n%s\n+\n%s\n" % (rid, (" " + cmt) if cmt else "", seq, qual))
    for f_ in (fa, tout):
        try: os.remove(f_)
        except OSError: pass


buf = []
for line in sys.stdin:
    buf.append(line.rstrip("\n"))
    if len(buf) == 4:
        pp = buf[0][1:].split(None, 1)
        chunk.append((pp[0], buf[1], buf[3], pp[1] if len(pp) > 1 else ""))
        n_read += 1; n_bp += len(buf[1]); buf = []
        if len(chunk) >= a.chunk:
            flush(chunk); chunk = []
            el = time.time() - t0
            sys.stderr.write("[%6.0fs] read=%d (%.2f Gbp) hit=%d  %.0f Mbp/s\n"
                             % (el, n_read, n_bp / 1e9, nhit, n_bp / 1e6 / max(el, 1)))
            sys.stderr.flush()
flush(chunk)
fq.close()

# ---- 통과 read 에만 blastn (5s.tsv 필드용) ----
hf = a.out + ".hits.fastq"
fa = a.out + ".hits.fa"
with open(fa, "w") as fh, open(hf) as f:
    i = 0
    for L in f:
        i += 1
        if i % 4 == 1: fh.write(">" + L[1:].split()[0] + "\n")
        elif i % 4 == 2: fh.write(L)
o = a.out + ".gene.tsv"
r = subprocess.run(["blastn", "-task", "blastn", "-word_size", "9", "-dust", "no",
                    "-evalue", "1e-3", "-db", a.gene_db, "-query", fa,
                    "-num_threads", str(a.threads), "-max_target_seqs", "10000",
                    "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send",
                    "-out", o])
if r.returncode != 0:
    sys.stderr.write("FATAL: 후처리 blastn 실패\n"); raise SystemExit(9)
hits = collections.defaultdict(list)
for L in open(o):
    f_ = L.rstrip("\n").split("\t")
    if int(f_[3]) < 40: continue
    hits[f_[0]].append((float(f_[2]), int(f_[3]), int(f_[4]), int(f_[5])))
mmflag = {}
lens = {}
i = 0
with open(hf) as f:
    for L in f:
        i += 1
        if i % 4 == 1:
            pp = L[1:].rstrip("\n").split(None, 1)
            cur = pp[0]; mmflag[cur] = "Y" if (len(pp) > 1 and "MM:Z" in pp[1]) else "N"
        elif i % 4 == 2:
            lens[cur] = len(L.strip())
with open(a.out + ".5s.tsv", "w") as t5:
    t5.write("read_id\tread_len\tn_hsp\tn_gene\tbest_pid\tspacings\thas_MM\n")
    for rid in lens:
        v = hits.get(rid, [])
        if v:
            st = sorted(min(x[2], x[3]) for x in v)
            m = [st[0]]
            for s_ in st[1:]:
                if s_ - m[-1] >= 150: m.append(s_)
            sp = [m[j + 1] - m[j] for j in range(len(m) - 1)]
            t5.write("%s\t%d\t%d\t%d\t%.1f\t%s\t%s\n"
                     % (rid, lens[rid], len(v), len(m), max(x[0] for x in v),
                        ",".join(map(str, sp)) if sp else "-", mmflag[rid]))
        else:
            t5.write("%s\t%d\t0\t0\t0.0\t-\t%s\n" % (rid, lens[rid], mmflag[rid]))
for f_ in (fa, o):
    try: os.remove(f_)
    except OSError: pass
el = time.time() - t0
sys.stderr.write("\n=== 완료 ===\nread=%d  %.3f Gbp  %.0fs  %.0f Mbp/s  hit=%d (%.4f%%)  chunk=%d\n"
                 % (n_read, n_bp / 1e9, el, n_bp / 1e6 / max(el, 1), nhit,
                    100.0 * nhit / max(n_read, 1), nchunk))
MIN_READ_FOR_ZERO_HIT = 50000
if n_read >= MIN_READ_FOR_ZERO_HIT and nhit == 0:
    sys.stderr.write("FATAL: read %d 인데 히트 0 — 조용한 실패\n" % n_read); raise SystemExit(9)
