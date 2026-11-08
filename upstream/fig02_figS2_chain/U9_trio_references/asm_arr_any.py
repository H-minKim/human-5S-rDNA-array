# -*- coding: utf-8 -*-
# Copy of asm_arr_all.py with only input and output changed (2026-09-18, extraction of trio canonical candidates)
#   original: ~/asm_arr_all.py (HPRC release2 index CSV, faidx with remote .fai/.gzi)
#   changed (unavoidable): the public assemblies (1kGP-ONT hifiasm asm_1/2.fa, HPRC year1 genbank .fa.gz) have
#     no remote .fai/.gzi, and the node disk (64 GB) cannot hold 3 GB x 20
#     -> ① finding motif contigs is the same streaming grep as the original
#       ② only the extraction of contig sequences uses awk in a second streaming pass to save just those contigs, instead of faidx
#   unchanged: motifs (4 5S gene 25-mers + reverse complements), top 4 contigs, probes minimap2 asm20 --secondary=no
#     , probe alignment 20 kb or more, both PROX50k/DIST50k, array = between the inner ends of the probes, unit = 2100~2400 between blastn HSP starts
#     , output format (arr_<tag>.fa header, units_<tag>.txt, asmarr_<TAG>.tsv 9 columns)
#   input: list TSV (sample, hap, name, url), env vars W (output folder) TAG TH
import io, os, sys, subprocess, collections, time
W = os.path.expanduser(os.environ.get("W", "~/asmarr_trio"))
os.makedirs(W, exist_ok=True)
GENE = os.path.expanduser("~/s5s_tuc/ref/5S_gene_query.fa")
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


def units_of(seq, tag):
    """Chain of unit start coordinates (0-based within the array), same as the original"""
    fa = os.path.join(W, "a_%s.fa" % tag)
    io.open(fa, "w", newline=NL).write(">arr" + NL + seq + NL)
    db = os.path.join(W, "d_%s" % tag)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db], check=True, capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send", "-dust", "no", "-word_size", "9",
                        "-evalue", "1e-5", "-max_target_seqs", "100000000",
                        "-num_threads", TH], capture_output=True, text=True, check=True)
    st = set(); pl = mi = 0
    for L in r.stdout.splitlines():
        f = L.split("\t")
        if len(f) >= 2:
            a, b = int(f[0]), int(f[1]); st.add(min(a, b))
            if a < b: pl += 1
            else: mi += 1
    K = sorted(x - 1 for x in st)
    ok = sum(1 for i in range(len(K) - 1) if 2100 <= K[i + 1] - K[i] <= 2400)
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try: os.remove(db + e)
        except OSError: pass
    os.remove(fa)
    return ok, K, (1 if pl >= mi else -1)


def stream(u):
    src = u.replace("s3://human-pangenomics/", "https://s3-us-west-2.amazonaws.com/human-pangenomics/") \
           .replace("s3://1000g-ont/", "https://1000g-ont.s3.amazonaws.com/")
    return "curl -sfL " + repr(src) + (" | zcat" if src.endswith(".gz") else "")


OUT = io.open(os.path.join(W, "asmarr_%s.tsv" % TAG), "a", encoding="utf8", newline=NL)
for L0 in io.open(sys.argv[1], encoding="utf8"):
    if not L0.strip() or L0.startswith("#"): continue
    smp, hap, anm, u = L0.rstrip(NL).split("\t")[:4]
    tag = "%s_%s" % (smp, hap)
    if os.path.exists(os.path.join(W, "arr_%s.fa" % tag)): continue
    t0 = time.time()
    # ① count motif contigs (same as the original)
    pat = "|".join(MOT)
    pr = subprocess.Popen(stream(u) + " | grep -E " + repr("^>|" + pat), shell=True, stdout=subprocess.PIPE,
                          text=True, stderr=subprocess.DEVNULL)
    cur = None; cnt = collections.Counter()
    for L in pr.stdout:
        if L[0] == ">": cur = L[1:].split()[0]
        elif cur: cnt[cur] += 1
    pr.wait()
    if not cnt:
        OUT.write("\t".join([smp, hap, anm, "-", "-", "-", "-", "-", "NO_MOTIF"]) + NL); OUT.flush(); continue
    top = [c for c, n in cnt.most_common(4)]
    # ② save only those contigs by streaming (replaces faidx)
    keep = os.path.join(W, tag + ".keep.txt"); io.open(keep, "w", newline=NL).write(NL.join(top) + NL)
    regall = os.path.join(W, tag + ".top.fa")
    aw = "awk 'NR==FNR{k[$1]=1;next} /^>/{split(substr($0,2),a,\" \"); p=(a[1] in k)} p' " + repr(keep) + " -"
    subprocess.run(stream(u) + " | " + aw + " > " + repr(regall), shell=True, check=True)
    SEQ = {}; cur = None; buf = []
    for L in io.open(regall, encoding="utf8"):
        if L[0] == ">":
            if cur: SEQ[cur] = "".join(buf)
            cur = L[1:].split()[0]; buf = []
        else: buf.append(L.strip())
    if cur: SEQ[cur] = "".join(buf)
    os.remove(regall); os.remove(keep)
    found = None
    for c in top:
        if c not in SEQ or len(SEQ[c]) < 10000: continue
        reg = os.path.join(W, tag + ".reg.fa")
        io.open(reg, "w", newline=NL).write(">" + c + NL + SEQ[c] + NL)
        paf = subprocess.run(["minimap2", "-x", "asm20", "-t", TH, "--secondary=no", PP, reg], capture_output=True, text=True)
        d = {}
        for LL in paf.stdout.splitlines():
            f = LL.split("\t")
            if len(f) < 12 or int(f[9]) < 20000: continue
            if f[5] not in d or int(f[9]) > d[f[5]][2]: d[f[5]] = (int(f[2]), int(f[3]), int(f[9]))
        os.remove(reg)
        if "PROX50k" in d and "DIST50k" in d:
            a = sorted([d["PROX50k"][0], d["PROX50k"][1], d["DIST50k"][0], d["DIST50k"][1]])
            seq = SEQ[c]
            pmid = (d["PROX50k"][0] + d["PROX50k"][1]) // 2
            side = "PROX_FIRST" if pmid < a[1] else "DIST_FIRST"
            found = (c, a[1], a[2], seq[a[1]:a[2]], side)
            break
    if not found:
        t2 = ",".join("%s:%d" % (c, n) for c, n in cnt.most_common(2))
        OUT.write("\t".join([smp, hap, anm, t2, "-", "-", "-", "-", "NO_BOTH_PROBE"]) + NL); OUT.flush()
        sys.stderr.write("  %s h%s 탐침 양쪽 없음 (%.0fs)\n" % (smp, hap, time.time() - t0)); continue
    c, a1, a2, seq, side = found
    ok, K, sg = units_of(seq, tag)
    io.open(os.path.join(W, "arr_%s.fa" % tag), "w", newline=NL).write(
        ">%s|%s|%d|%d|%s|%s%s%s%s" % (tag, c, a1, a2, side, "+" if sg > 0 else "-", NL, seq, NL))
    io.open(os.path.join(W, "units_%s.txt" % tag), "w", newline=NL).write(NL.join(str(x) for x in K) + NL)
    OUT.write("\t".join([smp, hap, anm, c, str(a1), str(a2), str(a2 - a1), str(ok), "OK"]) + NL); OUT.flush()
    sys.stderr.write("  %s h%s  %s · %d-%d (%d bp) · 유닛 %d · %s %s (%.0fs)\n"
                     % (smp, hap, c[:24], a1, a2, a2 - a1, ok, side, "+" if sg > 0 else "-", time.time() - t0))
OUT.close()
print("done")
