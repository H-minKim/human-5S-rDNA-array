# -*- coding: utf-8 -*-
# . GM12878: the target, input and WGS denominator counts are recounted on the same ruler as the K562 panel G, that is, base quality Q >= 20.
#   The gate, the reference and the counting are a copy of run() in mixgm_run.py with only the quality column (f2[10]) added. Output goes to new files, so READLEVEL_GM_MIX.tsv is not overwritten.
#   Encoding: TBP, POLR3G and INPUT are Phred+64 (Q = ord - 64), WGS is Phred+33 (Q = ord - 33). The quality characters of every file are counted in full and written to the log.
#   Two versions: A = the current alignment (aligned with --phred33, the same reads as the published gate) with the true Q used only in the counting.
#           C = aligned from a fastq rewritten with true Phred+33 quality (ord < 64 becomes Q0) and the true Q, so the ruler matches K562 down to the alignment. This is the sensitivity check.
#   Checks: alt_all/dep_all of version A against READLEVEL_GM_MIX.tsv filtered, WGS against denom_GM_WGS_ERR3239334.tsv, final read count against mix_gm.log.
import os, sys, io, gzip, subprocess, collections, re, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
W = os.path.expanduser("~/p11_encode")
ARRAY = ("chr1", 227743598, 228025872)
IDX_G = os.path.join(W, "idx", "chm13")
IDX_AD = os.path.join(W, "ref", "GM_ARRAY_plus_DECOY")
assert os.path.exists(IDX_AD + ".1.bt2"), "GM 색인이 없다"
canon = "".join(l.strip() for l in open(os.path.join(W, "ref", "canon_v2.fa")) if not l.startswith(">")).upper()
L = len(canon)
MK = []
for _pl in ("c3_poslist_3zone_GM12878.tsv", "c3_poslist_mix_new_GM12878.tsv"):
    for ln in open(os.path.join(W, _pl), encoding="utf-8"):
        if not ln.strip(): continue
        v, tag, nu = ln.rstrip(chr(10)).split("\t")
        MK.append((int("".join(c for c in v if c.isdigit())), v[-1], v, tag))
POS = {p for p, _, _, _ in MK}
QMIN = 20
B2 = ["-N", "0", "-L", "20", "--score-min", "L,-12,0", "--no-unal", "-p", "12"]
print("GM12878 표지 %d 개 · Q>=%d · 색인 %s" % (len(MK), QMIN, IDX_AD), flush=True)
HIST = {}


def fq_hist(fp, n_only=False):
    """Every quality character of one file, which is the evidence for the encoding call."""
    c = collections.Counter(); nr = 0; low = 0; ones = 0
    with gzip.open(fp, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 != 3: continue
            q = line.rstrip(chr(10)); nr += 1; c.update(q)
            o = [ord(x) for x in q]
            if max(o) < 74: low += 1
            if min(o) == 33: ones += 1
    return dict(reads=nr, reads_max_lt74=low, reads_with_ord33=ones, hist={str(ord(k)): v for k, v in sorted(c.items())})


def run(files, lab, off, recode):
    tmp = os.path.join(W, "gmq20_%s.fq" % lab); nin = 0
    with open(tmp, "w") as o:
        for fp in files:
            with gzip.open(os.path.join(W, fp), "rt") as f:
                buf = []
                for line in f:
                    buf.append(line.rstrip(chr(10)))
                    if len(buf) < 4: continue
                    if recode: buf[3] = "".join(chr(max(ord(x) - off, 0) + 33) for x in buf[3])
                    o.write(chr(10).join(buf) + chr(10)); buf = []; nin += 1
    off_cnt = 33 if recode else off      # if the quality was rewritten, the counting uses +33 as well
    samA = os.path.join(W, "gmq20_%s_A.sam" % lab)
    subprocess.run(["bowtie2", "-x", IDX_AD, "-U", tmp, "-a"] + B2 + ["-S", samA], stderr=subprocess.DEVNULL, check=False)
    best = collections.defaultdict(lambda: [-999, set()]); rec = {}
    for ln in open(samA):
        if ln[0] == "@": continue
        f2 = ln.rstrip(chr(10)).split("\t")
        if int(f2[1]) & 4: continue
        AS = next((int(t[5:]) for t in f2[11:] if t.startswith("AS:i:")), -999)
        qn, rf = f2[0], f2[2]
        if AS > best[qn][0]: best[qn] = [AS, {rf}]
        elif AS == best[qn][0]: best[qn][1].add(rf)
        if rf == "ARRAY" and (qn not in rec or AS >= rec[qn][0]):
            rec[qn] = (AS, int(f2[3]), f2[5], f2[9], f2[10])
    os.remove(samA)
    passA = {q for q, (a, rs) in best.items() if rs == {"ARRAY"}}
    sub = os.path.join(W, "gmq20_%s_sub.fq" % lab)
    with open(tmp) as f, open(sub, "w") as o:
        buf = []
        for line in f:
            buf.append(line.rstrip(chr(10)))
            if len(buf) < 4: continue
            if buf[0][1:].split()[0] in passA: o.write(chr(10).join(buf) + chr(10))
            buf = []
    samB = os.path.join(W, "gmq20_%s_B.sam" % lab)
    subprocess.run(["bowtie2", "-x", IDX_G, "-U", sub, "-k", "20"] + B2 + ["-S", samB], stderr=subprocess.DEVNULL, check=False)
    ob = collections.defaultdict(lambda: -999)
    for ln in open(samB):
        if ln[0] == "@": continue
        f2 = ln.rstrip(chr(10)).split("\t")
        if int(f2[1]) & 4: continue
        if f2[2] == ARRAY[0] and ARRAY[1] <= int(f2[3]) <= ARRAY[2]: continue
        AS = next((int(t[5:]) for t in f2[11:] if t.startswith("AS:i:")), -999)
        if AS > ob[f2[0]]: ob[f2[0]] = AS
    os.remove(samB); os.remove(sub); os.remove(tmp)
    keep = {q for q in passA if rec.get(q) and rec[q][0] > ob[q]}
    pall = collections.defaultdict(collections.Counter); phq = collections.defaultdict(collections.Counter)
    mq = collections.Counter(); nd = nt = 0
    for qn, (AS, rs, cig, seq, qual) in rec.items():
        if qn not in keep: continue
        rp, qp = rs, 0
        for num, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
            num = int(num)
            if op in "M=X":
                for k in range(num):
                    g = ((rp + k - 1) % L) + 1
                    if g in POS:
                        ch = qual[qp + k]; nt += 1; mq[ch] += 1; pall[g][seq[qp + k]] += 1
                        if (ord(ch) - off_cnt) >= QMIN: phq[g][seq[qp + k]] += 1
                        else: nd += 1
                rp += num; qp += num
            elif op in "DN": rp += num
            elif op in "IS": qp += num
    print("  %-16s 입력 %s · 붙음 %s · A통과 %s · 최종 %s · 표지관측 %s · Q<%d 탈락 %s (%.1f%%)"
          % (lab, format(nin, ","), format(len(best), ","), format(len(passA), ","), format(len(keep), ","),
             format(nt, ","), QMIN, format(nd, ","), 100.0 * nd / max(1, nt)), flush=True)
    return pall, phq, dict(nin=nin, n_aligned=len(best), pass_a=len(passA), kept=len(keep), n_marker_obs=nt, n_lowq=nd,
                           marker_qual={str(ord(k)): v for k, v in sorted(mq.items())})


GM = "gm12878/extract2/%s.5S.fastq.gz"
JOBS = [("TBP", 64, [GM % x for x in ("ENCFF000WFP", "ENCFF000WFQ")]),
        ("POLR3G", 64, [GM % x for x in ("ENCFF000WCO", "ENCFF000WCU")]),
        ("INPUT", 64, [GM % x for x in ("ENCFF000VWV", "ENCFF000VWD", "ENCFF000VWK")]),
        ("WGS_ERR3239334", 33, ["gmwgs/ERR3239334.5S.fastq.gz"])]
for lab, off, fs in JOBS:
    if off == 64:
        for fp in fs:
            HIST[fp] = fq_hist(os.path.join(W, fp))
            print("  %s reads %d · 최대 ord<74 인 read %d · ord 33 낀 read %d" % (fp, HIST[fp]["reads"], HIST[fp]["reads_max_lt74"], HIST[fp]["reads_with_ord33"]), flush=True)
out = os.path.join(W, "GM_QUALGATE_20261101.tsv")
ST = {}
with open(out, "w") as o:
    o.write("version\tlabel\tvariant\tpos\talt\ttag\talt_all\tdep_all\talt_hq\tdep_hq" + chr(10))
    for ver, jobs in (("A", JOBS), ("C", [j for j in JOBS if j[1] == 64])):
        for lab, off, fs in jobs:
            pall, phq, st = run(fs, lab, off, ver == "C")
            ST["%s|%s" % (ver, lab)] = st
            for p, alt, v, tag in MK:
                ca, ch = pall[p], phq[p]
                o.write("%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d%s" % (ver, lab, v, p, alt, tag, ca[alt], sum(ca.values()), ch[alt], sum(ch.values()), chr(10)))
            o.flush()
json.dump(dict(file_hist=HIST, runs=ST), open(os.path.join(W, "GM_QUALGATE_20261101_stats.json"), "w"), indent=0)
print("→", out, flush=True)
