# -*- coding: utf-8 -*-
# K562 zone-marker counts behind the short-read alignment gate, over all kept bases and over bases of quality 20 or more.
#   The reads of each ENCODE track in JOBS are aligned with Bowtie2.
#   Gate A : the best score (-a) against the array plus the dispersed 5S decoys is on the array alone.
#   Gate B : the array score is higher than the best score outside the array in CHM13 (-k 20). Ties are dropped.
#   At each marker position the marker allele and the depth are counted over the kept reads, once for all bases and once
#   for bases with Q >= 20 (the K562 files are Phred+33, Q = ord - 33).
import os, sys, io, gzip, subprocess, collections, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
W = os.path.expanduser("~/p11_encode")
ARRAY = ("chr1", 227743598, 228025872)
IDX_G = os.path.join(W, "idx", "chm13")
IDX_AD = os.path.join(W, "ref", "ARRAY_plus_DECOY")
if not os.path.exists(IDX_AD + ".1.bt2"):
    for cand in ("ARRAY_plus_DECOY", "K562_ARRAY_plus_DECOY", "GM_ARRAY_plus_DECOY"):
        if os.path.exists(os.path.join(W, "ref", cand + ".1.bt2")):
            IDX_AD = os.path.join(W, "ref", cand); break
print("combined reference index: %s" % os.path.basename(IDX_AD), flush=True)
canon = "".join(l.strip() for l in open(os.path.join(W, "ref", "canon_v2.fa")) if not l.startswith(">")).upper()
L = len(canon)
PL = os.path.join(W, "c3_poslist_v4.tsv")
if not os.path.exists(PL): PL = os.path.join(W, "c3_poslist_3zone.tsv")
MK = []
for ln in open(PL, encoding="utf-8"):
    parts = ln.rstrip(chr(10)).split(chr(9))
    if len(parts) < 2: continue
    v, tag = parts[0], parts[1]
    nu = int(parts[2]) if len(parts) > 2 and parts[2].strip().isdigit() else 0
    MK.append((int("".join(c for c in v if c.isdigit())), v[-1], v, tag, nu))
POS = {p for p, _, _, _, _ in MK}
QMIN = 20
B2 = ["-N", "0", "-L", "20", "--score-min", "L,-12,0", "--no-unal", "-p", "12"]
print("position list %s, markers %d, quality Q>=%d" % (os.path.basename(PL), len(MK), QMIN), flush=True)

RG = "regate/%s.5S.fastq.gz"
JOBS = [
    ("WGSNAT",  ["cat_WGS.fq.gz"]),
    ("TBP",     [RG % x for x in ("ENCFF000ZBZ", "ENCFF000ZCA")]),
    ("POLR3A",  [RG % x for x in ("ENCFF000YYV", "ENCFF000YYW")]),
    ("BDP1",    [RG % x for x in ("ENCFF000YGF", "ENCFF000YGH", "ENCFF000YGI")]),
    ("BRF1",    [RG % x for x in ("ENCFF000YGP", "ENCFF000YGR", "ENCFF000YGS")]),
    ("GTF3C2",  [RG % x for x in ("ENCFF000ZBY", "ENCFF000ZCN", "ENCFF000ZCO")]),
    ("ATAC",    [RG % x for x in ("ENCFF391BFJ", "ENCFF350ZZR", "ENCFF186CQZ", "ENCFF440UAD")]),
    ("SNAPC5",  [RG % x for x in ("ENCFF921IAW", "ENCFF056GLE", "ENCFF044LPQ", "ENCFF881ZCU")]),
    ("INPUT",   [RG % x for x in ("ENCFF000YRC", "ENCFF000YRN", "ENCFF000YQC")]),
]
OUT = open(os.path.join(W, "K562_QUALGATE.tsv"), "w")
OUT.write("label\tvariant\tpos\talt\tn_unit\ttag\talt_all\tdep_all\talt_hq\tdep_hq" + chr(10))

for lab, files in JOBS:
    ex = [f for f in files if os.path.exists(os.path.join(W, f))]
    if not ex:
        print("  %-8s no input, skipped" % lab, flush=True); continue
    tmp = os.path.join(W, "jb_%s.fq" % lab); nin = 0
    with open(tmp, "w") as o:
        for fp in ex:
            with gzip.open(os.path.join(W, fp), "rt") as f:
                buf = []
                for line in f:
                    buf.append(line.rstrip(chr(10)))
                    if len(buf) < 4: continue
                    o.write(chr(10).join(buf) + chr(10)); buf = []; nin += 1
    samA = os.path.join(W, "jb_A.sam")
    subprocess.run(["bowtie2", "-x", IDX_AD, "-U", tmp, "-a"] + B2 + ["-S", samA], stderr=subprocess.DEVNULL, check=False)
    best = collections.defaultdict(lambda: [-999, set()]); rec = {}
    for ln in open(samA):
        if ln[0] == "@": continue
        f2 = ln.rstrip(chr(10)).split(chr(9))
        if int(f2[1]) & 4: continue
        AS = next((int(t[5:]) for t in f2[11:] if t.startswith("AS:i:")), -999)
        qn, rf = f2[0], f2[2]
        if AS > best[qn][0]: best[qn] = [AS, {rf}]
        elif AS == best[qn][0]: best[qn][1].add(rf)
        if rf == "ARRAY" and (qn not in rec or AS >= rec[qn][0]):
            rec[qn] = (AS, int(f2[3]), f2[5], f2[9], f2[10])
    os.remove(samA)
    passA = {q for q, (a, rs) in best.items() if rs == {"ARRAY"}}
    sub = os.path.join(W, "jb_sub.fq")
    with open(tmp) as f, open(sub, "w") as o:
        buf = []
        for line in f:
            buf.append(line.rstrip(chr(10)))
            if len(buf) < 4: continue
            if buf[0][1:].split()[0] in passA: o.write(chr(10).join(buf) + chr(10))
            buf = []
    samB = os.path.join(W, "jb_B.sam")
    subprocess.run(["bowtie2", "-x", IDX_G, "-U", sub, "-k", "20"] + B2 + ["-S", samB], stderr=subprocess.DEVNULL, check=False)
    ob = collections.defaultdict(lambda: -999)
    for ln in open(samB):
        if ln[0] == "@": continue
        f2 = ln.rstrip(chr(10)).split(chr(9))
        if int(f2[1]) & 4: continue
        if f2[2] == ARRAY[0] and ARRAY[1] <= int(f2[3]) <= ARRAY[2]: continue
        AS = next((int(t[5:]) for t in f2[11:] if t.startswith("AS:i:")), -999)
        if AS > ob[f2[0]]: ob[f2[0]] = AS
    os.remove(samB); os.remove(sub); os.remove(tmp)
    keep = {q for q in passA if rec.get(q) and rec[q][0] > ob[q]}
    pall = collections.defaultdict(collections.Counter); phq = collections.defaultdict(collections.Counter)
    nd = nt = 0
    for qn, (AS, rs, cig, seq, qual) in rec.items():
        if qn not in keep: continue
        rp, qp = rs, 0
        for num, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
            num = int(num)
            if op in "M=X":
                for k in range(num):
                    g = ((rp + k - 1) % L) + 1
                    if g in POS:
                        nt += 1; pall[g][seq[qp + k]] += 1
                        if (ord(qual[qp + k]) - 33) >= QMIN: phq[g][seq[qp + k]] += 1
                        else: nd += 1
                rp += num; qp += num
            elif op in "DN": rp += num
            elif op in "IS": qp += num
    print("  %-8s input %s, aligned %s, kept %s (%.1f%%), marker observations %s, Q<20 dropped %s (%.1f%%)"
          % (lab, format(nin, ","), format(len(best), ","), format(len(keep), ","),
             100.0*len(keep)/max(1, len(best)), format(nt, ","), format(nd, ","), 100.0*nd/max(1, nt)), flush=True)
    for p, alt, v, tg, nu in MK:
        ca, ch = pall[p], phq[p]
        OUT.write("%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d%s"
                  % (lab, v, p, alt, nu, tg, ca[alt], sum(ca.values()), ch[alt], sum(ch.values()), chr(10)))
    OUT.flush()
OUT.close()
print("→ ~/p11_encode/K562_QUALGATE.tsv", flush=True)
