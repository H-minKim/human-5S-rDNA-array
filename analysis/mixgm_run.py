# -*- coding: utf-8 -*-
# GM12878 zone-marker counts behind the same short-read alignment gate as K562.
#   Markers : c3_poslist_3zone_GM12878.tsv and c3_poslist_mix_new_GM12878.tsv. Decoys : the GM12878 dispersed 5S consensus,
#   in the Bowtie2 index ref/GM_ARRAY_plus_DECOY.
#   Gate A : the best Bowtie2 score (-a) against the array plus decoys is on the array alone.
#   Gate B : the array score is higher than the best score outside the array in CHM13 (-k 20). Ties are dropped.
#   Each ENCODE track in JOBS is counted behind the gate (filtered) and without it (all).
#   A reference count from GIAB HG001 HiFi reads of the region is aligned with minimap2 to two copies of the unit, without the gate.
import os, sys, io, gzip, subprocess, collections, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
W = os.path.expanduser("~/p11_encode")
ARRAY = ("chr1", 227743598, 228025872)
IDX_G = os.path.join(W, "idx", "chm13")
canon = "".join(l.strip() for l in open(os.path.join(W, "ref", "canon_v2.fa")) if not l.startswith(">")).upper()
L = len(canon)
MK = []
for _pl in ("c3_poslist_3zone_GM12878.tsv", "c3_poslist_mix_new_GM12878.tsv"):
    for ln in open(os.path.join(W, _pl), encoding="utf-8"):
        if not ln.strip(): continue
        v, tag, nu = ln.rstrip(chr(10)).split("\t")
        MK.append((int("".join(c for c in v if c.isdigit())), v[-1], v, tag))
POS = {p for p, _, _, _ in MK}
print("GM12878 markers %d" % len(MK), flush=True)

# ---- the existing combined index is used as it is ; the script stops if it is missing
IDX_AD = os.path.join(W, "ref", "GM_ARRAY_plus_DECOY")
assert os.path.exists(IDX_AD + ".1.bt2"), "GM index is missing"
print("existing GM index:", IDX_AD, flush=True)
B2 = ["-N", "0", "-L", "20", "--score-min", "L,-12,0", "--no-unal", "-p", "12"]


def run(files, lab):
    tmp = os.path.join(W, "gmmix_%s.fq" % lab); nin = 0
    with open(tmp, "w") as o:
        for fp in files:
            with gzip.open(os.path.join(W, fp), "rt") as f:
                buf = []
                for line in f:
                    buf.append(line.rstrip(chr(10)))
                    if len(buf) < 4:
                        continue
                    o.write(chr(10).join(buf) + chr(10)); buf = []; nin += 1
    samA = os.path.join(W, "gmmix_%s_A.sam" % lab)
    subprocess.run(["bowtie2", "-x", IDX_AD, "-U", tmp, "-a"] + B2 + ["-S", samA],
                   stderr=subprocess.DEVNULL, check=False)
    best = collections.defaultdict(lambda: [-999, set()]); rec = {}
    for ln in open(samA):
        if ln[0] == "@":
            continue
        f2 = ln.rstrip(chr(10)).split("\t")
        if int(f2[1]) & 4:
            continue
        AS = next((int(t[5:]) for t in f2[11:] if t.startswith("AS:i:")), -999)
        qn, rf = f2[0], f2[2]
        if AS > best[qn][0]:
            best[qn] = [AS, {rf}]
        elif AS == best[qn][0]:
            best[qn][1].add(rf)
        if rf == "ARRAY" and (qn not in rec or AS >= rec[qn][0]):
            rec[qn] = (AS, int(f2[3]), f2[5], f2[9])
    os.remove(samA)
    passA = {q for q, (a, rs) in best.items() if rs == {"ARRAY"}}
    sub = os.path.join(W, "gmmix_%s_sub.fq" % lab)
    with open(tmp) as f, open(sub, "w") as o:
        buf = []
        for line in f:
            buf.append(line.rstrip(chr(10)))
            if len(buf) < 4:
                continue
            if buf[0][1:].split()[0] in passA:
                o.write(chr(10).join(buf) + chr(10))
            buf = []
    samB = os.path.join(W, "gmmix_%s_B.sam" % lab)
    subprocess.run(["bowtie2", "-x", IDX_G, "-U", sub, "-k", "20"] + B2 + ["-S", samB],
                   stderr=subprocess.DEVNULL, check=False)
    ob = collections.defaultdict(lambda: -999)
    for ln in open(samB):
        if ln[0] == "@":
            continue
        f2 = ln.rstrip(chr(10)).split("\t")
        if int(f2[1]) & 4:
            continue
        if f2[2] == ARRAY[0] and ARRAY[1] <= int(f2[3]) <= ARRAY[2]:
            continue
        AS = next((int(t[5:]) for t in f2[11:] if t.startswith("AS:i:")), -999)
        if AS > ob[f2[0]]:
            ob[f2[0]] = AS
    os.remove(samB); os.remove(sub); os.remove(tmp)
    keep = {q for q in passA if rec.get(q) and rec[q][0] > ob[q]}
    pile = collections.defaultdict(collections.Counter)
    pall = collections.defaultdict(collections.Counter)
    for qn, (AS, rs, cig, seq) in rec.items():
        rp, qp, kp = rs, 0, qn in keep
        for num, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
            num = int(num)
            if op in "M=X":
                for k in range(num):
                    g = ((rp + k - 1) % L) + 1
                    if g in POS:
                        pall[g][seq[qp + k]] += 1
                        if kp:
                            pile[g][seq[qp + k]] += 1
                rp += num; qp += num
            elif op in "DN":
                rp += num
            elif op in "IS":
                qp += num
    print("  %-10s input %9s, aligned %8s, gate A %8s, kept %8s (%.1f%%)"
          % (lab, format(nin, ","), format(len(best), ","), format(len(passA), ","),
             format(len(keep), ","), 100.0 * len(keep) / max(1, len(best))), flush=True)
    return pile, pall


AC = "access/extract/%s.5S.fastq.gz"
GM = "gm12878/extract2/%s.5S.fastq.gz"
JOBS = [
    ("ATAC", [AC % x for x in ("ENCFF948NNH", "ENCFF453JCV", "ENCFF661MZA",
                               "ENCFF733YOH", "ENCFF856KEW", "ENCFF670PJU")]),
    ("DNase", [AC % x for x in ("ENCFF138WVS", "ENCFF895CGN", "ENCFF454HGM",
                                "ENCFF353HJH", "ENCFF556WJB", "ENCFF465NQE")]),
    ("TBP", [GM % x for x in ("ENCFF000WFP", "ENCFF000WFQ")]),
    ("POLR3G", [GM % x for x in ("ENCFF000WCO", "ENCFF000WCU")]),
    ("INPUT", [GM % x for x in ("ENCFF000VWV", "ENCFF000VWD", "ENCFF000VWK")]),
]
res = {}
for lab, fs in JOBS:
    res[lab] = run([f for f in fs if os.path.exists(os.path.join(W, f))], lab)

# ---- HiFi reference count : long reads cannot pass the short-read gate, so they are aligned to two copies of the unit
x2 = os.path.join(W, "ref", "canon_v2_x2_forDISP.fa")
sam = os.path.join(W, "gmmix_hifi.sam")
with open(sam, "w") as o:
    subprocess.run(["minimap2", "-a", "-x", "map-hifi", "-t", "12", "--secondary=no",
                    x2, os.path.join(W, "gm12878", "hifi", "hg001_region.fastq.gz")],
                   stdout=o, stderr=subprocess.DEVNULL)
hp = collections.defaultdict(collections.Counter); nseg = 0
for ln in open(sam):
    if ln[0] == "@":
        continue
    f = ln.rstrip(chr(10)).split("\t")
    if int(f[1]) & 0x904 or f[9] == "*":
        continue
    nseg += 1
    rp, qp, seq = int(f[3]), 0, f[9]
    for num, op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
        num = int(num)
        if op in "M=X":
            for k in range(num):
                g = ((rp + k - 1) % L) + 1
                if g in POS:
                    hp[g][seq[qp + k]] += 1
            rp += num; qp += num
        elif op in "DN":
            rp += num
        elif op in "IS":
            qp += num
os.remove(sam)
print("  %-10s aligned HiFi segments %d (reference count, no gate)" % ("HIFI", nseg), flush=True)

out = os.path.join(W, "READLEVEL_GM_MIX.tsv")
with open(out, "w") as o:
    o.write("label\tfilter\tvariant\ttag\tpos\talt\talt_n\tdepth" + chr(10))
    for lab in res:
        pile, pall = res[lab]
        for p, alt, v, tag in MK:
            for nm, c in (("filtered", pile[p]), ("all", pall[p])):
                o.write("%s\t%s\t%s\t%s\t%d\t%s\t%d\t%d%s"
                        % (lab, nm, v, tag, p, alt, c[alt], sum(c.values()), chr(10)))
    for p, alt, v, tag in MK:
        c = hp[p]
        for nm in ("filtered", "all"):
            o.write("HIFI\t%s\t%s\t%s\t%d\t%s\t%d\t%d%s"
                    % (nm, v, tag, p, alt, c[alt], sum(c.values()), chr(10)))
print("→", out)
