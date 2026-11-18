#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Is a marker contig significant, that is derived from a real locus, or did it arise by chance? The call is made from the evidence in the reads that make up the contig alone.
#   Absence in the other dataset (PacBio against ONT) is not used as counter-evidence: they come from different groups and different cell samples.
#   A contig is a cluster from mkmarkerasm: at least 2 reads with the marker unit in the same bin. Read classes come from mkreadcls.
#   The evidence, per contig:
#     1) independent reads: the number of distinct reads, that is molecules
#     2) variant base quality: the Phred Q of the base read as the marker allele in each read, from the original FASTQ, and the number of reads at Q >= 20
#                              plus the difference from the median Q at the same canonical position in that read's other units, which shows whether the quality is locally poor for that read
#     3) agreement between reads: a local alignment score between the constituent windows with the marker bin held fixed (similarity 2J-1, gap -0.8), reported as the median over pairs and the fraction of negative pairs
#   The call: significant when 1) is at least 2, 2) gives at least 2 reads at Q >= 20, and 3) has a pair median of at least SMIN = 3; otherwise not significant, with the reason recorded.
#   Usage: python3 mkcontigsig.py <CELL> <TAG> <fa>:<PB|ONT>:<dataset> ... -- <dataset>=<fastq.gz> ...
import sys, os, re, gzip, json, collections, statistics as st, random
argv = sys.argv; cut = argv.index("--"); FQ = dict(x.split("=", 1) for x in argv[cut + 1:]); sys.argv = argv[:cut]
exec(open(os.path.expanduser("~/p11_encode/mkcheck/mkmarkerasm.py")).read().split("RES = []; DIRECT = []")[0])
random.seed(20261014)
SPEC = {x.split(":")[2]: x.split(":")[0] for x in sys.argv[3:]}
# ---- paf index: (read, qs, qe) -> (strand, ts, cigar)
PIDX = {}
for src, fa in SPEC.items():
    for L in open(fa + ".paf"):
        f = L.rstrip("\n").split("\t")
        cg = [x for x in f[12:] if x.startswith("cg:Z:")]
        if cg: PIDX[(f[0], int(f[2]), int(f[3]))] = (f[4], int(f[7]), cg[0][5:])
# ---- FASTQ quality, for reads with a unit chain only
need = collections.defaultdict(set)
for q, r in READS.items(): need[r["src"]].add(q.split("|", 1)[1])
QUAL = {}
for src, path in FQ.items():
    with gzip.open(os.path.expanduser(path), "rt") as f:
        while True:
            h = f.readline()
            if not h: break
            f.readline(); f.readline(); qq = f.readline().rstrip("\n")
            nm = h[1:].split()[0]
            if nm in need[src]: QUAL[(src, nm)] = qq
log("quality values read", {s: sum(1 for k in QUAL if k[0] == s) for s in FQ})
def base_q(q, r, u, tpos):
    """Q of the base at canonical coordinate tpos in unit u of read q (None if absent)"""
    rec = PIDX.get((q, u["qs"], u["qe"]))
    qq = QUAL.get((r["src"], q.split("|", 1)[1]))
    if not rec or qq is None: return None
    stn, ts, cg = rec; rp, tp = 0, ts
    for n, op in re.findall(r"(\d+)([MIDNSH=X])", cg):
        n = int(n)
        if op in "M=X":
            if tp < tpos <= tp + n:
                off = rp + (tpos - tp - 1); ro = u["qs"] + off if stn == "+" else u["qe"] - 1 - off
                return ord(qq[ro]) - 33 if 0 <= ro < len(qq) else None
            rp += n; tp += n
        elif op == "I": rp += n
        elif op in "DN":
            if tp < tpos <= tp + n: return None
            tp += n
    return None
# ---- contigs and read classes
G5 = json.load(open("%s/%s_markerasm_v5.json" % (W, TAG)))["results"]
RC = {(o["variant"], o["read"], o["k"]): o for o in json.load(open("%s/%s_readcls_v4.json" % (W, TAG)))}
VK = {m["variant"]: (p, a) for (p, a), m in MK.items()}
OUT = []
for gi, g in enumerate(G5):
    reads = {q for q, k in g["idx"]}
    if len(reads) < 2: continue
    v = g["variant"]; p, a = VK[v]; mz = g["mzone"]
    items = []
    for q, k in g["idx"]:
        r = READS[q]; u = r["ch"][k]
        qa = base_q(q, r, u, p)
        qref = [base_q(q, r, w, p) for i, w in enumerate(r["ch"]) if i != k and (p, a) not in w["raw"] and w["lo"] <= p <= w["hi"]]
        qref = [x for x in qref if x is not None]
        o = RC.get((v, q, k), {})
        lo, hi = max(0, k - HALF), min(len(r["ch"]), k + HALF + 1)
        items.append(dict(read=q, plat=r["plat"], k=k, q_alt=qa, q_ref_med=st.median(qref) if qref else None, cls=o.get("cls", ""), zone=o.get("zone", ""),
                          win=r["ch"][lo:hi], c=k - lo))
    # 3) pairwise alignment, up to 60 pairs
    prs = [(i, j) for i in range(len(items)) for j in range(i + 1, len(items)) if items[i]["read"] != items[j]["read"]]
    random.shuffle(prs); sc = []
    for i, j in prs[:60]:
        A, Bw = items[i], items[j]
        sc.append(center_profile(A["win"], A["c"], Bw["win"], sim_f)[Bw["c"]])
    # contig class: A if any constituent read carries the original locus
    cls = {it["cls"][:2].strip() for it in items}
    if "A" in cls: gcls = "A 원래 자리"
    else:
        zs = set()
        for it in items:
            if it["zone"] and it["zone"] != "미정":
                for z in it["zone"].split("|"): zs.add({"notD": "M", "notP": "M"}.get(z, z) if mz == "M" else z)
        zs -= {"D?", "P?"}
        gcls = "같은 구역 · 다른 자리/조립본 밖/모호" if (not zs or zs <= {mz} or mz in zs) else "다른 구역"
    # independent reads, by read name, and for PacBio down to the ZMW
    mol = {re.sub(r"/ccs.*$", "", it["read"]) for it in items}
    qa = [it["q_alt"] for it in items if it["q_alt"] is not None]
    # where one read carries several units, the highest Q per read is taken
    best = {}
    for it in items:
        if it["q_alt"] is not None: best[it["read"]] = max(best.get(it["read"], -1), it["q_alt"])
    nq20 = sum(1 for x in best.values() if x >= 20)
    dq = [it["q_alt"] - it["q_ref_med"] for it in items if it["q_alt"] is not None and it["q_ref_med"] is not None]
    pmed = st.median(sc) if sc else None
    why = []
    if len(mol) < 2: why.append("독립 read <2")
    if nq20 < 2: why.append("변이 염기 Q>=20 read %d 개" % nq20)
    if pmed is None or pmed < SMIN: why.append("read 끼리 맞음 약함(쌍 중앙 %s)" % (None if pmed is None else round(pmed, 2)))
    OUT.append(dict(contig=gi, variant=v, mzone=mz, cls=gcls, n_read=len(reads), n_mol=len(mol), plat=dict(collections.Counter(READS[q]["plat"] for q in reads)),
                    q_alt=sorted(qa), n_q20=nq20, dq_med=st.median(dq) if dq else None, pair_med=None if pmed is None else round(pmed, 2),
                    pair_neg=round(sum(1 for x in sc if x < 0) / len(sc), 2) if sc else None, n_pair=len(sc),
                    sig="유의" if not why else "유의 아님", why="; ".join(why),
                    per_read=[dict(read=it["read"], plat=it["plat"], k=it["k"], q_alt=it["q_alt"], q_ref_med=it["q_ref_med"], cls=it["cls"], zone=it["zone"]) for it in items]))
json.dump(OUT, open("%s/%s_contigsig.json" % (W, TAG), "w"), ensure_ascii=False, indent=0)
c = collections.Counter((o["cls"], o["sig"]) for o in OUT)
log("contig", len(OUT), {"%s|%s" % k: n for k, n in sorted(c.items())})
