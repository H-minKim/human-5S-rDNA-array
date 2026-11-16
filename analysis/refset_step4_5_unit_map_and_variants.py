# -*- coding: utf-8 -*-
# The script that produced the ruler the figures use: g8canonv11f1/CONSENSUS_MAP_G8_full.tsv (bb68fe3a2a23) and CONSENSUS_SUM_G8_full.tsv (4c8252462824), with all 616 slots re-cut.
#   Steps 4 and 5: unify the orientation, cut the units, call the variants.

# -*- coding: utf-8 -*-
# The gene boundary.
#   A blastn gene HSP stops one or two bases short where the last base of the gene carries a variant (measured: 50 of 57,878 HSPs over the 616 canonical sequences).
#   The boundary is therefore set by extending the HSP to the full gene length using the query coordinates.
#     plus strand: start = sstart - (qstart - 1), end = send + (119 - qend)
#     minus strand: start = send - (119 - qend), end = sstart + (qstart - 1)
# Inputs and settings, following METHOD_CONSENSUS_MAP_G8_REBUILD.md:
#   1) SEQF: hap_seq_node, from the environment variable HSN. If it is unset the script stops.
#   2) .fa: the folder CMPD.
#   3) The orientation comes from the majority strand alone (method section 4). The revcomp column of SUM records whether the majority strand required a flip.
#   4) W: the working folder.
#   5) The output carries a slot_space column (the assembly, method section 3-3).
#   6) A missing canonical sequence, a missing file, and a minus-strand majority persisting after the flip are counted rather than passed over silently.
# Consensus map G8: the assembly arrays are cut again from the beginning.
#
#  Why the cut is made in the canonical phase:
#    G6 cuts the unit at the gene start, giving [gene 119 bp][AluY .. NTS2].
#         The canonical unit is [AluY .. NTS2][gene], so the two are rotated by 119 bp.
#         minimap2 cannot wrap around, so in that frame canonical positions 2122-2240 are never covered,
#         the variants stop at 2101 and no variant is called in the 5S gene.
#
#  G8 cuts in the canonical phase from the start.
#    The canonical unit is [AluY 1-318][NTS1 319-1178][GA-rich 1179-1551][NTS2 1552-2121][gene 2122-2240].
#    So the canonical unit begins immediately after the gene ends,
#    and the unit boundary is therefore taken at the end of the gene HSP:
#         unit i = (gene_end[i] + 1) .. gene_end[i+1].
#       No rotation is then needed: the unit as cut is already in the canonical phase.
#
#  The flank has already been handled upstream.
#    The .fa files are flank-to-flank sequences cut with the 50 kb flank probes,
#    so the edges of a .fa are the junctions. The flank is not touched again here.
#
#  The first partial unit is excluded.
#    From the junction to the end of the first gene is not a complete unit (median 1,516 bp measured).
#    The remainder at the far end of the array is the same. Both are marked PARTIAL and
#    no variants are called for them.
#
#  The checks are made against sequence:
#    1) the unit sequence is reconstructed from the variant list and compared with the original character by character;
#    2) the alignment is made against a context-flanked reference with slack, so that it covers the ends of the unit
#       (minimap2 does not extend through a terminal mismatch and so swallows that variant;
#        measured, 237 of 400 such cases are real variants);
#    3) the unit boundaries are confirmed from the gene HSPs, never by dividing bp by 2240.
import io
import os
import re
import sys
import subprocess
import collections

NL = "\n"
W = os.path.expanduser(os.environ["WDIR"])
os.makedirs(W, exist_ok=True)
R = os.path.expanduser("~/s5s_tuc/ref")
CANON = os.path.join(R, "unit_canon.fa")
GENE = os.path.join(R, "5S_gene_query.fa")
TH = os.environ.get("TH", "4")
TAG = os.environ.get("TAG", "x")
PRESET = os.environ.get("PRESET", "asm20")
LO = int(os.environ.get("LO", "2100"))          # lower bound on an acceptable unit gap
HI = int(os.environ.get("HI", "2400"))          # upper bound
FRAG = int(os.environ.get("FRAG", "5"))         # an extension larger than this marks the HSP as a fragment
CMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
CS = re.compile(r"([:*+\-=])([A-Za-z0-9]+)")
NUM = re.compile(r"\d+")


def rc(s):
    return s.translate(CMP)[::-1]


def rd(p):
    return "".join(x.strip() for x in io.open(p, encoding="utf8").read().split(NL)[1:])


CANSEQ = rd(CANON).upper()
PAD = int(os.environ.get("PAD", "100"))           # the slack left on each side when cutting a unit
CTX = int(os.environ.get("CTX", "200"))           # the context on each side of the canonical unit; the array being tandem, these are the end and the start of the same unit
REFSEQ = CANSEQ[-CTX:] + CANSEQ + CANSEQ[:CTX]
CANREF = os.path.join(W, "canon_ctx.fa")
io.open(CANREF, "w", newline=NL).write(">canon_ctx" + NL + REFSEQ + NL)
if len(CANSEQ) != 2240:
    sys.stderr.write("정본유닛 길이가 2240 이 아니다: %d%s" % (len(CANSEQ), NL))
    raise SystemExit(1)
GENESEQ = rd(GENE).upper()
if len(GENESEQ) != 119:
    sys.stderr.write("유전자 질의 길이가 119 가 아니다: %d%s" % (len(GENESEQ), NL))
    raise SystemExit(1)
# whether the gene occupies 2122..2240 of the canonical unit is confirmed from the sequence rather than assumed
if CANSEQ[2121:2240] != GENESEQ:
    sys.stderr.write("정본 2122-2240 이 유전자와 다르다 — 중단%s" % NL)
    raise SystemExit(1)

SEQF = {}
for x in io.open(os.path.expanduser(os.environ["HSN"]), encoding="utf8"):
    f = x.rstrip(NL).split("\t")
    if len(f) >= 3 and f[0] != "sample":
        SEQF[(f[0], f[1])] = f[2]
RCF = {}
for i, x in enumerate([]):   # RCF is not used (method section 4)
    if i == 0:
        continue
    f = x.rstrip(NL).split("\t")
    if len(f) >= 6:
        RCF[(f[0], f[1])] = (f[5] == "Y")
FLIP = set()


def gene_hits(seq, tag):
    """(start, end) of the gene HSP and the dominant strand. Coordinates are 1-based, inclusive."""
    fa = os.path.join(W, "g_%s.fa" % tag)
    db = os.path.join(W, "gd_%s" % tag)
    io.open(fa, "w", newline=NL).write(">x" + NL + seq + NL)
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
                   capture_output=True)
    r = subprocess.run(["blastn", "-task", "blastn", "-query", GENE, "-db", db,
                        "-outfmt", "6 sstart send length pident qstart qend", "-dust", "no",
                        "-word_size", "9", "-evalue", "1e-5",
                        "-max_target_seqs", "100000000", "-num_threads", TH],
                       capture_output=True, text=True)
    plus = minus = 0
    hits = []
    for L in r.stdout.splitlines():
        g = L.split("\t")
        if len(g) < 6:
            continue
        a, b = int(g[0]), int(g[1])
        qs_, qe_ = int(g[4]), int(g[5])
        glen = len(GENESEQ)
        if a <= b:
            plus += 1
            s_, e_ = a - (qs_ - 1), b + (glen - qe_)
        else:
            minus += 1
            s_, e_ = b - (glen - qe_), a + (qs_ - 1)
        hits.append((max(1, s_), min(len(seq), e_), int(g[2]),
                     (qs_ - 1) + (glen - qe_)))   # fourth field = the amount of extension
    for e in (".ndb", ".nhr", ".nin", ".not", ".nsq", ".ntf", ".nto"):
        try:
            os.remove(db + e)
        except OSError:
            pass
    os.remove(fa)
    # Method F: interval overlap plus a per-fragment call, not 1 kb buckets.
    #   The flaw in bucketing: a phantom end at 145074 and the real end at 145147 differ by 73 bp,
    #                 but fall into buckets 144 and 145, so both survive and produce a spurious break and a unit whose boundary is 73 bp out.
    #   The extension rule applies as well, because a real end can fall one or two bases short where its last base carries a variant.
    keep = []
    for h in sorted(hits):
        if keep and h[0] <= keep[-1][1]:          # 1) where intervals overlap, keep the longer
            if h[2] > keep[-1][2]:
                keep[-1] = h
        else:
            keep.append(h)
    keep.sort(key=lambda x: x[1])
    i = 0
    while i + 1 < len(keep):                      # 2) for neighbours closer than LO, discard only the fragment
        if keep[i + 1][1] - keep[i][1] < LO:
            fa_, fb_ = keep[i][3] > FRAG, keep[i + 1][3] > FRAG
            if fa_ and not fb_:
                keep.pop(i)
                i = max(0, i - 1)
                continue
            if fb_ and not fa_:
                keep.pop(i + 1)
                i = max(0, i - 1)
                continue
        i += 1
    return [(a_, b_, c_) for (a_, b_, c_, d_) in keep], (minus > plus)


def call_variants(seqs, tag):
    q = os.path.join(W, "q_%s.fa" % tag)
    with io.open(q, "w", newline=NL) as g:
        for i, s in enumerate(seqs):
            g.write(">U%d%s%s%s" % (i, NL, s, NL))
    r = subprocess.run(["minimap2", "-x", PRESET, "-c", "--cs=long", "-t", TH,
                        "--secondary=no", CANREF, q], capture_output=True, text=True)
    OUT = {}
    for L in r.stdout.splitlines():
        f2 = L.split("\t")
        if len(f2) < 12:
            continue
        i = int(f2[0][1:])
        if i in OUT:
            continue
        cs = ""
        for x in f2[12:]:
            if x.startswith("cs:Z:"):
                cs = x[5:]
        qs, qe, ts, te = int(f2[2]), int(f2[3]), int(f2[7]), int(f2[8])
        toks = []
        reb = []
        tp = ts
        for m in CS.finditer(cs):
            op, v = m.group(1), m.group(2)
            if op == "=":
                reb.append(v.upper())
                tp += len(v)
            elif op == ":":
                n = int(v)
                reb.append(REFSEQ[tp:tp + n])
                tp += n
            elif op == "*":
                if CTX <= tp < CTX + 2240:
                    toks.append("%d%s" % (tp - CTX + 1, v[1].upper()))
                reb.append(v[1].upper())
                tp += 1
            elif op == "-":
                if CTX <= tp < CTX + 2240:
                    toks.append("%dD%d" % (tp - CTX + 1, len(v)))
                tp += len(v)
            elif op == "+":
                if CTX <= tp < CTX + 2240:
                    toks.append("%dI%d" % (tp - CTX + 1, len(v)))
                reb.append(v.upper())
        # No end filling: the alignment is made against a context-flanked reference with slack, so the alignment alone is used.
        OUT[i] = dict(toks=sorted(set(toks)), qs=qs, qe=qe, ts=ts, te=te,
                      reb="".join(reb), fill=0, cov_e=(2240 if te >= CTX + 2240 else 0), cov_s=(0 if ts <= CTX else 1))
    os.remove(q)
    return OUT


def verify(seq, res):
    if res is None:
        return "NOALN", "no alignment"
    if res["reb"] != seq[res["qs"]:res["qe"]]:
        return "SEQ_MISMATCH", "rebuilt %d vs original %d" % (len(res["reb"]), res["qe"] - res["qs"])
    if res["cov_e"] < 2240 or res["cov_s"] > 0:
        return "NO_FULL", "canonical %d..%d (+fill %d)" % (res["ts"], res["te"], res["fill"])
    return ("OK_FILLED" if res["fill"] else "OK"), (
        "end %d bp filled from sequence" % res["fill"] if res["fill"] else "")


HDR = ["sample", "slot", "N", "array_pos", "seg_type", "unit_start", "unit_end",
       "unit_len", "gap_to_next", "chain_ok", "n_var", "n_var_gene",
       "varlist", "verify", "verify_note", "slot_space"]
MAP = io.open(os.path.join(W, "CONSENSUS_MAP_G8_%s.tsv" % TAG), "w",
              encoding="utf8", newline=NL)
MAP.write("\t".join(HDR) + NL)
SUM = io.open(os.path.join(W, "CONSENSUS_SUM_G8_%s.tsv" % TAG), "w",
              encoding="utf8", newline=NL)
SUM.write("sample\tslot\tbp\tN_units\trevcomp\toff_prox\toff_dist\tn_break\n")
STAT = collections.Counter()

for line in io.open(sys.argv[1], encoding="utf8"):
    f = line.strip().split("\t")
    if len(f) < 2:
        continue
    smp, sl = f[0], f[1]
    fn = SEQF.get((smp, sl))
    if not fn:
        STAT["정본 없음"] += 1
        continue
    p = os.path.join(os.path.expanduser(os.environ["CMPD"]), fn)
    if not os.path.exists(p):
        STAT["파일 없음"] += 1
        sys.stderr.write("  파일 없음 %s_%s %s%s" % (smp, sl, p, NL))
        continue
    a = rd(p).upper()
    if RCF.get((smp, sl)):
        a = rc(a)
    hits, isrc = gene_hits(a, "%s_%s_%s" % (TAG, smp, sl))
    if isrc:
        # where the orientation is reversed the sequence is reverse-complemented and searched again, rather than the coordinates being flipped by hand
        a = rc(a)
        FLIP.add((smp, sl))
        hits, isrc = gene_hits(a, "%s_%s_%s_r" % (TAG, smp, sl))
        if isrc:
            STAT["뒤집은 뒤에도 - 가닥 우세"] += 1
    if len(hits) < 3:
        STAT["유전자 HSP 부족"] += 1
        continue
    # unit boundary = gene end to next gene end, which is the canonical phase
    ends = [e for (s, e, ln) in hits]
    bounds = []
    nbreak = 0
    for i in range(len(ends) - 1):
        g = ends[i + 1] - ends[i]
        if LO <= g <= HI:
            bounds.append((ends[i] + 1, ends[i + 1]))     # 1-based inclusive
        else:
            nbreak += 1
    if not bounds:
        STAT["사슬 없음"] += 1
        continue
    N = len(bounds)
    seqs = [a[max(0, s - 1 - PAD):min(len(a), e + PAD)] for (s, e) in bounds]   # slack on each side
    V = call_variants(seqs, "%s_%s_%s" % (TAG, smp, sl))
    off_prox = ends[0]                       # from the junction to the end of the first gene
    off_dist = len(a) - ends[-1]             # from the end of the last gene to the junction
    SUM.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\n"
              % (smp, sl, len(a), N, "Y" if (smp, sl) in FLIP else "N",
                 off_prox, off_dist, nbreak))
    # the proximal fragment, junction to first gene end; no variants are called
    MAP.write("\t".join([smp, sl, str(N), "0", "PARTIAL_PROX", "1", str(ends[0]),
                         str(ends[0]), "0", "-", "0", "0", "-", "SKIP", "not a complete unit", "assembly"]) + NL)
    STAT["PARTIAL_PROX"] += 1
    for i, (s, e) in enumerate(bounds):
        res = V.get(i)
        vd, note = verify(seqs[i], res)
        STAT[vd] += 1
        gap = (bounds[i + 1][0] - s) if i + 1 < N else 0
        if vd in ("OK", "OK_FILLED"):
            toks = res["toks"]
            ng = sum(1 for t in toks if int(NUM.match(t).group()) > 2121)
            STAT["유전자변이"] += ng
            vl = ",".join(toks) if toks else "-"
            MAP.write("\t".join([smp, sl, str(N), str(i + 1), "UNIT", str(s), str(e),
                                 str(e - s + 1), str(gap), "Y", str(len(toks)),
                                 str(ng), vl, vd, note, "assembly"]) + NL)
        else:
            MAP.write("\t".join([smp, sl, str(N), str(i + 1), "UNIT", str(s), str(e),
                                 str(e - s + 1), str(gap), "Y", "0", "0", "-",
                                 vd, note, "assembly"]) + NL)
    # the distal fragment, last gene end to junction
    if off_dist > 0:
        MAP.write("\t".join([smp, sl, str(N), str(N + 1), "PARTIAL_DIST",
                             str(ends[-1] + 1), str(len(a)), str(off_dist), "0", "-",
                             "0", "0", "-", "SKIP", "not a complete unit", "assembly"]) + NL)
        STAT["PARTIAL_DIST"] += 1
    sys.stderr.write("  %-10s 슬롯%s 유닛 %d · 끊김 %d%s" % (smp, sl, N, nbreak, NL))
    sys.stderr.flush()
MAP.close()
SUM.close()
print("consensus map G8")
for k, v in sorted(STAT.items()):
    print("  %-18s %d" % (k, v))
